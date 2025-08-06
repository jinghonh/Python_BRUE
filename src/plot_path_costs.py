import argparse
import os
import json
from dataclasses import dataclass
from typing import List, Tuple

import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat
from scipy.interpolate import interp1d


@dataclass
class PlotParams:
    """Configuration parameters for plotting."""
    save_path: str = 'results/'
    font_name: str = 'Arial'
    font_size: int = 10
    show_grid: bool = True
    figure_dpi: int = 300
    figure_size: Tuple[int, int] = (8, 6)


@dataclass
class Boundary:
    """Structure to hold boundary data."""
    left_x: np.ndarray
    left_y: np.ndarray
    right_x: np.ndarray
    right_y: np.ndarray


# 默认参数配置
DEFAULT_CONFIG = {
    'zeta': 15,              # 默认zeta值
    'subset_index': 0,       # 默认子集索引
    'num_flows': 50,         # 默认流量向量数量
    'cache_dir': 'matlab/cache',  # 默认缓存目录
    'results_dir': 'results'      # 默认结果目录
}


def calculate_travel_time(x: np.ndarray, t0: np.ndarray, c: np.ndarray) -> np.ndarray:
    """Calculate travel time using the BPR function."""
    return t0 * (1 + 0.15 * (x / c) ** 4)


def calculate_real_time(x: np.ndarray, m: np.ndarray, free_flow_time: np.ndarray, max_capacity: np.ndarray) -> np.ndarray:
    """Calculate real travel time including penalties."""
    index = np.where(np.sum(m, axis=0) != 0)[0]
    time = calculate_travel_time(x[:, index], free_flow_time[index], max_capacity[index])
    path_time = time @ m[:, index].T
    return path_time + 15 * (1 - np.exp(-0.02 * path_time))


def prepare_path_costs_data(total_valid_flow: np.ndarray, selected_indices: np.ndarray, relation_matrix: np.ndarray,
                            money_coeffs: np.ndarray, free_flow_time: np.ndarray, max_capacity: np.ndarray) -> Tuple:
    """Process and prepare path costs data for plotting."""
    q = len(selected_indices)
    money_costs_per_path = money_coeffs @ relation_matrix.T

    all_path_costs = []
    all_path_time_costs = []
    all_path_money_costs = []

    for i in range(q):
        current_flow = total_valid_flow[selected_indices[i], :]
        x = current_flow @ relation_matrix
        rt = calculate_real_time(np.atleast_2d(x), relation_matrix, free_flow_time, max_capacity)

        path_time_costs = rt.flatten()
        valid_paths = path_time_costs > 0
        
        path_money_costs_valid = money_costs_per_path[valid_paths]
        path_time_costs_valid = path_time_costs[valid_paths]

        all_path_time_costs.extend(path_time_costs_valid)
        all_path_money_costs.extend(path_money_costs_valid)

        costs = np.column_stack((path_time_costs_valid, path_money_costs_valid))
        sort_idx = np.argsort(path_money_costs_valid)
        all_path_costs.append(costs[sort_idx, :])

    color_variations = create_color_scheme(q)
    return all_path_costs, np.array(all_path_time_costs), np.array(all_path_money_costs), color_variations


def create_color_scheme(num_colors: int) -> np.ndarray:
    """Create a consistent color scheme for plots."""
    base_color = np.array([0.2, 0.4, 0.8])
    color_variations = np.zeros((num_colors, 3))
    for i in range(num_colors):
        color_shift = (i / (num_colors - 1 + np.finfo(float).eps)) * 0.4
        color_variations[i, :] = np.minimum(1, base_color + np.array([-0.15 + color_shift, color_shift, 0.1 - color_shift]))
    return color_variations


def calculate_feasible_region_boundary(all_path_time_costs: np.ndarray, all_path_money_costs: np.ndarray) -> Boundary:
    """Calculate the boundary of the feasible region based on all path costs."""
    if all_path_time_costs.size == 0 or all_path_money_costs.size == 0:
        return Boundary(np.array([]), np.array([]), np.array([]), np.array([]))
        
    unique_money_values = np.unique(all_path_money_costs)
    left_x, left_y, right_x, right_y = [], [], [], []

    for money_val in unique_money_values:
        indices = np.abs(all_path_money_costs - money_val) < 1e-3
        if np.any(indices):
            times_for_money = all_path_time_costs[indices]
            left_x.append(np.min(times_for_money))
            left_y.append(money_val)
            right_x.append(np.max(times_for_money))
            right_y.append(money_val)

    if not left_y:
        return Boundary(np.array([]), np.array([]), np.array([]), np.array([]))

    left_y, left_x = zip(*sorted(zip(left_y, left_x)))
    right_y, right_x = zip(*sorted(zip(right_y, right_x)))

    return Boundary(np.array(left_x), np.array(left_y), np.array(right_x), np.array(right_y))


def calculate_equilibrium_line(left_boundary_x: np.ndarray, upper_limit_x: np.ndarray, 
                               upper_limit_y: np.ndarray) -> np.ndarray:
    """
    Calculate equilibrium line position between left boundary and upper limit.
    
    Args:
        left_boundary_x: The left boundary X values (time costs)
        upper_limit_x: The upper limit X values (TTB_max_delta)
        upper_limit_y: The money costs corresponding to the boundaries
        
    Returns:
        The equilibrium line X values (T_eqm)
    """
    # Get normalized money costs
    money_max = np.max(upper_limit_y)
    money_min = np.min(upper_limit_y)
    money_range = money_max - money_min
    
    if money_range > 0:
        # Calculate weights: lower money costs get higher weights
        weights = 0.7 + 0.1 * ((upper_limit_y - money_min) / money_range)
    else:
        weights = 0.5 * np.ones_like(upper_limit_y)
    
    # Calculate equilibrium line position based on weights
    # With low money cost (large weight), eqm_limit_x will be closer to leftBoundaryX
    # With high money cost (small weight), eqm_limit_x will be closer to upperLimitX
    eqm_limit_x = left_boundary_x + (upper_limit_x - left_boundary_x) * weights
    
    return eqm_limit_x


def get_or_compute_TTB_max_delta_TTB_max(boundary: Boundary, zeta: int, results_dir: str = 'results') -> Tuple[np.ndarray, np.ndarray]:
    """
    获取或计算特定zeta值的TTB_max_delta和TTB_max，使用JSON格式存储以便手动编辑
    
    Args:
        boundary: 边界数据
        zeta: zeta值
        results_dir: 结果存储目录
        
    Returns:
        Tuple[np.ndarray, np.ndarray]: TTB_max_delta和TTB_max数组
    """
    # 创建存储目录
    os.makedirs(results_dir, exist_ok=True)
    cache_file = os.path.join(results_dir, f"tmax_teqm_zeta{zeta}.json")
    
    # 检查是否存在缓存文件
    if os.path.exists(cache_file):
        # print(f"Loading TTB_max_delta and TTB_max from {cache_file}")
        try:
            with open(cache_file, 'r') as f:
                data = json.load(f)
            
            # 从JSON加载数据
            money_values = np.array(data['money_values'])
            TTB_max_delta_values = np.array(data['TTB_max_delta'])
            TTB_max_values = np.array(data['TTB_max'])
            
            # 如果边界货币值与缓存不同，需要进行插值
            if len(money_values) != len(boundary.left_y) or not np.allclose(money_values, boundary.left_y):
                # 为TTB_max_delta创建插值函数
                f_tmax = interp1d(
                    money_values, TTB_max_delta_values,
                    kind='linear', bounds_error=False, fill_value="extrapolate"
                )
                TTB_max_delta = f_tmax(boundary.left_y)
                
                # 为TTB_max创建插值函数
                f_teqm = interp1d(
                    money_values, TTB_max_values,
                    kind='linear', bounds_error=False, fill_value="extrapolate"
                )
                TTB_max = f_teqm(boundary.left_y)
            else:
                # 直接使用缓存值
                TTB_max_delta = TTB_max_delta_values
                TTB_max = TTB_max_values
                
            # print(f"Successfully loaded TTB_max_delta and TTB_max for zeta={zeta}")
            return TTB_max_delta, TTB_max
        except Exception as e:
            print(f"Error loading TTB_max_delta and TTB_max from cache: {e}")
            print("Computing new values...")
    
    # 计算TTB_max_delta和TTB_max
    print(f"Computing new TTB_max_delta and TTB_max for zeta={zeta}")
    TTB_max_delta = boundary.left_x * 0.3 + boundary.right_x * 0.7
    TTB_max = calculate_equilibrium_line(boundary.left_x, TTB_max_delta, boundary.left_y)
    
    # 将数据保存为JSON格式（易于编辑）
    data = {
        'money_values': boundary.left_y.tolist(),
        'TTB_max_delta': TTB_max_delta.tolist(),
        'TTB_max': TTB_max.tolist(),
        'description': f"TTB_max_delta和T_eqm数据，用于zeta={zeta}的所有子集。可以手动编辑这些值以确保一致性。"
    }
    
    with open(cache_file, 'w') as f:
        json.dump(data, f, indent=2)
    
    print(f"Saved TTB_max_delta and TTB_max to {cache_file}")
    return TTB_max_delta, TTB_max


def configure_plot(ax: plt.Axes, params: PlotParams, title: str, xlabel: str, ylabel: str):
    """Apply common plot configurations."""
    # ax.set_title(title, fontsize=params.font_size + 4, fontweight='bold')
    ax.set_xlabel(xlabel, fontsize=params.font_size + 2, fontweight='bold')
    ax.set_ylabel(ylabel, fontsize=params.font_size + 2, fontweight='bold')
    if params.show_grid:
        ax.grid(True, which='major', alpha=0.2)
        ax.grid(True, which='minor', alpha=0.1)
    ax.tick_params(axis='both', which='major', labelsize=params.font_size)


def save_figure(fig: plt.Figure, base_filename: str, params: PlotParams, zeta: int, subset_index: int):
    """Save the figure with a descriptive name."""
    os.makedirs(params.save_path, exist_ok=True)
    filename = os.path.join(params.save_path, f"{base_filename}_zeta{zeta}_subset{subset_index}.pdf")
    fig.savefig(filename, dpi=params.figure_dpi, bbox_inches='tight')
    print(f"Figure saved to {filename}")


def plot_time_money_cost_relationship(all_path_costs: List[np.ndarray], boundary: Boundary,
                                      color_variations: np.ndarray, params: PlotParams, zeta: int, subset_index: int):
    """Plot the time-money cost relationship."""
    fig, ax = plt.subplots(figsize=params.figure_size)
    configure_plot(ax, params, 'Path Time-Money Cost Relationship', r'Time', r'Money Cost')

    if boundary.left_x.size > 0:
        boundary_x = np.concatenate([boundary.left_x, boundary.right_x[::-1]])
        boundary_y = np.concatenate([boundary.left_y, boundary.right_y[::-1]])
        ax.fill(boundary_x, boundary_y, color=[0.9, 0.95, 1], alpha=1, edgecolor='none', label='Feasible Region'); ax.plot(boundary.left_x, boundary.left_y, '-', color=[0.4, 0.5, 0.8], linewidth=0.5); ax.plot(boundary.right_x, boundary.right_y, '-', color=[0.4, 0.5, 0.8], linewidth=0.5)

    for i, costs in enumerate(all_path_costs):
        if costs.size > 0:
            ax.plot(costs[:, 0], costs[:, 1], '-', color=(*color_variations[i], 0.7), linewidth=1.2)
            ax.scatter(costs[:, 0], costs[:, 1], s=30, c=[color_variations[i]], marker='o', alpha=0.8,
                       edgecolors='none')

    handles, labels = ax.get_legend_handles_labels()
    if handles:
        unique_labels = dict(zip(labels, handles))
        ax.legend(unique_labels.values(), unique_labels.keys(), loc='best', fontsize=params.font_size - 1)
    
    save_figure(fig, 'path_time_money', params, zeta, subset_index)
    plt.close(fig)


def plot_teqm_flow_path_cost(all_path_costs: List[np.ndarray], color_variations: np.ndarray, boundary: Boundary, params: PlotParams, zeta: int, subset_index: int, relation_matrix: np.ndarray,
                               money_coeffs: np.ndarray, free_flow_time: np.ndarray, max_capacity: np.ndarray):
    """
    Plots the path cost line for the flow derived from the first T_eqm scatter point,
    with other sample flows also plotted.
    """
    # 1. Load scatter data for T_eqm point
    scatter_file = os.path.join('results', f'scatter_points_e{zeta}.json')
    try:
        with open(scatter_file, 'r') as f:
            scatter_data = json.load(f)
    except FileNotFoundError:
        print(f"Warning: Scatter file for T_eqm point not found: {scatter_file}")
        return

    # 2. Get the first T_eqm point
    teqm_points = scatter_data.get('tmax_constraint_points', [])
    if not teqm_points:
        print(f"Warning: No 'tmax_constraint_points' found in {scatter_file}")
        return
    
    f1, f2 = teqm_points[0]
    
    # 3. Create the flow vector for T_eqm
    f5 = 10000.0 - f1 - f2
    n_paths = relation_matrix.shape[0]
    teqm_flow = np.zeros(n_paths)
    teqm_flow[0] = f1  # Path 1
    teqm_flow[1] = f2  # Path 2

    if subset_index == 0 and n_paths >= 3:
        teqm_flow[2] = f5  # Path 5 is at index 2 for subset 0
    elif subset_index in [1, 2] and n_paths >= 5:
        teqm_flow[4] = f5  # Path 5 is at index 4 for subsets 1 and 2
            
    # 4. Calculate path costs for this single T_eqm flow
    x = teqm_flow @ relation_matrix
    rt = calculate_real_time(np.atleast_2d(x), relation_matrix, free_flow_time, max_capacity)
    path_time_costs = rt.flatten()
    money_costs_per_path = money_coeffs @ relation_matrix.T
    teqm_costs = np.column_stack((path_time_costs, money_costs_per_path))
    
    # 5. Plotting
    fig, ax = plt.subplots(figsize=params.figure_size)
    configure_plot(ax, params, f'Path Costs with Highlighted TTB_max Flow (Zeta={zeta})', r'Time', r'Money Cost')
    
    # Layer 1: Plot background feasible region
    if boundary.left_x.size > 0:
        boundary_x = np.concatenate([boundary.left_x, boundary.right_x[::-1]])
        boundary_y = np.concatenate([boundary.left_y, boundary.right_y[::-1]])
        ax.fill(boundary_x, boundary_y, color=[0.9, 0.95, 1], alpha=1, edgecolor='none', label='Feasible Region'); ax.plot(boundary.left_x, boundary.left_y, '-', color=[0.4, 0.5, 0.8], linewidth=0.5); ax.plot(boundary.right_x, boundary.right_y, '-', color=[0.4, 0.5, 0.8], linewidth=0.5)

    # Layer 2: Plot all other flows
    for i, costs in enumerate(all_path_costs):
        if costs.size > 0:
            ax.plot(costs[:, 0], costs[:, 1], '-', color=(*color_variations[i], 0.7), linewidth=1.2)
            ax.scatter(costs[:, 0], costs[:, 1], s=30, c=[color_variations[i]], marker='o', alpha=0.8,
                       edgecolors='none')

    # Layer 3: Plot the prominent TTB_max cost line
    sort_idx = np.argsort(teqm_costs[:, 1])
    ax.plot(teqm_costs[sort_idx, 0], teqm_costs[sort_idx, 1], '-', color='#984ea3', linewidth=2.5, marker='D', markersize=8, label=r'Flow from $TTB_{max}$ point', zorder=10)
    
    # 6. Add legend and save
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        unique_labels = dict(zip(labels, handles))
        ax.legend(unique_labels.values(), unique_labels.keys(), loc='best', fontsize=params.font_size - 1)

    save_figure(fig, 'teqm_flow_path_cost_with_others', params, zeta, subset_index)
    plt.close(fig)


def plot_path_costs_with_upper_limit(all_path_costs: List[np.ndarray], boundary: Boundary,
                                     color_variations: np.ndarray, params: PlotParams, zeta: int, subset_index: int):
    """Plot path costs with an upper limit and equilibrium line, with correct layering."""
    fig, ax = plt.subplots(figsize=params.figure_size)
    configure_plot(ax, params, 'Path Costs with Upper Limit', r'Time', r'Money Cost')

    # --- Data Processing ---
    if boundary.left_x.size == 0:
        print("Warning: Boundary data is empty. Skipping plot generation.")
        plt.close(fig)
        return

    # 使用共享函数获取TTB_max_delta和TTB_max
    TTB_max_delta, TTB_max = get_or_compute_TTB_max_delta_TTB_max(boundary, zeta, params.save_path)
    
    f_tmax = interp1d(boundary.left_y, TTB_max_delta, kind='nearest', bounds_error=False, fill_value="extrapolate")

    feasible_paths_data = []
    infeasible_paths_data = []
    feasible_costs_points = []

    for i, costs in enumerate(all_path_costs):
        if costs.size > 0:
            is_feasible = np.all(costs[:, 0] <= f_tmax(costs[:, 1]))
            if is_feasible:
                feasible_paths_data.append({'costs': costs, 'color': color_variations[i]})
                feasible_costs_points.append(costs)
            else:
                infeasible_paths_data.append({'costs': costs, 'color': [0.8, 0.8, 0.8]})

    # --- Drawing Layers ---

    # Layer 1: Total Feasible Region (Background)
    boundary_x = np.concatenate([boundary.left_x, boundary.right_x[::-1]])
    boundary_y = np.concatenate([boundary.left_y, boundary.right_y[::-1]])
    ax.fill(boundary_x, boundary_y, color=[0.9, 0.95, 1], alpha=1, edgecolor='none', label='Total Feasible Region')

    # Layer 2: Infeasible Paths
    infeasible_handle = None
    for path in infeasible_paths_data:
        costs = path['costs']
        color = path['color']
        ax.plot(costs[:, 0], costs[:, 1], '-', color=(*color[:3], 0.5), linewidth=1.0, zorder=2)
        h = ax.scatter(costs[:, 0], costs[:, 1], s=25, c=[color], marker='o', alpha=0.3, edgecolors='none', zorder=2)
        if infeasible_handle is None:
            infeasible_handle = h

    # Layer 3: Feasible Flow Sub-Region
    feasible_region_handle = None
    if feasible_costs_points:
        all_feasible_time = np.concatenate([c[:, 0] for c in feasible_costs_points])
        all_feasible_money = np.concatenate([c[:, 1] for c in feasible_costs_points])
        if all_feasible_time.size > 2:
            feasible_sub_boundary = calculate_feasible_region_boundary(all_feasible_time, all_feasible_money)
            if feasible_sub_boundary.left_x.size > 0:
                sub_boundary_x = np.concatenate([feasible_sub_boundary.left_x, feasible_sub_boundary.right_x[::-1]])
                sub_boundary_y = np.concatenate([feasible_sub_boundary.left_y, feasible_sub_boundary.right_y[::-1]])
                feasible_region_handle = ax.fill(sub_boundary_x, sub_boundary_y, color=[0.7, 0.95, 0.7], alpha=0.5,
                        edgecolor=[0.4, 0.6, 0.4], linewidth=0.5, label='Feasible Flow Region', zorder=3)[0]

    # Layer 4: Feasible Paths
    feasible_handle = None
    for path in feasible_paths_data:
        costs = path['costs']
        color = path['color']
        ax.plot(costs[:, 0], costs[:, 1], '-', color=(*color[:3], 0.7), linewidth=1.2, zorder=4)
        h = ax.scatter(costs[:, 0], costs[:, 1], s=30, c=[color], marker='o', alpha=0.8, edgecolors='none', zorder=4)
        if feasible_handle is None:
            feasible_handle = h

    # Layer 5: TTB_max_delta and T_eqm lines (Topmost)
    TTB_max_delta_handle = ax.plot(TTB_max_delta, boundary.left_y, '-', color=[0.8, 0.2, 0.2], linewidth=2.5, label=r'$TTB_{max} + \delta_k$', zorder=10)[0]
    TTB_max_handle = ax.plot(TTB_max, boundary.left_y, '-.', color=[0.5, 0.0, 0.8], linewidth=2.5, label=r'$TTB_{max}$', zorder=10)[0]

    # 添加图例
    legend_handles = []
    legend_labels = []
    
    # 总可行区域
    legend_handles.append(ax.fill([], [], color=[0.9, 0.95, 1], alpha=1, edgecolor='none')[0])
    legend_labels.append('Total Feasible Region')
    
    # 可行流量区域
    if feasible_region_handle is not None:
        legend_handles.append(feasible_region_handle)
        legend_labels.append('Feasible Flow Region')
    
    # 可行路径点
    if feasible_handle is not None:
        legend_handles.append(feasible_handle)
        legend_labels.append('Feasible Flow Paths')
    
    # 不可行路径点
    if infeasible_handle is not None:
        legend_handles.append(infeasible_handle)
        legend_labels.append('Infeasible Flow Paths')
    
    # TTB_max_delta线
    legend_handles.append(TTB_max_delta_handle)
    legend_labels.append(r'$TTB_{max} + \delta_k$')
    
    # T_eqm线
    legend_handles.append(TTB_max_handle)
    legend_labels.append(r'$TTB_{max}$')
    
    # 显示图例
    ax.legend(
        legend_handles,
        legend_labels,
        loc='best',
        fontsize=params.font_size - 1,
        framealpha=0.9,
        fancybox=True,
        shadow=True
    )

    save_figure(fig, 'path_costs_upper_limit', params, zeta, subset_index)
    plt.close(fig)


def plot_path_costs_below_equilibrium(all_path_costs: List[np.ndarray], boundary: Boundary,
                                      color_variations: np.ndarray, params: PlotParams, zeta: int, subset_index: int):
    """筛选并绘制T_eqm（紫色折线）下方的可行方案。"""
    fig, ax = plt.subplots(figsize=params.figure_size)
    configure_plot(ax, params, 'Path Costs Below Equilibrium Line', r'Time', r'Money Cost')

    # --- Data Processing ---
    if boundary.left_x.size == 0:
        print("Warning: Boundary data is empty. Skipping plot generation.")
        plt.close(fig)
        return

    # 使用共享函数获取TTB_max_delta和TTB_max
    TTB_max_delta, TTB_max = get_or_compute_TTB_max_delta_TTB_max(boundary, zeta, params.save_path)
    
    # 创建T_eqm的插值函数用于判断
    f_teqm = interp1d(boundary.left_y, TTB_max, kind='nearest', bounds_error=False, fill_value="extrapolate")

    # 筛选出所有点都在T_eqm下方的路径
    below_eqm_paths_data = []
    below_eqm_costs_points = []
    other_paths_data = []

    for i, costs in enumerate(all_path_costs):
        if costs.size > 0:
            # 判断所有点是否都在T_eqm下方
            is_below_eqm = np.all(costs[:, 0] <= f_teqm(costs[:, 1]))
            if is_below_eqm:
                # 使用更浅的紫色系颜色
                purple_color = np.array([0.7, 0.5, 0.9]) + np.array([0.05, 0.05, 0.0]) * (i % 3)
                purple_color = np.minimum(purple_color, 1.0)  # 确保颜色值不超过1
                below_eqm_paths_data.append({'costs': costs, 'color': purple_color})
                below_eqm_costs_points.append(costs)
            else:
                # 其他路径用浅灰色
                other_paths_data.append({'costs': costs, 'color': [0.85, 0.85, 0.85]})

    # --- Drawing Layers ---

    # Layer 1: Total Feasible Region (Background)
    boundary_x = np.concatenate([boundary.left_x, boundary.right_x[::-1]])
    boundary_y = np.concatenate([boundary.left_y, boundary.right_y[::-1]])
    ax.fill(boundary_x, boundary_y, color=[0.95, 0.95, 0.95], alpha=1, edgecolor='none', label='Total Feasible Region')

    # Layer 2: Other Paths (非平衡线下方的路径)
    for path in other_paths_data:
        costs = path['costs']
        color = path['color']
        ax.plot(costs[:, 0], costs[:, 1], '-', color=(*color[:3], 0.3), linewidth=0.8, zorder=2)
        ax.scatter(costs[:, 0], costs[:, 1], s=20, c=[color], marker='o', alpha=0.2, edgecolors='none', zorder=2)

    # Layer 3: Below Equilibrium Region
    if below_eqm_costs_points:
        all_below_eqm_time = np.concatenate([c[:, 0] for c in below_eqm_costs_points])
        all_below_eqm_money = np.concatenate([c[:, 1] for c in below_eqm_costs_points])
        if all_below_eqm_time.size > 2:
            below_eqm_boundary = calculate_feasible_region_boundary(all_below_eqm_time, all_below_eqm_money)
            if below_eqm_boundary.left_x.size > 0:
                below_eqm_x = np.concatenate([below_eqm_boundary.left_x, below_eqm_boundary.right_x[::-1]])
                below_eqm_y = np.concatenate([below_eqm_boundary.left_y, below_eqm_boundary.right_y[::-1]])
                ax.fill(below_eqm_x, below_eqm_y, color=[0.9, 0.8, 1.0], alpha=0.5,
                        edgecolor=[0.6, 0.4, 0.7], linewidth=0.5, label='Below Equilibrium Region', zorder=3)

    # Layer 4: Below Equilibrium Paths
    for path in below_eqm_paths_data:
        costs = path['costs']
        color = path['color']
        ax.plot(costs[:, 0], costs[:, 1], '-', color=(*color[:3], 0.5), linewidth=1.5, zorder=4)
        ax.scatter(costs[:, 0], costs[:, 1], s=35, c=[color], marker='o', alpha=0.6, edgecolors='none', zorder=4)

    # Layer 5: T_eqm line (Topmost)
    ax.plot(TTB_max, boundary.left_y, '-.', color=[0.5, 0.0, 0.8], linewidth=2.5, label=r'$TTB_{max}$', zorder=10)

    # 添加图例
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        unique_labels = dict(zip(labels, handles))
        ax.legend(unique_labels.values(), unique_labels.keys(), loc='best', fontsize=params.font_size - 1)

    save_figure(fig, 'path_costs_below_equilibrium', params, zeta, subset_index)
    plt.close(fig)


def plot_rbs_point_cost(all_path_costs: List[np.ndarray], boundary: Boundary,
                       color_variations: np.ndarray, params: PlotParams, zeta: int, subset_index: int,
                       relation_matrix: np.ndarray, money_coeffs: np.ndarray, free_flow_time: np.ndarray,
                       max_capacity: np.ndarray):
    """
    在plot_path_costs_with_upper_limit的基础上，额外绘制绿色区域(RBS)中一个点的双目标成本折线
    
    Args:
        all_path_costs: 所有路径成本数据
        boundary: 边界数据
        color_variations: 颜色方案
        params: 绘图参数
        zeta: zeta值
        subset_index: 子集索引
        relation_matrix: 关系矩阵
        money_coeffs: 金钱成本系数
        free_flow_time: 自由流时间
        max_capacity: 最大容量
    """
    # --- Data Processing ---
    if boundary.left_x.size == 0:
        print("Warning: Boundary data is empty. Skipping plot generation.")
        return
    
    # 使用共享函数获取TTB_max_delta和TTB_max
    TTB_max_delta, TTB_max = get_or_compute_TTB_max_delta_TTB_max(boundary, zeta, params.save_path)
    
    f_tmax = interp1d(boundary.left_y, TTB_max_delta, kind='nearest', bounds_error=False, fill_value="extrapolate")

    feasible_paths_data = []
    infeasible_paths_data = []
    feasible_costs_points = []

    for i, costs in enumerate(all_path_costs):
        if costs.size > 0:
            is_feasible = np.all(costs[:, 0] <= f_tmax(costs[:, 1]))
            if is_feasible:
                feasible_paths_data.append({'costs': costs, 'color': color_variations[i]})
                feasible_costs_points.append(costs)
            else:
                infeasible_paths_data.append({'costs': costs, 'color': [0.8, 0.8, 0.8]})
    
    # 1. 加载散点数据获取RBS点
    scatter_file = os.path.join('results', f'scatter_points_e{zeta}.json')
    try:
        with open(scatter_file, 'r') as f:
            scatter_data = json.load(f)
    except FileNotFoundError:
        print(f"警告: RBS点的散点文件未找到: {scatter_file}")
        return
    
    # 2. 获取RBS区域的第一个点
    rbs_points = scatter_data.get('path_constraint_points', [])
    if not rbs_points:
        print(f"警告: 未在 {scatter_file} 中找到'path_constraint_points'")
        return
    
    f1, f2 = rbs_points[0]
    
    # 3. 创建流量向量
    f5 = 10000.0 - f1 - f2
    n_paths = relation_matrix.shape[0]
    rbs_flow = np.zeros(n_paths)
    rbs_flow[0] = f1  # Path 1
    rbs_flow[1] = f2  # Path 2
    
    if subset_index == 0 and n_paths >= 3:
        rbs_flow[2] = f5  # Path 5 is at index 2 for subset 0
    elif subset_index in [1, 2] and n_paths >= 5:
        rbs_flow[4] = f5  # Path 5 is at index 4 for subsets 1 and 2
            
    # 4. 计算路径成本
    x = rbs_flow @ relation_matrix
    rt = calculate_real_time(np.atleast_2d(x), relation_matrix, free_flow_time, max_capacity)
    path_time_costs = rt.flatten()
    money_costs_per_path = money_coeffs @ relation_matrix.T
    rbs_costs = np.column_stack((path_time_costs, money_costs_per_path))
    
    # --- Drawing Layers ---
    fig, ax = plt.subplots(figsize=params.figure_size)
    configure_plot(ax, params, f'Path Costs with Upper Limit and RBS Point (Zeta={zeta})', r'Time', r'Money Cost')

    # Layer 1: Total Feasible Region (Background)
    boundary_x = np.concatenate([boundary.left_x, boundary.right_x[::-1]])
    boundary_y = np.concatenate([boundary.left_y, boundary.right_y[::-1]])
    ax.fill(boundary_x, boundary_y, color=[0.9, 0.95, 1], alpha=1, edgecolor='none', label='Total Feasible Region')

    # Layer 2: Infeasible Paths
    infeasible_handle = None
    for path in infeasible_paths_data:
        costs = path['costs']
        color = path['color']
        ax.plot(costs[:, 0], costs[:, 1], '-', color=(*color[:3], 0.5), linewidth=1.0, zorder=2)
        h = ax.scatter(costs[:, 0], costs[:, 1], s=25, c=[color], marker='o', alpha=0.3, edgecolors='none', zorder=2)
        if infeasible_handle is None:
            infeasible_handle = h

    # Layer 3: Feasible Flow Sub-Region
    feasible_region_handle = None
    if feasible_costs_points:
        all_feasible_time = np.concatenate([c[:, 0] for c in feasible_costs_points])
        all_feasible_money = np.concatenate([c[:, 1] for c in feasible_costs_points])
        if all_feasible_time.size > 2:
            feasible_sub_boundary = calculate_feasible_region_boundary(all_feasible_time, all_feasible_money)
            if feasible_sub_boundary.left_x.size > 0:
                sub_boundary_x = np.concatenate([feasible_sub_boundary.left_x, feasible_sub_boundary.right_x[::-1]])
                sub_boundary_y = np.concatenate([feasible_sub_boundary.left_y, feasible_sub_boundary.right_y[::-1]])
                feasible_region_handle = ax.fill(sub_boundary_x, sub_boundary_y, color=[0.7, 0.95, 0.7], alpha=0.5,
                        edgecolor=[0.4, 0.6, 0.4], linewidth=0.5, label='Feasible Flow Region', zorder=3)[0]

    # Layer 4: Feasible Paths
    feasible_handle = None
    for path in feasible_paths_data:
        costs = path['costs']
        color = path['color']
        ax.plot(costs[:, 0], costs[:, 1], '-', color=(*color[:3], 0.7), linewidth=1.2, zorder=4)
        h = ax.scatter(costs[:, 0], costs[:, 1], s=30, c=[color], marker='o', alpha=0.8, edgecolors='none', zorder=4)
        if feasible_handle is None:
            feasible_handle = h

    # Layer 5: TTB_max_delta and T_eqm lines
    TTB_max_delta_handle = ax.plot(TTB_max_delta, boundary.left_y, '-', color=[0.8, 0.2, 0.2], linewidth=2.5, label=r'$TTB_{max} + \delta_k$', zorder=10)[0]
    TTB_max_handle = ax.plot(TTB_max, boundary.left_y, '-.', color=[0.5, 0.0, 0.8], linewidth=2.5, label=r'$TTB_{max}$', zorder=10)[0]
    
    # Layer 6: RBS点的成本线 (最上层)
    sort_idx = np.argsort(rbs_costs[:, 1])
    rbs_handle = ax.plot(rbs_costs[sort_idx, 0], rbs_costs[sort_idx, 1], '-', color='yellow', linewidth=2.5, marker='D', markersize=8, label=r'Flow from $RBS_0^\varepsilon$ point', zorder=11)[0]

    # 添加图例
    legend_handles = []
    legend_labels = []
    
    # 总可行区域
    legend_handles.append(ax.fill([], [], color=[0.9, 0.95, 1], alpha=1, edgecolor='none')[0])
    legend_labels.append('Total Feasible Region')
    
    # 可行流量区域
    if feasible_region_handle is not None:
        legend_handles.append(feasible_region_handle)
        legend_labels.append('Feasible Flow Region')
    
    # 可行路径点
    if feasible_handle is not None:
        legend_handles.append(feasible_handle)
        legend_labels.append('Feasible Flow Paths')
    
    # 不可行路径点
    if infeasible_handle is not None:
        legend_handles.append(infeasible_handle)
        legend_labels.append('Infeasible Flow Paths')
    
    # TTB_max_delta线
    legend_handles.append(TTB_max_delta_handle)
    legend_labels.append(r'$TTB_{max} + \delta_k$')
    
    # T_eqm线
    legend_handles.append(TTB_max_handle)
    legend_labels.append(r'$TTB_{max}$')
    
    # RBS点的成本线
    legend_handles.append(rbs_handle)
    legend_labels.append(r'Flow from $RBS_0^\varepsilon$ point')
    
    # 显示图例
    ax.legend(
        legend_handles,
        legend_labels,
        loc='best',
        fontsize=params.font_size - 1,
        framealpha=0.9,
        fancybox=True,
        shadow=True
    )

    save_figure(fig, 'rbs_point_path_cost_with_limits', params, zeta, subset_index)
    plt.close(fig)


def plot_ttbmax_point_cost(all_path_costs: List[np.ndarray], boundary: Boundary,
                         color_variations: np.ndarray, params: PlotParams, zeta: int, subset_index: int,
                         relation_matrix: np.ndarray, money_coeffs: np.ndarray, free_flow_time: np.ndarray,
                         max_capacity: np.ndarray):
    """
    在plot_path_costs_with_upper_limit的基础上，额外绘制紫色区域(TTB_max)中一个点的双目标成本折线
    
    Args:
        all_path_costs: 所有路径成本数据
        boundary: 边界数据
        color_variations: 颜色方案
        params: 绘图参数
        zeta: zeta值
        subset_index: 子集索引
        relation_matrix: 关系矩阵
        money_coeffs: 金钱成本系数
        free_flow_time: 自由流时间
        max_capacity: 最大容量
    """
    # --- Data Processing ---
    if boundary.left_x.size == 0:
        print("Warning: Boundary data is empty. Skipping plot generation.")
        return
    
    # 使用共享函数获取TTB_max_delta和TTB_max
    TTB_max_delta, TTB_max = get_or_compute_TTB_max_delta_TTB_max(boundary, zeta, params.save_path)
    
    f_tmax = interp1d(boundary.left_y, TTB_max_delta, kind='nearest', bounds_error=False, fill_value="extrapolate")

    feasible_paths_data = []
    infeasible_paths_data = []
    feasible_costs_points = []

    for i, costs in enumerate(all_path_costs):
        if costs.size > 0:
            is_feasible = np.all(costs[:, 0] <= f_tmax(costs[:, 1]))
            if is_feasible:
                feasible_paths_data.append({'costs': costs, 'color': color_variations[i]})
                feasible_costs_points.append(costs)
            else:
                infeasible_paths_data.append({'costs': costs, 'color': [0.8, 0.8, 0.8]})
    
    # 1. 加载散点数据获取TTB_max点
    scatter_file = os.path.join('results', f'scatter_points_e{zeta}.json')
    try:
        with open(scatter_file, 'r') as f:
            scatter_data = json.load(f)
    except FileNotFoundError:
        print(f"警告: TTB_max点的散点文件未找到: {scatter_file}")
        return
    
    # 2. 获取TTB_max区域的第一个点
    ttbmax_points = scatter_data.get('tmax_constraint_points', [])
    if not ttbmax_points:
        print(f"警告: 未在 {scatter_file} 中找到'tmax_constraint_points'")
        return
    
    f1, f2 = ttbmax_points[0]
    
    # 3. 创建流量向量
    f5 = 10000.0 - f1 - f2
    n_paths = relation_matrix.shape[0]
    ttbmax_flow = np.zeros(n_paths)
    ttbmax_flow[0] = f1  # Path 1
    ttbmax_flow[1] = f2  # Path 2
    
    if subset_index == 0 and n_paths >= 3:
        ttbmax_flow[2] = f5  # Path 5 is at index 2 for subset 0
    elif subset_index in [1, 2] and n_paths >= 5:
        ttbmax_flow[4] = f5  # Path 5 is at index 4 for subsets 1 and 2
            
    # 4. 计算路径成本
    x = ttbmax_flow @ relation_matrix
    rt = calculate_real_time(np.atleast_2d(x), relation_matrix, free_flow_time, max_capacity)
    path_time_costs = rt.flatten()
    money_costs_per_path = money_coeffs @ relation_matrix.T
    ttbmax_costs = np.column_stack((path_time_costs, money_costs_per_path))
    
    # --- Drawing Layers ---
    fig, ax = plt.subplots(figsize=params.figure_size)
    configure_plot(ax, params, f'Path Costs with Upper Limit and TTB_max Point (Zeta={zeta})', r'Time', r'Money Cost')

    # Layer 1: Total Feasible Region (Background)
    boundary_x = np.concatenate([boundary.left_x, boundary.right_x[::-1]])
    boundary_y = np.concatenate([boundary.left_y, boundary.right_y[::-1]])
    ax.fill(boundary_x, boundary_y, color=[0.9, 0.95, 1], alpha=1, edgecolor='none', label='Total Feasible Region')

    # Layer 2: Infeasible Paths
    infeasible_handle = None
    for path in infeasible_paths_data:
        costs = path['costs']
        color = path['color']
        ax.plot(costs[:, 0], costs[:, 1], '-', color=(*color[:3], 0.5), linewidth=1.0, zorder=2)
        h = ax.scatter(costs[:, 0], costs[:, 1], s=25, c=[color], marker='o', alpha=0.3, edgecolors='none', zorder=2)
        if infeasible_handle is None:
            infeasible_handle = h

    # Layer 3: Feasible Flow Sub-Region
    feasible_region_handle = None
    if feasible_costs_points:
        all_feasible_time = np.concatenate([c[:, 0] for c in feasible_costs_points])
        all_feasible_money = np.concatenate([c[:, 1] for c in feasible_costs_points])
        if all_feasible_time.size > 2:
            feasible_sub_boundary = calculate_feasible_region_boundary(all_feasible_time, all_feasible_money)
            if feasible_sub_boundary.left_x.size > 0:
                sub_boundary_x = np.concatenate([feasible_sub_boundary.left_x, feasible_sub_boundary.right_x[::-1]])
                sub_boundary_y = np.concatenate([feasible_sub_boundary.left_y, feasible_sub_boundary.right_y[::-1]])
                feasible_region_handle = ax.fill(sub_boundary_x, sub_boundary_y, color=[0.7, 0.95, 0.7], alpha=0.5,
                        edgecolor=[0.4, 0.6, 0.4], linewidth=0.5, label='Feasible Flow Region', zorder=3)[0]

    # Layer 4: Feasible Paths
    feasible_handle = None
    for path in feasible_paths_data:
        costs = path['costs']
        color = path['color']
        ax.plot(costs[:, 0], costs[:, 1], '-', color=(*color[:3], 0.7), linewidth=1.2, zorder=4)
        h = ax.scatter(costs[:, 0], costs[:, 1], s=30, c=[color], marker='o', alpha=0.8, edgecolors='none', zorder=4)
        if feasible_handle is None:
            feasible_handle = h

    # Layer 5: TTB_max_delta and T_eqm lines
    TTB_max_delta_handle = ax.plot(TTB_max_delta, boundary.left_y, '-', color=[0.8, 0.2, 0.2], linewidth=2.5, label=r'$TTB_{max} + \delta_k$', zorder=10)[0]
    TTB_max_handle = ax.plot(TTB_max, boundary.left_y, '-.', color=[0.5, 0.0, 0.8], linewidth=2.5, label=r'$TTB_{max}$', zorder=10)[0]
    
    # Layer 6: TTB_max点的成本线 (最上层)
    sort_idx = np.argsort(ttbmax_costs[:, 1])
    ttbmax_handle = ax.plot(ttbmax_costs[sort_idx, 0], ttbmax_costs[sort_idx, 1], '-', color='orange', linewidth=2.5, marker='D', markersize=8, label=r'Flow from $TTB_{max}$ point', zorder=11)[0]

    # 添加图例
    legend_handles = []
    legend_labels = []
    
    # 总可行区域
    legend_handles.append(ax.fill([], [], color=[0.9, 0.95, 1], alpha=1, edgecolor='none')[0])
    legend_labels.append('Total Feasible Region')
    
    # 可行流量区域
    if feasible_region_handle is not None:
        legend_handles.append(feasible_region_handle)
        legend_labels.append('Feasible Flow Region')
    
    # 可行路径点
    if feasible_handle is not None:
        legend_handles.append(feasible_handle)
        legend_labels.append('Feasible Flow Paths')
    
    # 不可行路径点
    if infeasible_handle is not None:
        legend_handles.append(infeasible_handle)
        legend_labels.append('Infeasible Flow Paths')
    
    # TTB_max_delta线
    legend_handles.append(TTB_max_delta_handle)
    legend_labels.append(r'$TTB_{max} + \delta_k$')
    
    # T_eqm线
    legend_handles.append(TTB_max_handle)
    legend_labels.append(r'$TTB_{max}$')
    
    # TTB_max点的成本线
    legend_handles.append(ttbmax_handle)
    legend_labels.append(r'Flow from $TTB_{max}$ point')
    
    # 显示图例
    ax.legend(
        legend_handles,
        legend_labels,
        loc='best',
        fontsize=params.font_size - 1,
        framealpha=0.9,
        fancybox=True,
        shadow=True
    )

    save_figure(fig, 'ttbmax_point_path_cost_with_limits', params, zeta, subset_index)
    plt.close(fig)


def run_with_params(zeta=None, subset_index=None, num_flows=None, cache_dir=None, results_dir=None, plot_params=None):
    """
    使用内部参数运行路径成本绘图功能
    
    参数:
        zeta (int): zeta值，默认为None (使用DEFAULT_CONFIG中的值)
        subset_index (int): 子集索引，默认为None (使用DEFAULT_CONFIG中的值)
        num_flows (int): 流量向量数量，默认为None (使用DEFAULT_CONFIG中的值)
        cache_dir (str): 缓存目录，默认为None (使用DEFAULT_CONFIG中的值)
        results_dir (str): 结果目录，默认为None (使用DEFAULT_CONFIG中的值)
        plot_params (PlotParams): 绘图参数，默认为None (使用默认PlotParams)
    
    返回:
        dict: 包含边界和路径成本数据的字典
    """
    # 使用默认配置或传入的参数
    zeta = zeta if zeta is not None else DEFAULT_CONFIG['zeta']
    subset_index = subset_index if subset_index is not None else DEFAULT_CONFIG['subset_index']
    num_flows = num_flows if num_flows is not None else DEFAULT_CONFIG['num_flows']
    cache_dir = cache_dir if cache_dir is not None else DEFAULT_CONFIG['cache_dir']
    results_dir = results_dir if results_dir is not None else DEFAULT_CONFIG['results_dir']
    
    # 设置绘图样式
    plt.rcParams['text.usetex'] = True
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['font.serif'] = ['Computer Modern Roman', 'Times New Roman']
    plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'
    
    # 设置绘图参数
    params = plot_params if plot_params is not None else PlotParams(save_path=results_dir)
    
    # 加载数据
    cache_file = os.path.join(cache_dir, f'cache_zeta{zeta}_subset{subset_index}.mat')
    if not os.path.exists(cache_file):
        print(f"Error: Cache file not found at {cache_file}")
        return None
    
    # print(f"Loading data from {cache_file}...")
    data = loadmat(cache_file)
    total_valid_flow = data.get('totalValidFlow')
    relation_matrix = data.get('relationMatrix')
    
    if total_valid_flow is None or relation_matrix is None:
        print("Error: 'totalValidFlow' or 'relationMatrix' not found in cache file.")
        return None
    
    print("Data loaded successfully.")
    
    # 选择随机流量索引
    num_total_flows = total_valid_flow.shape[0]
    if num_total_flows == 0:
        print("No valid flow data found in the cache file.")
        return None
    
    selected_indices = np.random.permutation(num_total_flows)[:min(num_flows, num_total_flows)]
    
    # 定义常量
    money_coeffs = np.array([20, 15, 1, 0, 0, 0, 0, 1])
    free_flow_time = np.array([18, 22.5, 12, 24, 2.4, 6, 24, 12])
    max_capacity = np.array([3600, 3600, 1800, 1800, 1800, 1800, 1800, 1800])
    
    # 准备数据
    all_path_costs, all_path_time_costs, all_path_money_costs, color_variations = \
        prepare_path_costs_data(total_valid_flow, selected_indices, relation_matrix,
                              money_coeffs, free_flow_time, max_capacity)
    
    # 计算总体边界
    boundary = calculate_feasible_region_boundary(all_path_time_costs, all_path_money_costs)
    
    # 生成图表
    print("Generating plots...")
    plot_time_money_cost_relationship(all_path_costs, boundary, color_variations, params, zeta, subset_index)

    
    # Call the new functions to plot special points' path costs
    if subset_index == 0:
        plot_teqm_flow_path_cost(all_path_costs, color_variations, boundary, params, zeta, subset_index, relation_matrix,
                                money_coeffs, free_flow_time, max_capacity)
        # 新增：绘制RBS点和TTB_max点的双目标成本折线图
        plot_rbs_point_cost(all_path_costs, boundary, color_variations, params, zeta, subset_index, relation_matrix,
                          money_coeffs, free_flow_time, max_capacity)
        plot_ttbmax_point_cost(all_path_costs, boundary, color_variations, params, zeta, subset_index, relation_matrix,
                             money_coeffs, free_flow_time, max_capacity)

    plot_path_costs_with_upper_limit(all_path_costs, boundary, color_variations, params, zeta, subset_index)
    plot_path_costs_below_equilibrium(all_path_costs, boundary, color_variations, params, zeta, subset_index)
    print("Plots generated successfully.")
    
    # 返回数据以便进一步处理
    return {
        'boundary': boundary,
        'all_path_costs': all_path_costs,
        'all_path_time_costs': all_path_time_costs,
        'all_path_money_costs': all_path_money_costs
    }


def main():
    """Main function to load data and generate plots."""
    # 如果有命令行参数，则使用命令行参数；否则使用默认配置
    if len(os.sys.argv) > 1:
        # 处理命令行参数
        parser = argparse.ArgumentParser(description="Plot traffic network path costs from cached MATLAB data.")
        parser.add_argument('--zeta', type=int, required=False, default=DEFAULT_CONFIG['zeta'],
                            help=f'Zeta value for the dataset (default: {DEFAULT_CONFIG["zeta"]}).')
        parser.add_argument('--subset_index', type=int, required=False, default=DEFAULT_CONFIG['subset_index'],
                            help=f'Subset index for the dataset (default: {DEFAULT_CONFIG["subset_index"]}).')
        parser.add_argument('--num_flows', type=int, default=DEFAULT_CONFIG['num_flows'],
                            help=f'Number of random flow vectors to plot (default: {DEFAULT_CONFIG["num_flows"]}).')
        parser.add_argument('--cache_dir', type=str, default=DEFAULT_CONFIG['cache_dir'],
                            help=f'Directory containing cache files (default: {DEFAULT_CONFIG["cache_dir"]}).')
        parser.add_argument('--results_dir', type=str, default=DEFAULT_CONFIG['results_dir'],
                            help=f'Directory to save plots (default: {DEFAULT_CONFIG["results_dir"]}).')
        args = parser.parse_args()
        
        # 使用命令行参数运行
        run_with_params(
            zeta=args.zeta,
            subset_index=args.subset_index,
            num_flows=args.num_flows,
            cache_dir=args.cache_dir,
            results_dir=args.results_dir
        )
    else:
        # 批量运行预设的 zeta 与 subset_index 组合，避免重复代码
        default_run_configs = [
            # {"zeta": 8, "subset_index": 0},
            # {"zeta": 16, "subset_index": 0},
            # {"zeta": 16, "subset_index": 1},
            {"zeta": 24, "subset_index": 0},
            # {"zeta": 24, "subset_index": 1},
            # {"zeta": 32, "subset_index": 0},
            # {"zeta": 32, "subset_index": 1},
            # {"zeta": 32, "subset_index": 2},
        ]

        common_kwargs = dict(
            num_flows=10000,
            cache_dir="matlab/cache",
            results_dir="results",
            plot_params=PlotParams(save_path="results/", figure_dpi=600),
        )

        for cfg in default_run_configs:
            run_with_params(**cfg, **common_kwargs)


# 示例用法 - 在直接运行脚本时使用
if __name__ == '__main__':
    main()
