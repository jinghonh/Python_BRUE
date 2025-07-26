import argparse
import os
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


def configure_plot(ax: plt.Axes, params: PlotParams, title: str, xlabel: str, ylabel: str):
    """Apply common plot configurations."""
    ax.set_title(title, fontsize=params.font_size + 4, fontweight='bold')
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
    configure_plot(ax, params, 'Path Time-Money Cost Relationship', r'Time Cost ($T$)', r'Money Cost ($M$)')

    if boundary.left_x.size > 0:
        boundary_x = np.concatenate([boundary.left_x, boundary.right_x[::-1]])
        boundary_y = np.concatenate([boundary.left_y, boundary.right_y[::-1]])
        ax.fill(boundary_x, boundary_y, color=[0.9, 0.95, 1], alpha=1, edgecolor='none', label='Feasible Region')
        ax.plot(boundary.left_x, boundary.left_y, '-', color=[0.4, 0.5, 0.8], linewidth=0.5)
        ax.plot(boundary.right_x, boundary.right_y, '-', color=[0.4, 0.5, 0.8], linewidth=0.5)

    for i, costs in enumerate(all_path_costs):
        if costs.size > 0:
            ax.plot(costs[:, 0], costs[:, 1], '-', color=(*color_variations[i], 0.7), linewidth=1.2)
            ax.scatter(costs[:, 0], costs[:, 1], s=30, c=[color_variations[i]], marker='o', alpha=0.8,
                       edgecolors='none', label=f'Flow Vector {i+1}' if i < 5 else "")

    handles, labels = ax.get_legend_handles_labels()
    if handles:
        unique_labels = dict(zip(labels, handles))
        ax.legend(unique_labels.values(), unique_labels.keys(), loc='best', fontsize=params.font_size - 1)
    
    save_figure(fig, 'path_time_money', params, zeta, subset_index)
    plt.close(fig)


def plot_path_costs_with_upper_limit(all_path_costs: List[np.ndarray], boundary: Boundary,
                                     color_variations: np.ndarray, params: PlotParams, zeta: int, subset_index: int):
    """Plot path costs with an upper limit and equilibrium line, with correct layering."""
    fig, ax = plt.subplots(figsize=params.figure_size)
    configure_plot(ax, params, 'Path Costs with Upper Limit', r'Time Cost ($T$)', r'Money Cost ($M$)')

    # --- Data Processing ---
    if boundary.left_x.size == 0:
        print("Warning: Boundary data is empty. Skipping plot generation.")
        plt.close(fig)
        return

    t_max = (boundary.left_x + boundary.right_x) / 2
    t_eqm = t_max - 50. / (10 + boundary.left_y)
    f_tmax = interp1d(boundary.left_y, t_max, kind='nearest', bounds_error=False, fill_value="extrapolate")

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
    for path in infeasible_paths_data:
        costs = path['costs']
        color = path['color']
        ax.plot(costs[:, 0], costs[:, 1], '-', color=(*color[:3], 0.5), linewidth=1.0, zorder=2)
        ax.scatter(costs[:, 0], costs[:, 1], s=25, c=[color], marker='o', alpha=0.3, edgecolors='none', zorder=2)

    # Layer 3: Feasible Flow Sub-Region
    if feasible_costs_points:
        all_feasible_time = np.concatenate([c[:, 0] for c in feasible_costs_points])
        all_feasible_money = np.concatenate([c[:, 1] for c in feasible_costs_points])
        if all_feasible_time.size > 2:
            feasible_sub_boundary = calculate_feasible_region_boundary(all_feasible_time, all_feasible_money)
            if feasible_sub_boundary.left_x.size > 0:
                sub_boundary_x = np.concatenate([feasible_sub_boundary.left_x, feasible_sub_boundary.right_x[::-1]])
                sub_boundary_y = np.concatenate([feasible_sub_boundary.left_y, feasible_sub_boundary.right_y[::-1]])
                ax.fill(sub_boundary_x, sub_boundary_y, color=[0.7, 0.95, 0.7], alpha=0.5,
                        edgecolor=[0.4, 0.6, 0.4], linewidth=0.5, label='Feasible Flow Region', zorder=3)

    # Layer 4: Feasible Paths
    for path in feasible_paths_data:
        costs = path['costs']
        color = path['color']
        ax.plot(costs[:, 0], costs[:, 1], '-', color=(*color[:3], 0.7), linewidth=1.2, zorder=4)
        ax.scatter(costs[:, 0], costs[:, 1], s=30, c=[color], marker='o', alpha=0.8, edgecolors='none', zorder=4)

    # Layer 5: T_max and T_eqm lines (Topmost)
    ax.plot(t_max, boundary.left_y, '-.', color=[0.8, 0.2, 0.2], linewidth=2, label=r'$T_{max}$')
    ax.plot(t_eqm, boundary.left_y, '-.', color=[0.5, 0.0, 0.8], linewidth=2, label=r'$T_{eqm}$')

    save_figure(fig, 'path_costs_upper_limit', params, zeta, subset_index)
    plt.close(fig)


def plot_path_costs_below_equilibrium(all_path_costs: List[np.ndarray], boundary: Boundary,
                                      color_variations: np.ndarray, params: PlotParams, zeta: int, subset_index: int):
    """筛选并绘制T_eqm（紫色折线）下方的可行方案。"""
    fig, ax = plt.subplots(figsize=params.figure_size)
    configure_plot(ax, params, 'Path Costs Below Equilibrium Line', r'Time Cost ($T$)', r'Money Cost ($M$)')

    # --- Data Processing ---
    if boundary.left_x.size == 0:
        print("Warning: Boundary data is empty. Skipping plot generation.")
        plt.close(fig)
        return

    # 计算T_max和T_eqm
    t_max = (boundary.left_x + boundary.right_x) / 2
    t_eqm = t_max - 50. / (10 + boundary.left_y)
    
    # 创建T_eqm的插值函数用于判断
    f_teqm = interp1d(boundary.left_y, t_eqm, kind='nearest', bounds_error=False, fill_value="extrapolate")

    # 筛选出所有点都在T_eqm下方的路径
    below_eqm_paths_data = []
    below_eqm_costs_points = []
    other_paths_data = []

    for i, costs in enumerate(all_path_costs):
        if costs.size > 0:
            # 判断所有点是否都在T_eqm下方
            is_below_eqm = np.all(costs[:, 0] <= f_teqm(costs[:, 1]))
            if is_below_eqm:
                # 使用紫色系颜色（接近原T_eqm的颜色但略有变化）
                purple_color = np.array([0.5, 0.0, 0.8]) + np.array([0.1, 0.1, 0.0]) * (i % 3)
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
        ax.plot(costs[:, 0], costs[:, 1], '-', color=(*color[:3], 0.7), linewidth=1.5, zorder=4)
        ax.scatter(costs[:, 0], costs[:, 1], s=35, c=[color], marker='o', alpha=0.8, edgecolors='none', zorder=4)

    # Layer 5: T_eqm line (Topmost)
    ax.plot(t_eqm, boundary.left_y, '-.', color=[0.5, 0.0, 0.8], linewidth=2.5, label=r'$T_{eqm}$', zorder=5)

    # 添加图例
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        unique_labels = dict(zip(labels, handles))
        ax.legend(unique_labels.values(), unique_labels.keys(), loc='best', fontsize=params.font_size - 1)

    save_figure(fig, 'path_costs_below_equilibrium', params, zeta, subset_index)
    plt.close(fig)


def main():
    """Main function to load data and generate plots."""
    parser = argparse.ArgumentParser(description="Plot traffic network path costs from cached MATLAB data.")
    parser.add_argument('--zeta', type=int, required=True, help='Zeta value for the dataset.')
    parser.add_argument('--subset_index', type=int, required=True, help='Subset index for the dataset.')
    parser.add_argument('--num_flows', type=int, default=1000, help='Number of random flow vectors to plot.')
    parser.add_argument('--cache_dir', type=str, default='matlab/cache', help='Directory containing cache files.')
    parser.add_argument('--results_dir', type=str, default='results', help='Directory to save plots.')
    args = parser.parse_args()

    # Setup plotting style
    plt.rcParams['text.usetex'] = True
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['font.serif'] = ['Computer Modern Roman', 'Times New Roman']
    plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'

    params = PlotParams(save_path=args.results_dir)

    # Load data from .mat file
    cache_file = os.path.join(args.cache_dir, f'cache_zeta{args.zeta}_subset{args.subset_index}.mat')
    if not os.path.exists(cache_file):
        print(f"Error: Cache file not found at {cache_file}")
        return

    print(f"Loading data from {cache_file}...")
    data = loadmat(cache_file)
    total_valid_flow = data.get('totalValidFlow')
    relation_matrix = data.get('relationMatrix')

    if total_valid_flow is None or relation_matrix is None:
        print("Error: 'totalValidFlow' or 'relationMatrix' not found in cache file.")
        return
        
    print("Data loaded successfully.")

    # Select random flow indices
    num_total_flows = total_valid_flow.shape[0]
    if num_total_flows == 0:
        print("No valid flow data found in the cache file.")
        return
        
    selected_indices = np.random.permutation(num_total_flows)[:min(args.num_flows, num_total_flows)]

    # Define constants
    money_coeffs = np.array([20, 15, 1, 0, 0, 0, 0, 1])
    free_flow_time = np.array([18, 22.5, 12, 24, 2.4, 6, 24, 12])
    max_capacity = np.array([3600, 3600, 1800, 1800, 1800, 1800, 1800, 1800])

    # Prepare data
    all_path_costs, all_path_time_costs, all_path_money_costs, color_variations = \
        prepare_path_costs_data(total_valid_flow, selected_indices, relation_matrix,
                                money_coeffs, free_flow_time, max_capacity)

    # Calculate overall boundary
    boundary = calculate_feasible_region_boundary(all_path_time_costs, all_path_money_costs)

    # Generate plots
    print("Generating plots...")
    plot_time_money_cost_relationship(all_path_costs, boundary, color_variations, params, args.zeta, args.subset_index)
    plot_path_costs_with_upper_limit(all_path_costs, boundary, color_variations, params, args.zeta, args.subset_index)
    plot_path_costs_below_equilibrium(all_path_costs, boundary, color_variations, params, args.zeta, args.subset_index)
    print("Plots generated successfully.")


if __name__ == '__main__':
    main()