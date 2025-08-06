import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
# --- Matplotlib helpers
from matplotlib.lines import Line2D
import matplotlib.patheffects as pe  # 用于给文字添加描边效果，避免覆盖标记
from matplotlib.colors import to_rgb, to_hex

# 统一样式配置
from plot_styles import REGION_STYLES
import json
import os
from scipy.io import loadmat  # 添加导入.mat文件的库

# 读取CSV数据
def load_data(file_path):
    return pd.read_csv(file_path)

# 读取MATLAB的.mat文件
def load_mat_data(file_path):
    return loadmat(file_path)


# ----------------------------------------------------------------------
# 辅助函数：根据给定颜色生成浅/深色变体
# ----------------------------------------------------------------------


def _lighten_color(color: str, amount: float = 0.4) -> str:
    """Blend color towards white by given amount (0~1)."""
    r, g, b = to_rgb(color)
    r = r + (1 - r) * amount
    g = g + (1 - g) * amount
    b = b + (1 - b) * amount
    return to_hex((r, g, b))


def _darken_color(color: str, amount: float = 0.25) -> str:
    """Blend color towards black by given amount (0~1)."""
    r, g, b = to_rgb(color)
    r *= (1 - amount)
    g *= (1 - amount)
    b *= (1 - amount)
    return to_hex((r, g, b))
    

# 计算可行域边界
def calculate_feasible_region_boundary(time_costs, money_costs):
    """
    计算可行域边界，基于路径成本
    
    参数:
        time_costs: 时间成本数组
        money_costs: 金钱成本数组
        
    返回:
        boundary: 包含左右边界的字典
    """
    # 如果没有成本数据，返回默认边界
    if len(time_costs) == 0 or len(money_costs) == 0:
        return {
            'leftX': np.array([0, 10]),
            'leftY': np.array([10, 0]),
            'rightX': np.array([20, 10]),
            'rightY': np.array([0, 10])
        }
    
    # 找到唯一的金钱成本值
    unique_money_values = np.unique(money_costs)
    left_boundary_x = []
    left_boundary_y = []
    right_boundary_x = []
    right_boundary_y = []
    
    # 对于每个金钱成本值，找到最小和最大时间成本
    for money in unique_money_values:
        same_money_idx = np.abs(money_costs - money) < 0.001
        
        if np.sum(same_money_idx) > 0:
            times_for_money = time_costs[same_money_idx]
            min_time = np.min(times_for_money)
            max_time = np.max(times_for_money)
            
            # 添加到边界数组
            left_boundary_x.append(min_time)
            left_boundary_y.append(money)
            right_boundary_x.append(max_time)
            right_boundary_y.append(money)
    
    # 按金钱成本排序边界点
    sort_idx = np.argsort(left_boundary_y)
    left_boundary_y = np.array(left_boundary_y)[sort_idx]
    left_boundary_x = np.array(left_boundary_x)[sort_idx]
    
    sort_idx = np.argsort(right_boundary_y)
    right_boundary_y = np.array(right_boundary_y)[sort_idx]
    right_boundary_x = np.array(right_boundary_x)[sort_idx]
    
    # 创建边界字典
    boundary = {
        'leftX': left_boundary_x,
        'leftY': left_boundary_y,
        'rightX': right_boundary_x,
        'rightY': right_boundary_y
    }
    
    return boundary

# 计算路径成本
def calculate_all_path_costs(flow, relation_matrix):
    """
    计算所有路径成本
    
    参数:
        flow: 流量数组
        relation_matrix: 关系矩阵
        
    返回:
        path_costs: 路径成本数组 [时间成本, 金钱成本]
    """
    # 获取相关常量
    money = np.array([20, 15, 1, 0, 0, 0, 0, 1])
    free_flow_time = np.array([18, 22.5, 12, 24, 2.4, 6, 24, 12])
    max_capacity = np.array([3600, 3600, 1800, 1800, 1800, 1800, 1800, 1800])
    
    # 计算链路流量
    x = np.dot(flow, relation_matrix)
    
    # 计算实际行驶时间
    index = np.where(np.sum(relation_matrix, axis=0) != 0)[0]
    time = free_flow_time[index] * (1 + 0.15 * (x[index] / max_capacity[index])**4)
    
    # 计算每条路径的时间成本
    path_time = np.dot(time, relation_matrix[:, index].T)
    RT = path_time + 15 * (1 - np.exp(-0.02 * path_time))
    
    # 计算每条路径的货币成本
    path_money = np.dot(money, relation_matrix.T)
    
    # 组合时间和货币成本
    path_costs = np.column_stack((RT, path_money))
    return path_costs

# 从.mat文件计算TTB_max_delta上限
def calculate_tmax_from_mat(mat_file_path):
    """
    从.mat文件中读取数据并计算TTB_max_delta上限
    
    参数:
        mat_file_path: .mat文件路径
        
    返回:
        upper_limit_x: TTB_max_delta上限的x坐标
        upper_limit_y: TTB_max_delta上限的y坐标
    """
    try:
        # 读取.mat文件
        mat_data = load_mat_data(mat_file_path)
        
        # 提取数据
        total_valid_flow = mat_data['totalValidFlow']
        relation_matrix = mat_data['relationMatrix']
        
        # 收集所有路径成本
        all_time_costs = []
        all_money_costs = []
        
        # 对所有流量计算路径成本
        for i in range(total_valid_flow.shape[0]):
            flow = total_valid_flow[i, :]
            path_costs = calculate_all_path_costs(flow, relation_matrix)
            
            # 只收集非零流量路径的成本
            non_zero_flow_idx = np.where(flow > 0)[0]
            if len(non_zero_flow_idx) > 0:
                # 提取非零流量路径的成本
                non_zero_costs = path_costs[non_zero_flow_idx, :]
                # 排除金钱成本为0的路径
                non_zero_money_idx = non_zero_costs[:, 1] > 0
                if np.any(non_zero_money_idx):
                    all_time_costs.extend(non_zero_costs[non_zero_money_idx, 0])
                    all_money_costs.extend(non_zero_costs[non_zero_money_idx, 1])
        
        # 计算可行域边界
        all_time_costs = np.array(all_time_costs)
        all_money_costs = np.array(all_money_costs)
        boundary = calculate_feasible_region_boundary(all_time_costs, all_money_costs)
        
        # 计算TTB_max_delta上限（使用边界中点）
        upper_limit_x = (boundary['rightX'] + boundary['leftX']) / 2
        upper_limit_y = boundary['rightY']
        
        # 确保长度一致
        min_len = min(len(upper_limit_x), len(upper_limit_y))
        upper_limit_x = upper_limit_x[:min_len]
        upper_limit_y = upper_limit_y[:min_len]

        print("upper_limit_x:", upper_limit_x)
        print("T_eqm:", upper_limit_x-(50/(10+upper_limit_y)))
        print("upper_limit_y:", upper_limit_y)
        
        return upper_limit_x, upper_limit_y
        
    except Exception as e:
        print(f"从.mat文件计算TTB_max_delta时出错: {e}")
        return np.array([]), np.array([])

# 绘制三区域路径成本对比图
def plot_three_regions_comparison(df, zeta_value, figsize=(10, 8), use_mat_file=False):
    # 过滤指定zeta值的数据
    df_zeta = df[df['zeta值'] == zeta_value]
    
    # 创建图形
    fig, ax = plt.subplots(figsize=figsize)
    
    # 使用统一样式（根据需要映射区域编号 -> 样式键）
    region_key_map = {1: "RBS", 2: "BS", 3: "S"}
    colors = {k: REGION_STYLES[v].color for k, v in region_key_map.items()}
    markers = {k: REGION_STYLES[v].marker for k, v in region_key_map.items()}

    legend_handles = []
    
    # 计算图表的x轴范围，用于后续添加标签
    x_min, x_max = float('inf'), float('-inf')
    
    # 获取所有路径编号
    unique_paths = df_zeta['路径编号'].unique()
    
    # 处理每个区域的数据
    for region in [1, 2, 3]:
        for scheme in [1, 2]:
            # 获取当前方案的数据
            scheme_data = df_zeta[(df_zeta['区域'] == region) & (df_zeta['方案编号'] == scheme)]
            
            if not scheme_data.empty:
                # 提取路径的时间成本和金钱成本
                time_costs = scheme_data['时间成本'].values
                money_costs = scheme_data['金钱成本'].values
                flows = scheme_data['流量分配'].values
                path_ids = scheme_data['路径编号'].values
                
                # 更新x轴范围
                x_min = min(x_min, np.min(time_costs))
                x_max = max(x_max, np.max(time_costs))
                
                # 绘制散点
                scatter = ax.scatter(
                    time_costs,
                    money_costs,
                    facecolor=_lighten_color(colors[region], 0.4),
                    edgecolor=_darken_color(colors[region], 0.25),
                    marker=markers[region],
                    s=60,
                    alpha=0.9,
                    linewidths=0.8,
                )
                
                # 对于同一方案中的每个唯一金钱成本（即同一路径），连接这些点
                for path_id in unique_paths:
                    # 获取当前路径编号的数据点
                    path_points = scheme_data[scheme_data['路径编号'] == path_id]
                    
                    if not path_points.empty:
                        path_time_costs = path_points['时间成本'].values
                        path_money_cost = path_points['金钱成本'].values[0]  # 同一路径的金钱成本应该相同
                        
                        # 只绘制水平线，如果有多个点
                        if len(path_time_costs) > 1:
                            # 按时间成本排序
                            sorted_indices = np.argsort(path_time_costs)
                            sorted_time_costs = path_time_costs[sorted_indices]
                            
                            # 绘制水平连接线
                            ax.plot([sorted_time_costs[0], sorted_time_costs[-1]], 
                                    [path_money_cost, path_money_cost],
                                    '-', color=colors[region], linewidth=1.0, alpha=0.6)
                
                # 连接流量大于0的路径点
                pos_flow_indices = np.where(flows > 0)[0]
                if len(pos_flow_indices) > 0:
                    # 提取流量大于0的路径成本
                    pos_time_costs = time_costs[pos_flow_indices]
                    pos_money_costs = money_costs[pos_flow_indices]
                    
                    # 按金钱成本排序
                    sorted_indices = np.argsort(pos_money_costs)
                    sorted_time_costs = pos_time_costs[sorted_indices]
                    sorted_money_costs = pos_money_costs[sorted_indices]
                    
                    # 绘制连接线
                    ax.plot(sorted_time_costs, sorted_money_costs, '-', color=colors[region], linewidth=1.5)
                
                # 为图例创建标签
                if region == 1:
                    label = f"$S_0^{{\\varepsilon}}$ {scheme}"
                elif region == 2:
                    label = f"$BS_0^{{\\varepsilon}}$ {scheme}"
                else:
                    label = f"$RBS_0^{{\\varepsilon}}$ {scheme}"
                
                # 添加到图例句柄
                legend_handles.append(mpatches.Patch(color=colors[region], label=label))
    
    # 绘制 TTB_max_delta 上限线
    if use_mat_file:
        # 从.mat文件计算TTB_max_delta上限
        mat_file_path = f'matlab/cache/cache_zeta{zeta_value}_subset0.mat'
        upper_limit_x, upper_limit_y = calculate_tmax_from_mat(mat_file_path)
        
        if len(upper_limit_x) > 0 and len(upper_limit_y) > 0:
            # 排序以确保正确连线
            sort_idx = np.argsort(upper_limit_y)
            upper_limit_x = upper_limit_x[sort_idx]
            upper_limit_y = upper_limit_y[sort_idx]
            
            # 绘制TTB_max_delta线
            ax.plot(upper_limit_x, upper_limit_y, '-.', color=[0.8, 0.2, 0.2], linewidth=2)
        else:
            print(f"警告: 无法从.mat文件计算TTB_max_delta上限")
    else:
        # 使用CSV数据中的TTB_max_delta上限
        tmax_data = df_zeta[df_zeta['区域'] == 4]  # 使用区域4表示 TTB_max_delta 数据
        if not tmax_data.empty:
            # 按金钱成本排序
            tmax_data = tmax_data.sort_values('金钱成本')
            # 绘制实际 TTB_max_delta 线
            ax.plot(tmax_data['时间成本'].values, tmax_data['金钱成本'].values, 
                    '-.', color=[0.8, 0.2, 0.2], linewidth=2)
        else:
            x_range = np.linspace(min(df_zeta['时间成本'])*0.8, max(df_zeta['时间成本'])*1.2, 100)
            y_range = np.linspace(max(df_zeta['金钱成本'])*0.8, min(df_zeta['金钱成本'])*1.2, 100)
            ax.plot(x_range, y_range, '-.', color=[0.8, 0.2, 0.2], linewidth=2)
    
    # 在右侧添加路径标签
    margin = (x_max - x_min) * 0.05  # 设置一个边距
    # 获取唯一的金钱成本值并排序
    unique_money_costs = sorted(df_zeta['金钱成本'].unique(), reverse=True)
    
    for money_cost in unique_money_costs:
        # 对于每个金钱成本，找到对应的路径编号
        path_data = df_zeta[df_zeta['金钱成本'] == money_cost]
        if not path_data.empty:
            path_id = path_data['路径编号'].values[0]
            # 在右侧添加路径标签
            ax.text(x_max + margin, money_cost, f'$P_{path_id}$', 
                    fontsize=10, va='center', ha='left', color='#4D4D4D')
    
    # 设置轴标签和图例
    ax.set_xlabel('Time', fontsize=12)
    ax.set_ylabel('Money Cost', fontsize=12)
    ax.grid(True)
    
    # 扩展x轴范围以容纳路径标签
    ax.set_xlim(x_min - margin, x_max + margin * 5)
    
    # 添加图例
    handles = legend_handles + [Line2D([0], [0], linestyle='-.', color=[0.8, 0.2, 0.2], label=r'$TTB_{max} + \delta_k$')]
    ax.legend(handles=handles, loc='best', fontsize=10, framealpha=1)
    
    
    plt.tight_layout()
    plt.savefig(f'three_regions_comparison_zeta{zeta_value}.pdf', format='pdf', dpi=300)
    plt.show()

# Plot path cost comparison using scatter JSON and compute all path costs.
def plot_comparison_path_costs(zeta_value, figsize=(10, 8)):
    """Plot path cost comparison using scatter JSON and compute all path costs."""
    # Load scatter JSON
    scatter_file = os.path.join('results', f'scatter_points_e{zeta_value}.json')
    try:
        with open(scatter_file, 'r') as f:
            scatter_data = json.load(f)
    except FileNotFoundError:
        print(f"Warning: Scatter file not found: {scatter_file}")
        return
    # Load relation matrix from MATLAB cache for full path set
    mat_file = os.path.join('matlab', 'cache', f'cache_zeta{32}_subset2.mat')
    try:
        mat_data = load_mat_data(mat_file)
        relation_matrix = mat_data['relationMatrix']
        
    except Exception as e:
        print(f"Warning: Failed to load relation matrix: {e}")
        return

    fig, ax = plt.subplots(figsize=figsize)
    # Define regions: JSON key, color, legend label, and marker for scattering other paths
    regions = [
        ('s_constraint_points', REGION_STYLES['S'].color, REGION_STYLES['S'].label, REGION_STYLES['S'].marker),
        ('all_constraint_points', REGION_STYLES['BS'].color, REGION_STYLES['BS'].label, REGION_STYLES['BS'].marker),
        ('path_constraint_points', REGION_STYLES['RBS'].color, REGION_STYLES['RBS'].label, REGION_STYLES['RBS'].marker)
    ]
    total_flow = 10000.0
    # storage for csv rows
    csv_rows = []
    region_short = {
        's_constraint_points': 'S',
        'all_constraint_points': 'BS',
        'path_constraint_points': 'RBS'
    }
    # For each region: draw two sample polyline for P1,P2,P5 and scatter other paths
    for key, color, label, marker in regions:
        pts_list = scatter_data.get(key, [])
        for idx_pt, pt in enumerate(pts_list):
            f1, f2 = pt
            # Determine flows on paths: f1, f2, f5, others zero
            f5 = total_flow - f1 - f2
            n_paths = relation_matrix.shape[0]
            flow = np.zeros(n_paths)
            flow[0] = f1  # path 1
            flow[1] = f2  # path 2
            # path5 index = 4 in full 6-path set
            if n_paths >= 5:
                flow[4] = f5  # path 5
            # Compute costs for all paths
            costs = calculate_all_path_costs(flow, relation_matrix)
            times = costs[:, 0]
            monies = costs[:, 1]
            # Label only the first sample
            lbl = label if idx_pt == 0 else None
            # Connect only paths 1,2,5
            line_idx = [0, 1, 4]
            times_line = times[line_idx]
            monies_line = monies[line_idx]
            ax.plot(times_line, monies_line, '-', color=color, linewidth=1.5, marker=marker, markersize=8, label=lbl)
            # Scatter other paths
            other_idx = [i for i in range(len(times)) if i not in line_idx]
            ax.scatter(
                times[other_idx],
                monies[other_idx],
                facecolor=_lighten_color(color, 0.4),
                edgecolor=_darken_color(color, 0.25),
                marker=marker,
                s=50,
                linewidths=0.6,
            )

            # --- collect csv rows ---
            region_code = region_short.get(key, key)
            for p_idx in range(n_paths):
                csv_rows.append({
                    'region': region_code,
                    'sample': idx_pt + 1,
                    'path_id': p_idx + 1,
                    'flow': flow[p_idx],
                    'time_cost': times[p_idx],
                    'money_cost': monies[p_idx]
                })

    # Plot TTB_max_delta upper bound from JSON
    tmax_file = os.path.join('results', f'tmax_teqm_zeta{zeta_value}.json')
    try:
        with open(tmax_file, 'r') as f:
            tmax_data = json.load(f)
        mvals = np.array(tmax_data['money_values'], dtype=float)
        tmax_arr = np.array(tmax_data['TTB_max_delta'], dtype=float)
        idx_t = np.argsort(mvals)
        ax.plot(tmax_arr[idx_t], mvals[idx_t], '-', color=[0.8, 0.2, 0.2], linewidth=3, label=r'$TTB_{max} + \delta_k$')
    except FileNotFoundError:
        pass
    # --- annotate path labels at rightmost occurrence ---
    path_points = {}  # 存储每个路径ID对应的所有点
    
    # 收集每个路径ID对应的所有点
    for row in csv_rows:
        pid = row['path_id']
        t = row['time_cost']
        m = row['money_cost']
        if pid not in path_points:
            path_points[pid] = []
        path_points[pid].append((t, m))
    
    if path_points:
        # 找出所有点的最大时间成本用于放置标签
        all_time_costs = [max(points, key=lambda x: x[0])[0] for points in path_points.values()]
        max_time = max(all_time_costs) if all_time_costs else 0
        
        for pid, points in path_points.items():
            # 对每个路径ID找到时间成本最大和最小的点
            max_point = max(points, key=lambda x: x[0])
            min_point = min(points, key=lambda x: x[0])
            
            # 获取金钱成本和时间成本
            max_t, m = max_point  # 金钱成本使用最右侧点的值
            min_t, _ = min_point  # 最小时间成本点
            
            # 绘制水平虚线（从最左侧点到标签位置，与标签保持相同高度）
            ax.plot([min_t, max_time + 0.5], [m, m], color='gray', linestyle=':', linewidth=3)
            
            # 添加路径标签
            ax.text(max_time + 0.6, m, f'$P_{pid}$', fontsize=24, va='center', ha='left')

    # Finalize and save
    ax.set_xlabel('Time', fontsize=12)
    ax.set_ylabel('Money Cost', fontsize=12)
    ax.grid(True)
    ax.legend(loc='best', fontsize=10, framealpha=1)
    # 
    plt.tight_layout()
    plt.savefig(f'results/path_costs_comparison_zeta{zeta_value}.pdf', format='pdf', dpi=300)
    # plt.show()

    # Save CSV
    if csv_rows:
        import pandas as _pd
        df_out = _pd.DataFrame(csv_rows)
        os.makedirs('results', exist_ok=True)
        csv_path = os.path.join('results', f'flow_costs_zeta{zeta_value}.csv')
        df_out.to_csv(csv_path, index=False)
        print(f"Flow and cost data saved to {csv_path}")

# 绘制前两个区域的成本折线图，不包括Tmax
def plot_two_regions_path_costs(zeta_value, figsize=(10, 8)):
    """
    绘制仅包含前两个区域(S和BS)的成本折线图，不绘制Tmax上限
    
    参数:
        zeta_value: epsilon值
        figsize: 图形大小
    """
    # 加载散点JSON数据
    scatter_file = os.path.join('results', f'scatter_points_e{zeta_value}.json')
    try:
        with open(scatter_file, 'r') as f:
            scatter_data = json.load(f)
    except FileNotFoundError:
        print(f"警告: 未找到散点文件: {scatter_file}")
        return
    
    # 从MATLAB缓存加载关系矩阵
    mat_file = os.path.join('matlab', 'cache', f'cache_zeta{32}_subset2.mat')
    try:
        mat_data = load_mat_data(mat_file)
        relation_matrix = mat_data['relationMatrix']
    except Exception as e:
        print(f"警告: 加载关系矩阵失败: {e}")
        return

    fig, ax = plt.subplots(figsize=figsize)
    
    # 仅定义前两个区域：JSON键、颜色、图例标签和散点标记
    regions = [
        ('s_constraint_points', REGION_STYLES['S'].color, REGION_STYLES['S'].label, REGION_STYLES['S'].marker),
        ('all_constraint_points', REGION_STYLES['BS'].color, REGION_STYLES['BS'].label, REGION_STYLES['BS'].marker)
    ]
    total_flow = 10000.0
    
    # 存储CSV行的列表
    csv_rows = []
    region_short = {
        's_constraint_points': 'S',
        'all_constraint_points': 'BS'
    }

    # 存储每个区域每个方案的时间成本范围
    region_time_ranges = {}
    
    # 存储图例句柄
    legend_handles = []
    
    # 定义字母标签和存储点的位置
    point_labels = ['D', 'C', 'A', 'B']  # 点的标签
    point_label_idx = 0  # 用于跟踪标签索引
    labeled_points = []  # 存储带标签的点: [(x, y, label, color), ...]

    # 对于每个区域
    for key, color, label, marker in regions:
        pts_list = scatter_data.get(key, [])
        region_time_ranges[key] = []  # 用于存储该区域的所有方案的时间成本范围
        region_legend_added = False   # 跟踪该区域的图例是否已添加
        
        for idx_pt, pt in enumerate(pts_list):
            f1, f2 = pt
            # 确定路径上的流量：f1、f2、f5，其他为零
            f5 = total_flow - f1 - f2
            n_paths = relation_matrix.shape[0]
            flow = np.zeros(n_paths)
            flow[0] = f1  # 路径1
            flow[1] = f2  # 路径2
            # 路径5的索引在完整的6路径集中是4
            if n_paths >= 5:
                flow[4] = f5  # 路径5
            
            # 计算所有路径的成本
            costs = calculate_all_path_costs(flow, relation_matrix)
            times = costs[:, 0]
            monies = costs[:, 1]
            
            # 只处理有流量的路径
            flow_paths_idx = np.where(flow > 0)[0]
            
            if len(flow_paths_idx) > 0:
                # 提取有流量路径的成本
                flow_times = times[flow_paths_idx]
                flow_monies = monies[flow_paths_idx]
                
                # 按金钱成本排序，用于绘制连接线
                sorted_indices = np.argsort(flow_monies)
                sorted_times = flow_times[sorted_indices]
                sorted_monies = flow_monies[sorted_indices]
                
                # 存储这个方案的时间成本范围
                time_min = np.min(sorted_times)
                time_max = np.max(sorted_times)
                # 获取当前字母标签（若超出范围则为空字符串）
                current_letter = point_labels[point_label_idx - 1] if point_label_idx - 1 < len(point_labels) else ''
                region_time_ranges[key].append((time_min, time_max, sorted_monies[0], idx_pt, current_letter))
                
                # 仅为每个区域的第一个有效方案添加图例
                line_label = label if not region_legend_added else None
                if line_label:
                    region_legend_added = True
                
                # 绘制有流量路径的连接线和散点
                light_c = _lighten_color(color, 0.4)
                dark_c = _darken_color(color, 0.25)
                line = ax.plot(
                    sorted_times,
                    sorted_monies,
                    '-',
                    color=color,  # 线条保持基准色
                    linewidth=1.5,
                    marker=marker,
                    markersize=25,
                    label=line_label,
                    alpha=0.85,
                    markeredgecolor=dark_c,
                    markeredgewidth=1.2,
                    markerfacecolor=light_c,
                    zorder=99,
                )
                
                # 如果是第一个有标签的线，保存图例句柄
                if line_label:
                    legend_handles.append(line[0])
                
                # 为这条折线的所有散点添加同一个字母标签（例如A、B、C、D）
                # 只有当还有可用标签时才添加
                if point_label_idx < len(point_labels):
                    current_letter = point_labels[point_label_idx]
                    for t_val, m_val in zip(sorted_times, sorted_monies):
                        labeled_points.append((t_val, m_val, current_letter, dark_c))
                    # 一个折线使用一个字母，之后切换到下一个字母
                    point_label_idx += 1
 
                # --- 收集CSV行 ---
                region_code = region_short.get(key, key)
                for p_idx in flow_paths_idx:
                    csv_rows.append({
                        'region': region_code,
                        'sample': idx_pt + 1,
                        'path_id': p_idx + 1,
                        'flow': flow[p_idx],
                        'time_cost': times[p_idx],
                        'money_cost': monies[p_idx]
                    })
    
    # 绘制带有字母标签的点
    for time, money, label_text, point_color in labeled_points:
        txt = ax.text(
            time,
            money-0.1,
            label_text,
            fontsize=22,
            fontweight='bold',
            ha='center',
            va='center',
            color=point_color,  # 与 marker 边色一致
            zorder=120,  # 确保文字在最上层，但不遮挡标记边界
        )
        # 添加白色描边使文字在彩色标记上清晰可见
        # txt.set_path_effects([pe.withStroke(linewidth=2, foreground='white')])

    # 使用更好的视觉效果显示时间成本区间
    ax_top = ax.twiny()  # 创建一个共享y轴的新x轴
    ax_top.set_xlabel('')
    ax_top.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False, labeltop=False)  # 隐藏刻度和标签
    
    # 获取y轴范围，用于绘制区间
    y_max = ax.get_ylim()[1]
    y_min = ax.get_ylim()[0]
    y_range = y_max - y_min
    
    region_colors = {
        's_constraint_points': REGION_STYLES['S'].color,
        'all_constraint_points': REGION_STYLES['BS'].color,
    }
    
    # 区间位置的相对高度（S区域放在上方，BS区域放在下方）
    interval_positions = {
        's_constraint_points': {0: 0.80, 1: 0.75},  # S区域的两个区间放在图表上方，稍微向下移动
        'all_constraint_points': {0: 0.15, 1: 0.10}   # BS区域的两个区间放在图表下方
    }
    
    # 获取x轴范围，用于区间宽度的缩放
    x_min, x_max = ax.get_xlim()
    x_range = x_max - x_min
    
    # 绘制区间并添加标签
    for key, ranges in region_time_ranges.items():
        for i, (t_min, TTB_max_delta, money, idx, letter_tag) in enumerate(ranges):
            if i < 2:  # 只处理每个区域的前两个方案
                y_pos = y_min + y_range * interval_positions[key][i]  # 从下往上计算位置
                
                # 绘制箭头式区间 |<----->|
                # 1. 左侧竖线
                ax.plot([t_min, t_min], [y_pos - 0.01*y_range, y_pos + 0.01*y_range], 
                       color=region_colors[key], linewidth=2)
                
                # 2. 右侧竖线
                ax.plot([TTB_max_delta, TTB_max_delta], [y_pos - 0.01*y_range, y_pos + 0.01*y_range], 
                       color=region_colors[key], linewidth=2)
                
                # 3. 中间连接线和箭头
                ax.annotate('', xy=(TTB_max_delta, y_pos), xytext=(t_min, y_pos),
                          arrowprops=dict(arrowstyle='<->', color=region_colors[key], linewidth=1.5))
                
                # 4. epsilon标签，带字母上标
                mid_point = (t_min + TTB_max_delta) / 2
                # 计算区间宽度并保留两位小数
                interval_width = round(TTB_max_delta - t_min, 2)
                if letter_tag:
                    eps_label = rf'$\varepsilon^{{{letter_tag}}} = \mathbf{{{interval_width}}}$'
                else:
                    eps_label = rf'$\varepsilon = \mathbf{{{interval_width}}}$'
                ax.text(mid_point, y_pos + 0.01*y_range, eps_label,
                       fontsize=18, ha='center', va='bottom', color=region_colors[key], zorder=100)
                
                # 5. 方案标签
                ax.text(TTB_max_delta + 0.01 * x_range, y_pos, f'{region_short[key]}{idx+1}', 
                      fontsize=12, va='center', color=region_colors[key], fontweight='bold')

    # --- 在最右侧位置标注路径标签，只标注有流量的路径 ---
    path_points = {}  # 存储每个路径ID对应的所有点
    
    # 首先收集每个路径ID对应的所有点
    for row in csv_rows:
        pid = row['path_id']
        t = row['time_cost']
        m = row['money_cost']
        flow = row['flow']
        # 只为有流量的路径添加标签
        if flow > 0:
            if pid not in path_points:
                path_points[pid] = []
            path_points[pid].append((t, m))
    
    if path_points:
        # 找出所有点的最大时间成本用于放置标签
        all_time_costs = [max(points, key=lambda x: x[0])[0] for points in path_points.values()]
        max_time = max(all_time_costs) if all_time_costs else 0
        
        for pid, points in sorted(path_points.items()):
            # 对每个路径ID找到时间成本最大和最小的点
            max_point = max(points, key=lambda x: x[0])
            min_point = min(points, key=lambda x: x[0])
            
            # 获取金钱成本和时间成本
            max_t, m = max_point  # 金钱成本使用最右侧点的值
            min_t, _ = min_point  # 最小时间成本点
            
            # 绘制水平虚线（从最左侧点到标签位置，与标签保持相同高度）
            ax.plot([min_t, max_time + 0.5], [m, m], color='gray', linestyle=':', linewidth=2.5)
            
            # 添加路径标签
            ax.text(max_time + 0.8, m, f'$P_{pid}$', fontsize=20, va='center', ha='left')

    # 完成并保存
    ax.set_xlabel('Time', fontsize=14)
    ax.set_ylabel('Money Cost', fontsize=14)
    ax.grid(True, linestyle='--', alpha=0.7)
    
    # 只使用收集的图例句柄，确保每个区域只有一个图例项
    ax.legend(handles=legend_handles, loc='best', fontsize=30, framealpha=1)
    
    plt.tight_layout()
    plt.savefig(f'results/two_regions_path_costs_zeta{zeta_value}.pdf', format='pdf', dpi=300)
    # plt.show()

    # 保存CSV，只保存有流量的路径
    if csv_rows:
        import pandas as _pd
        df_out = _pd.DataFrame(csv_rows)
        os.makedirs('results', exist_ok=True)
        csv_path = os.path.join('results', f'two_regions_flow_costs_zeta{zeta_value}.csv')
        df_out.to_csv(csv_path, index=False)
        print(f"流量和成本数据已保存至 {csv_path}")

# 主函数
def main():
    # Set LaTeX rendering
    plt.rcParams.update({
        "text.usetex": True,
        "font.family": "serif",
        "font.serif": ["Computer Modern Roman"],
    })
    # Plot comparison for JSON-based data
    for zeta in [8, 16, 24, 32]:
        plot_comparison_path_costs(zeta)
        # 添加调用新的两区域绘图函数
        plot_two_regions_path_costs(zeta)
    
    # Optionally, existing three regions comparison calls can remain or be removed.

if __name__ == "__main__":
    main()