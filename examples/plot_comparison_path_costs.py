import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from scipy.io import loadmat  # 添加导入.mat文件的库

# 读取CSV数据
def load_data(file_path):
    return pd.read_csv(file_path)

# 读取MATLAB的.mat文件
def load_mat_data(file_path):
    return loadmat(file_path)

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

# 从.mat文件计算T_max上限
def calculate_tmax_from_mat(mat_file_path):
    """
    从.mat文件中读取数据并计算T_max上限
    
    参数:
        mat_file_path: .mat文件路径
        
    返回:
        upper_limit_x: T_max上限的x坐标
        upper_limit_y: T_max上限的y坐标
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
        
        # 计算T_max上限（使用边界中点）
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
        print(f"从.mat文件计算T_max时出错: {e}")
        return np.array([]), np.array([])

# 绘制三区域路径成本对比图
def plot_three_regions_comparison(df, zeta_value, figsize=(10, 8), use_mat_file=False):
    # 过滤指定zeta值的数据
    df_zeta = df[df['zeta值'] == zeta_value]
    
    # 创建图形
    fig, ax = plt.subplots(figsize=figsize)
    
    # 设置颜色和标记样式
    colors = {
        1: [0.0, 0.6, 0.3],  # 绿色 - 区域1：只满足路径约束
        2: [0.85, 0.4, 0.0], # 橙色 - 区域2：满足全部约束但违反T_max
        3: [0.0, 0.4, 0.8]   # 蓝色 - 区域3：满足全部约束且满足T_max
    }
    markers = {
        1: 'o',  # 区域1：圆形
        2: 's',  # 区域2：方形
        3: 'd'   # 区域3：菱形
    }
    
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
                scatter = ax.scatter(time_costs, money_costs, c=[colors[region]], 
                                     marker=markers[region], s=60, alpha=0.8)
                
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
                    label = f"$RS_0^{{\\varepsilon}}$ {scheme}"
                
                # 添加到图例句柄
                legend_handles.append(mpatches.Patch(color=colors[region], label=label))
    
    # 绘制 T_max 上限线
    if use_mat_file:
        # 从.mat文件计算T_max上限
        mat_file_path = f'matlab/cache/cache_zeta{zeta_value}_subset0.mat'
        upper_limit_x, upper_limit_y = calculate_tmax_from_mat(mat_file_path)
        
        if len(upper_limit_x) > 0 and len(upper_limit_y) > 0:
            # 排序以确保正确连线
            sort_idx = np.argsort(upper_limit_y)
            upper_limit_x = upper_limit_x[sort_idx]
            upper_limit_y = upper_limit_y[sort_idx]
            
            # 绘制T_max线
            ax.plot(upper_limit_x, upper_limit_y, '-.', color=[0.8, 0.2, 0.2], linewidth=2)
        else:
            print(f"警告: 无法从.mat文件计算T_max上限")
    else:
        # 使用CSV数据中的T_max上限
        tmax_data = df_zeta[df_zeta['区域'] == 4]  # 使用区域4表示 T_max 数据
        if not tmax_data.empty:
            # 按金钱成本排序
            tmax_data = tmax_data.sort_values('金钱成本')
            # 绘制实际 T_max 线
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
    ax.set_xlabel('Time Cost', fontsize=12)
    ax.set_ylabel('Money Cost', fontsize=12)
    ax.grid(True)
    
    # 扩展x轴范围以容纳路径标签
    ax.set_xlim(x_min - margin, x_max + margin * 5)
    
    # 添加图例
    handles = legend_handles + [Line2D([0], [0], linestyle='-.', color=[0.8, 0.2, 0.2], label='$T_{max}$')]
    ax.legend(handles=handles, loc='best', fontsize=10, framealpha=1)
    
    plt.title(f'Path Cost Comparison of Three Regions ($\\varepsilon = {zeta_value}$)', fontsize=14)
    plt.tight_layout()
    plt.savefig(f'three_regions_comparison_zeta{zeta_value}.pdf', format='pdf', dpi=300)
    plt.show()

# 绘制路径成本对比图（区域1和区域3）
def plot_comparison_path_costs(df, zeta_value, figsize=(10, 8)):
    # 过滤指定zeta值的数据
    df_zeta = df[df['zeta值'] == zeta_value]
    
    # 只保留区域1和区域3
    df_filtered = df_zeta[(df_zeta['区域'] == 1) | (df_zeta['区域'] == 3)]
    
    # 创建图形
    fig, ax = plt.subplots(figsize=figsize)
    
    # 设置颜色和标记
    colors = {
        1: [0.0, 0.6, 0.3],  # 绿色 - 区域1：只满足路径约束
        3: [0.0, 0.4, 0.8]   # 蓝色 - 区域3：满足全部约束且满足T_max
    }
    markers = {
        1: 'o',  # 区域1：圆形
        3: 'd'   # 区域3：菱形
    }
    
    legend_handles = []
    
    # 计算图表的x轴范围，用于后续添加标签
    x_min, x_max = float('inf'), float('-inf')
    
    # 获取所有路径编号
    unique_paths = df_filtered['路径编号'].unique()
    
    # 处理每个区域的数据
    for region in [1, 3]:
        for scheme in [1, 2]:
            # 获取当前方案的数据
            scheme_data = df_filtered[(df_filtered['区域'] == region) & (df_filtered['方案编号'] == scheme)]
            
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
                scatter = ax.scatter(time_costs, money_costs, c=[colors[region]], 
                                     marker=markers[region], s=60, alpha=0.8)
                
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
                else:
                    label = f"$BS_0^{{\\varepsilon}}$ {scheme}"
                
                # 添加到图例句柄
                legend_handles.append(mpatches.Patch(color=colors[region], label=label))
    
    # 在右侧添加路径标签
    margin = (x_max - x_min) * 0.05  # 设置一个边距
    # 获取唯一的金钱成本值并排序
    unique_money_costs = sorted(df_filtered['金钱成本'].unique(), reverse=True)
    
    for money_cost in unique_money_costs:
        # 对于每个金钱成本，找到对应的路径编号
        path_data = df_filtered[df_filtered['金钱成本'] == money_cost]
        if not path_data.empty:
            path_id = path_data['路径编号'].values[0]
            # 在右侧添加路径标签
            ax.text(x_max + margin, money_cost, f'$P_{path_id}$', 
                    fontsize=10, va='center', ha='left', color='#4D4D4D')
    
    # 设置轴标签和图例
    ax.set_xlabel('Time Cost', fontsize=12)
    ax.set_ylabel('Money Cost', fontsize=12)
    ax.grid(True)
    
    # 扩展x轴范围以容纳路径标签
    ax.set_xlim(x_min - margin, x_max + margin * 5)
    
    # 添加图例
    ax.legend(handles=legend_handles, loc='best', fontsize=10, framealpha=1)
    
    plt.title(f'Path Cost Comparison ($\\varepsilon = {zeta_value}$)', fontsize=14)
    plt.tight_layout()
    plt.savefig(f'path_costs_comparison_zeta{zeta_value}.pdf', format='pdf', dpi=300)
    plt.show()

# 主函数
def main():
    # 读取数据
    df = load_data('matlab/results/result.csv')
    
    # 设置LaTeX渲染
    plt.rcParams.update({
        "text.usetex": True,
        "font.family": "serif",
        "font.serif": ["Computer Modern Roman"],
    })
    
    # 绘制四张图
    plot_comparison_path_costs(df, 15)     # zeta=15的路径成本对比图
    plot_comparison_path_costs(df, 31)     # zeta=31的路径成本对比图

    plot_three_regions_comparison(df, 15, use_mat_file=True)  # 使用.mat文件的zeta=15三区域对比图
    plot_three_regions_comparison(df, 31, use_mat_file=True)  # 使用.mat文件的zeta=31三区域对比图

if __name__ == "__main__":
    main()