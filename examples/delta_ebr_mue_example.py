#!/usr/bin/env python3
"""
δ-EBR-MUE 求解器示例

该示例展示了如何使用δ-EBR-MUE求解器对不同的网络配置进行求解，并比较不同参数的影响。
"""

import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import time # Added for time.time()

# 添加父目录到模块搜索路径
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.delta_ebr_mue_solver import DeltaEBRMUESolver
from src.ue_solver import UESolver
from src.traffic_network_config import TrafficNetworkConfig

# mac 中文显示
plt.rcParams['font.family'] = ['Arial Unicode MS']


def ensure_output_dir():
    """确保输出目录存在"""
    os.makedirs('output', exist_ok=True)


def solve_multi_od_network():
    """求解多OD对网络的δ-EBR-MUE问题"""
    print("\n=== 求解多OD对网络的δ-EBR-MUE问题 ===")
    config = TrafficNetworkConfig.create_multi_od_network()
    
    # 设置δ和ε参数
    delta = 20
    epsilon = (0.5, 0.5)  # (时间容忍度, 金钱容忍度)
    
    # 初始化求解器
    solver = DeltaEBRMUESolver(config, delta=delta, epsilon=epsilon)
    
    # 运行求解
    result = solver.run(max_iter=20, convergence_threshold=0.01, use_epsilon=True)
    
    # 显示结果
    solver.display_results(result)
    solver.display_convergence()
    solver.plot_cost_analysis(result)
    
    return solver, result


def compare_with_ue(network_type='multi_od'):
    """比较δ-EBR-MUE和传统UE的区别"""
    print("\n=== 比较δ-EBR-MUE和传统UE ===")
    
    if network_type == 'multi_od':
        config = TrafficNetworkConfig.create_multi_od_network()
    elif network_type == 'single_od':
        config = TrafficNetworkConfig.create_single_od_network()
    else:
        config = TrafficNetworkConfig.create_basic_network()
    
    # 求解UE
    print("\n求解传统UE...")
    ue_solver = UESolver(config)
    ue_solver.run()
    
    # 获取UE结果
    ue_flows = {i: ue_solver.model.flow[i].value for i in range(1, config.num_paths + 1)}
    ue_times = {i: ue_solver.model.path_time[i].value for i in range(1, config.num_paths + 1)}
    
    # 计算UE的金钱成本
    ue_money = {}
    for i in range(1, config.num_paths + 1):
        ue_money[i] = sum(
            config.path_link_matrix.get((i, j), 0) * config.link_money_cost[j]
            for j in range(1, config.num_links + 1)
        )
    
    # 求解δ-EBR-MUE，尝试不同的δ值
    deltas = [0.5, 2.0, 5.0]
    delta_results = {}
    
    for delta in deltas:
        print(f"\n求解δ-EBR-MUE, δ = {delta}...")
        delta_solver = DeltaEBRMUESolver(config, delta=delta, epsilon=(0.5, 0.5))
        result = delta_solver.run(max_iter=20, convergence_threshold=0.01, use_epsilon=True)
        delta_results[delta] = result
    
    # 绘制比较图
    plt.figure(figsize=(12, 8))
    
    # 对比使用的路径数量
    ax1 = plt.subplot(211)
    models = ['UE'] + [f'δ-EBR-MUE (δ={d})' for d in deltas]
    
    used_paths_count = []
    # UE中使用的路径数（流量 > 0.01）
    ue_used_paths = sum(1 for flow in ue_flows.values() if flow > 0.01)
    used_paths_count.append(ue_used_paths)
    
    # δ-EBR-MUE中使用的路径数
    for delta, result in delta_results.items():
        delta_flows = result['flows']
        delta_used_paths = sum(1 for flow in delta_flows.values() if flow > 0.01)
        used_paths_count.append(delta_used_paths)
    
    ax1.bar(models, used_paths_count)
    ax1.set_title('使用的路径数量比较')
    ax1.set_ylabel('使用的路径数')
    
    # 对比平均时间和金钱成本
    ax2 = plt.subplot(212)
    
    # 计算UE的平均时间和金钱成本
    total_flow_ue = sum(ue_flows.values())
    avg_time_ue = sum(ue_flows[i] * ue_times[i] for i in ue_flows) / total_flow_ue
    avg_money_ue = sum(ue_flows[i] * ue_money[i] for i in ue_flows) / total_flow_ue
    
    avg_times = [avg_time_ue]
    avg_money = [avg_money_ue]
    
    # 计算各δ-EBR-MUE模型的平均时间和金钱成本
    for delta, result in delta_results.items():
        delta_flows = result['flows']
        delta_times = result['time_costs']
        delta_money = result['money_costs']
        
        total_flow_delta = sum(delta_flows.values())
        avg_time_delta = sum(delta_flows[i] * delta_times[i] for i in delta_flows if i in delta_times) / total_flow_delta
        avg_money_delta = sum(delta_flows[i] * delta_money[i] for i in delta_flows if i in delta_money) / total_flow_delta
        
        avg_times.append(avg_time_delta)
        avg_money.append(avg_money_delta)
    
    x = np.arange(len(models))
    width = 0.35
    
    ax2.bar(x - width/2, avg_times, width, label='平均时间成本')
    ax2.bar(x + width/2, avg_money, width, label='平均金钱成本')
    ax2.set_title('平均成本比较')
    ax2.set_xlabel('模型')
    ax2.set_ylabel('成本值')
    ax2.set_xticks(x)
    ax2.set_xticklabels(models, rotation=15)
    ax2.legend()
    
    plt.tight_layout()
    plt.savefig('output/model_comparison.png')
    plt.close()
    
    print("\n模型比较图已保存至: output/model_comparison.png")
    
    # 绘制收敛比较图
    plt.figure(figsize=(10, 6))
    
    for delta, result in delta_results.items():
        # 从结果列表中提取相对间隙
        gaps = []
        for iteration_result in result['results']:
            gaps.append(iteration_result['relative_gap'])
        
        iterations = list(range(1, len(gaps) + 1))
        plt.semilogy(iterations, gaps, 'o-', linewidth=2, label=f'δ = {delta}')
    
    plt.xlabel('迭代次数')
    plt.ylabel('相对间隙 (对数刻度)')
    plt.title('不同δ值下的收敛过程对比')
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.legend()
    
    plt.tight_layout()
    plt.savefig('output/convergence_comparison.png')
    plt.close()
    
    print("收敛对比图已保存至: output/convergence_comparison.png")
    
    return ue_solver, delta_results


def explore_epsilon_effect():
    """探索ε参数对结果的影响"""
    print("\n=== 探索ε参数对结果的影响 ===")
    config = TrafficNetworkConfig.create_basic_network()
    
    # 固定δ值，尝试不同的ε值
    delta = 20
    epsilons = [(0.0, 0.0), (0.5, 0.5), (1.0, 1.0), (2.0, 2.0)]
    
    epsilon_results = {}
    
    for eps in epsilons:
        print(f"\n求解δ-EBR-MUE, ε = {eps}...")
        solver = DeltaEBRMUESolver(config, delta=delta, epsilon=eps)
        result = solver.run(max_iter=20, convergence_threshold=0.01, use_epsilon=True)
        epsilon_results[eps] = result
    
    # TODO: 添加ε敏感性分析的可视化
    
    return epsilon_results


def validate_large_single_od_network():
    """验证大型单OD对网络的δ-EBR-MUE求解效果"""
    print("\n=== 验证大型单OD对网络的δ-EBR-MUE求解效果 ===")
    
    # 创建大型单OD对网络配置
    config = TrafficNetworkConfig.create_large_single_od_network()
    
    # 配置参数，根据网络特性调整delta和epsilon
    delta = 20  # 考虑到free_flow_time范围是10-28，设置一个适中的delta值
    epsilon = (0, 0)  # 设置较大的epsilon，使得更多路径成为候选
    
    print(f"\n使用参数: delta={delta}, epsilon={epsilon}")
    
    # 比较两种不同的求解方法
    solver_methods = ['convex', 'vi']
    solver_results = {}
    
    for method in solver_methods:
        print(f"\n=== 使用 {method} 方法求解 ===")
        
        # 初始化求解器
        solver = DeltaEBRMUESolver(config, delta=delta, epsilon=epsilon)
        
        # 运行求解
        print(f"\n运行求解器...")
        start_time = time.time()
        result = solver.run(max_iter=30, convergence_threshold=0.01, use_epsilon=True, solver_method=method)
        solve_time = time.time() - start_time
        
        # 显示结果
        solver.display_results(result)
        solver.plot_cost_analysis(result, fig_name=f'{method}_cost_analysis.png')
        
        # 存储结果
        solver_results[method] = {
            'result': result,
            'solver': solver,
            'solve_time': solve_time
        }
    
    # 比较两种方法的结果
    print("\n=== 两种求解方法的对比 ===")
    
    # 比较求解时间
    convex_time = solver_results['convex']['solve_time']
    vi_time = solver_results['vi']['solve_time']
    print(f"凸优化方法求解时间: {convex_time:.2f} 秒")
    print(f"变分不等式方法求解时间: {vi_time:.2f} 秒")
    print(f"VI方法相对凸优化方法的求解时间比: {vi_time/convex_time:.2f}")
    
    # 比较有效路径数量
    convex_result = solver_results['convex']['result']
    vi_result = solver_results['vi']['result']
    
    convex_feasible_paths = sum(len(paths) for paths in convex_result['feasible_paths'].values())
    vi_feasible_paths = sum(len(paths) for paths in vi_result['feasible_paths'].values())
    
    convex_used_paths = sum(1 for flow in convex_result['flows'].values() if flow > 0.01)
    vi_used_paths = sum(1 for flow in vi_result['flows'].values() if flow > 0.01)
    
    print(f"\n凸优化方法可行路径数: {convex_feasible_paths}")
    print(f"VI方法可行路径数: {vi_feasible_paths}")
    print(f"\n凸优化方法实际使用路径数: {convex_used_paths}")
    print(f"VI方法实际使用路径数: {vi_used_paths}")
    
    # 比较平均成本
    total_paths = config.num_paths
    
    # 计算凸优化方法的平均成本
    convex_flows = convex_result['flows']
    convex_time_costs = convex_result['time_costs']
    convex_money_costs = convex_result['money_costs']
    
    convex_total_flow = sum(convex_flows.values())
    convex_avg_time = sum(convex_flows[i] * convex_time_costs[i] 
                         for i in convex_flows if i in convex_time_costs) / convex_total_flow if convex_total_flow > 0 else 0
    convex_avg_money = sum(convex_flows[i] * convex_money_costs[i] 
                          for i in convex_flows if i in convex_money_costs) / convex_total_flow if convex_total_flow > 0 else 0
    
    # 计算VI方法的平均成本
    vi_flows = vi_result['flows']
    vi_time_costs = vi_result['time_costs']
    vi_money_costs = vi_result['money_costs']
    
    vi_total_flow = sum(vi_flows.values())
    vi_avg_time = sum(vi_flows[i] * vi_time_costs[i] 
                     for i in vi_flows if i in vi_time_costs) / vi_total_flow if vi_total_flow > 0 else 0
    vi_avg_money = sum(vi_flows[i] * vi_money_costs[i] 
                      for i in vi_flows if i in vi_money_costs) / vi_total_flow if vi_total_flow > 0 else 0
    
    print(f"\n凸优化方法平均时间成本: {convex_avg_time:.2f}")
    print(f"VI方法平均时间成本: {vi_avg_time:.2f}")
    print(f"凸优化方法平均金钱成本: {convex_avg_money:.2f}")
    print(f"VI方法平均金钱成本: {vi_avg_money:.2f}")
    
    # 绘制两种方法的收敛过程对比图
    plt.figure(figsize=(12, 6))
    
    # 提取收敛数据
    convex_gaps = [r['relative_gap'] for r in convex_result['results']]
    vi_gaps = [r['relative_gap'] for r in vi_result['results']]
    
    convex_iterations = range(1, len(convex_gaps) + 1)
    vi_iterations = range(1, len(vi_gaps) + 1)
    
    plt.semilogy(convex_iterations, convex_gaps, 'o-', linewidth=2, label='凸优化方法')
    plt.semilogy(vi_iterations, vi_gaps, 's-', linewidth=2, label='VI方法')
    
    plt.xlabel('迭代次数')
    plt.ylabel('相对间隙 (对数刻度)')
    plt.title('不同求解方法的收敛过程对比')
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.legend()
    plt.tight_layout()
    
    plt.savefig('output/solver_methods_comparison.png')
    plt.close()
    
    print("\n收敛过程对比图已保存至: output/solver_methods_comparison.png")
    
    # 绘制路径流量分布对比图
    plt.figure(figsize=(12, 8))
    
    # 选择流量最高的前20条路径进行可视化
    path_indices = sorted(range(1, total_paths + 1), 
                        key=lambda x: max(convex_flows.get(x, 0), vi_flows.get(x, 0)), 
                        reverse=True)[:20]
    
    x = np.arange(len(path_indices))
    width = 0.35
    
    convex_values = [convex_flows.get(i, 0) for i in path_indices]
    vi_values = [vi_flows.get(i, 0) for i in path_indices]
    
    plt.bar(x - width/2, convex_values, width, label='凸优化方法')
    plt.bar(x + width/2, vi_values, width, label='VI方法')
    
    plt.xlabel('路径ID')
    plt.ylabel('流量')
    plt.title('两种方法的路径流量分布对比(流量最高的20条路径)')
    plt.xticks(x, [str(i) for i in path_indices])
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.tight_layout()
    
    plt.savefig('output/flow_distribution_comparison.png')
    plt.close()
    
    print("流量分布对比图已保存至: output/flow_distribution_comparison.png")
    
    return solver_results


def main():
    """主函数"""
    # 确保输出目录存在
    ensure_output_dir()
    
    # 根据命令行参数选择要运行的示例
    if len(sys.argv) > 1:
        example = sys.argv[1]
        if example == "multi":
            solve_multi_od_network()
        elif example == "compare":
            compare_with_ue()
        elif example == "epsilon":
            explore_epsilon_effect()
        elif example == "large":
            validate_large_single_od_network()
        else:
            print(f"未知的示例: {example}")
            print("可用示例: multi, compare, epsilon, large")
    else:
        # 默认运行指定示例
        print("运行大型单OD对网络的δ-EBR-MUE求解效果")
        validate_large_single_od_network()


if __name__ == "__main__":
    main() 