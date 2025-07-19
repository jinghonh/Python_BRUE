#!/usr/bin/env python3
"""
用户均衡(UE)求解器示例

该示例展示了如何使用UE求解器对不同的网络配置进行求解，并可视化结果。
"""

import os
import sys
import matplotlib.pyplot as plt

# 添加父目录到模块搜索路径
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.ue_solver import UESolver
from src.traffic_network_config import TrafficNetworkConfig

# mac 中文显示
plt.rcParams['font.family'] = ['Arial Unicode MS']

def ensure_output_dir():
    """确保输出目录存在"""
    os.makedirs('output', exist_ok=True)


def solve_basic_network():
    """求解基础网络"""
    print("\n=== 求解基础网络 ===")
    config = TrafficNetworkConfig.create_basic_network()
    solver = UESolver(config)
    plots = solver.run()
    return plots


def solve_single_od_network():
    """求解单一OD对网络"""
    print("\n=== 求解单一OD对网络 ===")
    config = TrafficNetworkConfig.create_single_od_network()
    solver = UESolver(config)
    plots = solver.run()
    return plots


def solve_multi_od_network():
    """求解多OD对网络"""
    print("\n=== 求解多OD对网络 ===")
    config = TrafficNetworkConfig.create_multi_od_network()
    solver = UESolver(config)
    plots = solver.run()
    return plots


def solve_two_od_network():
    """求解两OD对网络"""
    print("\n=== 求解两OD对网络 ===")
    config = TrafficNetworkConfig.create_two_od_network()
    solver = UESolver(config)
    plots = solver.run()
    return plots


def solve_anti_example_network():
    """求解反例网络"""
    print("\n=== 求解反例网络 ===")
    config = TrafficNetworkConfig.create_anti_example_network()
    solver = UESolver(config)
    plots = solver.run()
    return plots


def compare_networks():
    """比较不同网络的求解结果"""
    print("\n=== 比较不同网络的求解结果 ===")
    
    # 创建不同的网络配置
    configs = {
        "基础网络": TrafficNetworkConfig.create_basic_network(),
        "单一OD对网络": TrafficNetworkConfig.create_single_od_network(),
        "多OD对网络": TrafficNetworkConfig.create_multi_od_network()
    }
    
    # 存储每个网络的平均时间成本
    avg_times = {}
    
    plt.figure(figsize=(12, 8))
    
    for idx, (name, config) in enumerate(configs.items(), 1):
        print(f"\n求解 {name}...")
        solver = UESolver(config)
        solver.run()
        
        # 计算平均时间成本
        total_flow = 0
        total_weighted_time = 0
        
        for group_name, group_pairs in config.od_groups.items():
            for i in group_pairs:
                flow = solver.model.flow[i].value
                time = solver.model.path_time[i].value
                total_flow += flow
                total_weighted_time += flow * time
        
        avg_time = total_weighted_time / total_flow if total_flow > 0 else 0
        avg_times[name] = avg_time
        
        # 绘制每个网络的路径流量分布
        plt.subplot(len(configs), 1, idx)
        
        all_paths = []
        all_flows = []
        
        for group_name, group_pairs in config.od_groups.items():
            paths = list(group_pairs)
            flows = [solver.model.flow[i].value for i in paths]
            all_paths.extend(paths)
            all_flows.extend(flows)
        
        plt.bar(all_paths, all_flows)
        plt.xlabel('路径ID')
        plt.ylabel('流量')
        plt.title(f'{name} (平均时间成本: {avg_time:.2f})')
        plt.grid(True, linestyle='--', alpha=0.6)
    
    plt.tight_layout()
    plt.savefig('output/network_comparison.png')
    plt.close()
    
    # 绘制平均时间成本比较图
    plt.figure(figsize=(10, 6))
    plt.bar(avg_times.keys(), avg_times.values())
    plt.ylabel('平均时间成本')
    plt.title('不同网络的平均时间成本比较')
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.tight_layout()
    plt.savefig('output/avg_time_comparison.png')
    plt.close()
    
    print("\n各网络平均时间成本:")
    for name, avg_time in avg_times.items():
        print(f"{name}: {avg_time:.2f}")
    
    return {
        "network_comparison": "output/network_comparison.png",
        "avg_time_comparison": "output/avg_time_comparison.png"
    }


def main():
    """主函数"""
    # 确保输出目录存在
    ensure_output_dir()
    
    # 根据命令行参数选择要运行的示例
    if len(sys.argv) > 1:
        example = sys.argv[1]
        if example == "basic":
            solve_basic_network()
        elif example == "single":
            solve_single_od_network()
        elif example == "multi":
            solve_multi_od_network()
        elif example == "two":
            solve_two_od_network()
        elif example == "anti":
            solve_anti_example_network()
        elif example == "compare":
            compare_networks()
        else:
            print(f"未知的示例: {example}")
            print("可用示例: basic, single, multi, two, anti, compare")
    else:
        # 默认运行所有示例
        solve_multi_od_network()
        compare_networks()


if __name__ == "__main__":
    main() 