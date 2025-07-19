"""
BRUE (Boundedly Rational User Equilibrium) 交通网络求解器
主入口点脚本
"""

import argparse
import sys
from src.brue_solver import BRUESolver
from src.brue_set_solver import BRUESetSolver
from src.brue_matlab_solver import BRUEMatlabSolver
from src.traffic_network_config import TrafficNetworkConfig

def run_basic_solver(args):
    """运行基础求解器"""
    # 创建网络配置
    if args.network_type == 'basic':
        config = TrafficNetworkConfig.create_basic_network()
    elif args.network_type == 'single':
        config = TrafficNetworkConfig.create_single_od_network()
    elif args.network_type == 'multi':
        config = TrafficNetworkConfig.create_multi_od_network()
    elif args.network_type == 'two':
        config = TrafficNetworkConfig.create_two_od_network()
    else:
        config = TrafficNetworkConfig.create_anti_example_network()

    # 初始化求解器
    solver = BRUESolver(config)
    
    # 运行求解
    results = solver.run_with_iterations()
    
    # 显示结果
    solver.display_results()
    
    # 分析路径成本
    solver.analyze_path_costs()
    
    # 可视化
    if args.visualize:
        solver.plot_initial_costs()
    
    return results

def run_set_solver(args):
    """运行集合求解器"""
    # 创建网络配置
    if args.network_type == 'basic':
        config = TrafficNetworkConfig.create_basic_network()
    elif args.network_type == 'single':
        config = TrafficNetworkConfig.create_single_od_network()
    elif args.network_type == 'multi':
        config = TrafficNetworkConfig.create_multi_od_network()
    elif args.network_type == 'two':
        config = TrafficNetworkConfig.create_two_od_network()
    else:
        config = TrafficNetworkConfig.create_anti_example_network()

    # 初始化求解器
    solver = BRUESetSolver(config)
    
    # 求解BRUE集合
    solutions = solver.find_brue_set_grid(args.zeta, args.num_points)
    
    # 可视化
    if args.visualize:
        solver.visualize_brue_set(solutions, args.zeta)
    
    return solutions

def main():
    parser = argparse.ArgumentParser(description='BRUE 交通网络求解器')
    subparsers = parser.add_subparsers(help='命令', dest='command')
    
    # 基础求解器命令
    basic_parser = subparsers.add_parser('basic', help='运行基础求解器')
    basic_parser.add_argument('--network', dest='network_type', choices=['basic', 'single', 'multi', 'two', 'anti'], 
                            default='basic', help='网络类型')
    basic_parser.add_argument('--visualize', action='store_true', help='显示可视化图表')
    
    # 集合求解器命令
    set_parser = subparsers.add_parser('set', help='运行集合求解器')
    set_parser.add_argument('--network', dest='network_type', choices=['basic', 'single', 'multi', 'two', 'anti'], 
                           default='basic', help='网络类型')
    set_parser.add_argument('--zeta', type=float, default=1.0, help='容忍度参数')
    set_parser.add_argument('--num-points', type=int, default=50, help='每个维度的采样点数')
    set_parser.add_argument('--visualize', action='store_true', help='显示可视化图表')
    
    args = parser.parse_args()
    
    if args.command == 'basic':
        run_basic_solver(args)
    elif args.command == 'set':
        run_set_solver(args)
    else:
        parser.print_help()
        sys.exit(1)

if __name__ == "__main__":
    main() 