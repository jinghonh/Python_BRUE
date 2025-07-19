"""
基本BRUE求解示例
展示如何创建交通网络并运行求解
"""
from src.brue_solver import BRUESolver
from src.traffic_network_config import TrafficNetworkConfig

def main():
    # 创建网络配置
    config = TrafficNetworkConfig.create_basic_network()
    
    # 初始化求解器
    solver = BRUESolver(config)
    
    # 运行迭代求解
    results = solver.run_with_iterations()
    
    # 显示结果
    solver.display_results()
    
    # 分析路径成本
    solver.analyze_path_costs()
    
    # 绘制初始成本分析
    solver.plot_initial_costs()
    
    return results

if __name__ == "__main__":
    main() 