from .brue_solver import BRUESolver
from .traffic_network_config import TrafficNetworkConfig
from pyomo.environ import *
import matplotlib.pyplot as plt
import numpy as np
import itertools
from rich.console import Console
from rich.table import Table
from typing import Dict, List, Tuple, Any
from tqdm import tqdm

class BRUESetSolver(BRUESolver):
    def __init__(self, config: TrafficNetworkConfig):
        super().__init__(config)
        self.console = Console()
        # 将路径-链接关系转换为矩阵，便于数值计算
        self.path_link_matrix = config.get_path_link_array()
        
    def calculate_real_time(self, flows: np.ndarray) -> np.ndarray:
        """计算实际旅行时间
        Args:
            flows: 路径流量数组
        Returns:
            实际旅行时间数组
        """
        # 确保path_link_matrix是正确的形状
        if self.path_link_matrix.shape[0] != len(flows):
            raise ValueError(f"路径数量不匹配：flows长度为{len(flows)}，但path_link_matrix行数为{self.path_link_matrix.shape[0]}")
        
        # 计算每条路径的旅行时间
        path_flows = flows @ self.path_link_matrix
        free_flow_times = np.array([
            self.config.free_flow_time[i]
            for i in range(1, self.config.num_links + 1)
        ])
        capacities = np.array([
            self.config.link_capacity[i]
            for i in range(1, self.config.num_links + 1)
        ])
        
        # BPR函数计算时间
        link_times = free_flow_times * (1 + 0.15 * (path_flows / capacities) ** 4)
        path_times = link_times @ self.path_link_matrix.T
        
        # 添加感知误差项
        return path_times + 15 * (1 - np.exp(-0.02 * path_times))
    
    def calculate_money_costs(self, flows: np.ndarray) -> np.ndarray:
        """计算金钱成本
        Args:
            flows: 路径流量数组
        Returns:
            金钱成本数组
        """
        money_costs = np.array([
            self.config.link_money_cost[i]
            for i in range(1, self.config.num_links + 1)
        ])
        # 计算每条路径的金钱成本
        path_money_costs = []
        for i in range(len(flows)):
            path_cost = 0
            for j in range(len(money_costs)):
                path_cost += self.path_link_matrix[i, j] * money_costs[j]
            path_money_costs.append(path_cost)
        return np.array(path_money_costs)
        
    def find_brue_set_grid(self, zeta: float, num_points: int = 50) -> List[Dict[str, Any]]:
        """使用网格搜索方法寻找BRUE集合
        Args:
            zeta: 容忍度参数
            num_points: 每个维度的采样点数
        Returns:
            满足条件的解集合
        """
        # 获取路径数量
        num_paths = self.config.num_links  # 直接使用配置中的路径数量
        total_demand = sum(self.config.od_demands.values())
        
        if num_paths > 6:
            self.console.print("[red]当前最多支持6条路径")
            return []
        
        # 使用Dirichlet分布生成满足需求约束的候选解
        num_samples = num_points ** 2  # 增加采样点数以获得更好的覆盖
        alpha = np.ones(num_paths)  # Dirichlet分布的参数
        candidates = np.random.dirichlet(alpha, num_samples) * total_demand
        
        # 验证每个候选解
        valid_solutions = []
        for flows in tqdm(candidates, desc="验证候选解"):
            # 计算时间和金钱成本
            time_costs = self.calculate_real_time(flows)
            money_costs = self.calculate_money_costs(flows)
            
            # 验证BRUE条件
            is_valid = True
            for i, j in itertools.combinations(range(num_paths), 2):
                time_diff = time_costs[i] - time_costs[j]
                money_diff = money_costs[i] - money_costs[j]
                
                # 检查时间成本差
                if abs(time_diff) > zeta:
                    is_valid = False
                    break
                    
                # 检查成本差乘积
                if time_diff * money_diff >= -1e-6:
                    is_valid = False
                    break
            
            if is_valid:
                solution = {
                    'flows': {i+1: flow for i, flow in enumerate(flows)},
                    'time_costs': {i+1: cost for i, cost in enumerate(time_costs)},
                    'money_costs': {i+1: cost for i, cost in enumerate(money_costs)}
                }
                valid_solutions.append(solution)
                
        self.console.print(f"[green]找到 {len(valid_solutions)} 个有效解")
        return valid_solutions
    
    def visualize_brue_set(self, solutions: List[Dict], zeta: float):
        """可视化BRUE集合的解"""
        if not solutions:
            self.console.print("[red]没有找到有效解")
            return
            
        num_paths = len(solutions[0]['flows'])
        
        # 1. 流量空间可视化（使用散点矩阵）
        if num_paths > 2:
            fig, axes = plt.subplots(num_paths, num_paths, figsize=(15, 15))
            plt.suptitle(f'BRUE Set Solutions in Flow Space (ζ = {zeta})')
            
            for i in range(num_paths):
                for j in range(num_paths):
                    if i != j:
                        flows_i = [sol['flows'][i+1] for sol in solutions]
                        flows_j = [sol['flows'][j+1] for sol in solutions]
                        axes[i, j].scatter(flows_i, flows_j, alpha=0.6, c='blue', s=5)
                        axes[i, j].set_xlabel(f'Path {i+1}')
                        axes[i, j].set_ylabel(f'Path {j+1}')
                    else:
                        axes[i, j].text(0.5, 0.5, f'Path {i+1}', 
                                      ha='center', va='center')
                        axes[i, j].axis('off')
        else:
            # 原有的两条路径可视化代码保持不变
            plt.figure(figsize=(12, 8))
            flows_x = [sol['flows'][1] for sol in solutions]
            flows_y = [sol['flows'][2] for sol in solutions]
            
            plt.scatter(flows_x, flows_y, alpha=0.6, c='blue', s=5)
            plt.title(f'BRUE Set Solutions in Flow Space (ζ = {zeta})')
            plt.xlabel('Flow on Path 1')
            plt.ylabel('Flow on Path 2')
            plt.grid(True, alpha=0.3)
            
            # 添加总需求约束线
            demand = sum(self.config.od_demands.values())
            x = np.linspace(0, demand, 100)
            plt.plot(x, demand - x, 'r--', alpha=0.5, label='Demand Constraint')
            plt.legend()
        
        plt.show()
        
        # 2. 成本空间可视化
        plt.figure(figsize=(12, 8))
        colors = plt.cm.rainbow(np.linspace(0, 1, num_paths))
        
        for sol in solutions:
            time_costs = list(sol['time_costs'].values())
            money_costs = list(sol['money_costs'].values())
            
            # 绘制每条路径的点
            for i in range(num_paths):
                plt.scatter(time_costs[i], money_costs[i], 
                           c=[colors[i]], alpha=0.6, s=5,
                           label=f'Path {i+1}' if sol == solutions[0] else "")
        
        plt.title(f'Cost Space (ζ = {zeta})')
        plt.xlabel('Travel Time')
        plt.ylabel('Money Cost')
        plt.grid(True, alpha=0.3)
        plt.legend()
        plt.show()
        
        # 3. 显示统计信息
        self.display_solution_stats(solutions)
    
    def display_solution_stats(self, solutions: List[Dict]):
        """显示解的统计信息"""
        if not solutions:
            self.console.print("[red]没有有效解可供分析")
            return
        
        num_paths = len(solutions[0]['flows'])
        table = Table(title="解集合统计信息")
        table.add_column("统计项", style="cyan")
        
        # 为每条路径添加列
        for i in range(1, num_paths + 1):
            table.add_column(f"路径{i}", justify="right")
        
        # 计算统计量
        flows = {i: [sol['flows'][i] for sol in solutions] for i in range(1, num_paths + 1)}
        time_costs = {i: [sol['time_costs'][i] for sol in solutions] for i in range(1, num_paths + 1)}
        money_costs = {i: [sol['money_costs'][i] for sol in solutions] for i in range(1, num_paths + 1)}
        
        # 添加统计行
        stats_rows = [
            ("流量均值", *[f"{np.mean(flows[i]):.2f}" for i in range(1, num_paths + 1)]),
            ("流量标准差", *[f"{np.std(flows[i]):.2f}" for i in range(1, num_paths + 1)]),
            ("时间成本均值", *[f"{np.mean(time_costs[i]):.2f}" for i in range(1, num_paths + 1)]),
            ("时间成本标准差", *[f"{np.std(time_costs[i]):.2f}" for i in range(1, num_paths + 1)]),
            ("金钱成本均值", *[f"{np.mean(money_costs[i]):.2f}" for i in range(1, num_paths + 1)]),
            ("金钱成本标准差", *[f"{np.std(money_costs[i]):.2f}" for i in range(1, num_paths + 1)])
        ]
        
        for row in stats_rows:
            table.add_row(*row)
        
        self.console.print(table)

def main():
    # 测试代码
    config = TrafficNetworkConfig.create_basic_network()
    solver = BRUESetSolver(config)
    
    # 使用不同的zeta值测试
    zeta_values = [15, 31]
    for zeta in zeta_values:
        solver.console.print(f"\n[yellow]测试 ζ = {zeta}")
        solutions = solver.find_brue_set_grid(zeta, num_points=50)
        solver.visualize_brue_set(solutions, zeta)
        
if __name__ == "__main__":
    main() 