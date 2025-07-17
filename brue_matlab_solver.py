import numpy as np
import matplotlib.pyplot as plt
from typing import Tuple, List, Dict
from rich.console import Console
from tqdm import tqdm

class BRUEMatlabSolver:
    def __init__(self):
        self.console = Console()
        
    def pairs_to_matrix(self, pairs: np.ndarray) -> np.ndarray:
        """将二元组转换为矩阵
        Args:
            pairs: n x 2 的矩阵，每一行是一个二元组
        Returns:
            转换后的矩阵
        """
        num_paths = 5  # 根据 range_min/range_max 的维度
        max_links = int(np.max(pairs))  # 获取最大链接编号
        matrix = np.zeros((num_paths, max_links))
        
        for i, pair in enumerate(pairs):
            if i < num_paths and int(pair[1]) <= max_links:  # 确保索引在有效范围内
                matrix[i, int(pair[1]-1)] = 1
        return matrix
    
    def calculate_real_time(self, x: np.ndarray, M: np.ndarray, 
                          free_flow_time: np.ndarray, max_capacity: np.ndarray) -> np.ndarray:
        """计算实际旅行时间
        Args:
            x: 流量
            M: 路径-链接矩阵
            free_flow_time: 自由流时间
            max_capacity: 最大容量
        Returns:
            实际旅行时间
        """
        # 找到非零列的索引
        index = np.where(np.sum(M, axis=0) != 0)[0]
        
        # 计算BPR函数
        time = self.ta(x[:, index], free_flow_time[index], max_capacity[index])
        path_time = time @ M[:, index].T
        
        # 添加感知误差
        return path_time + 15 * (1 - np.exp(-0.02 * path_time))
    
    def ta(self, x: np.ndarray, t0: np.ndarray, C: np.ndarray) -> np.ndarray:
        """交通时间计算公式
        Args:
            x: 流量
            t0: 自由流时间
            C: 容量
        Returns:
            交通时间
        """
        return t0 * (1 + 0.15 * (x / C) ** 4)
    
    def PD(self, x: np.ndarray) -> np.ndarray:
        """惩罚函数
        Args:
            x: 输入值
        Returns:
            惩罚值
        """
        return np.maximum(0, x)
    
    def evaluate_solution(self, f: np.ndarray, M: np.ndarray, zeta: float) -> Tuple[np.ndarray, np.ndarray]:
        """评估解的质量
        Args:
            f: 流量分配方案
            M: 路径-链接矩阵
            zeta: 容忍度参数
        Returns:
            目标函数值和误差值
        """
        # 定义常数
        money = np.array([20, 15, 1, 0, 0, 0, 0, 1])
        free_flow_time = np.array([18, 22.5, 12, 24, 2.4, 6, 24, 12])
        max_capacity = np.array([3600, 3600, 1800, 1800, 1800, 1800, 1800, 1800])
        
        # 计算金钱成本
        money_cost = money @ M.T
        
        # 计算流量分配到链接上的值
        x = f @ M
        
        # 计算实际旅行时间
        RT = self.calculate_real_time(x, M, free_flow_time, max_capacity)
        
        # 计算误差
        err = np.zeros(len(f))
        
        # 获取两两路径的组合
        n = M.shape[0]
        combs = [(i, j) for i in range(n) for j in range(i+1, n)]
        
        # 计算时间差异的惩罚
        for i, j in combs:
            err += self.PD(np.abs(RT[:, i] - RT[:, j]) - zeta)
            
        # 计算基于金钱差异的惩罚
        for i, j in combs:
            err += self.PD((RT[:, i] - RT[:, j]) * (money_cost[i] - money_cost[j]))
        
        # 计算目标函数
        ff = np.column_stack([
            np.sum(f * RT, axis=1),  # 总旅行时间
            f @ money_cost           # 总金钱成本
        ])
        
        return ff, err
    
    def find_brue_set(self, zeta: float, ordered_pair: np.ndarray, 
                     range_min: np.ndarray, range_max: np.ndarray,
                     num_points: int = 50) -> List[Dict]:
        """寻找BRUE集合
        Args:
            zeta: 容忍度参数
            ordered_pair: 路径-链接关系
            range_min: 每个维度的最小值
            range_max: 每个维度的最大值
            num_points: 每个维度的采样点数
        Returns:
            满足条件的解集合
        """
        # 转换路径-链接关系为矩阵
        relation_matrix = self.pairs_to_matrix(ordered_pair)
        n = len(range_min)
        
        # 生成网格点
        dim_nums = np.ones(n, dtype=int) * num_points
        samples = [np.linspace(range_min[i], range_max[i], dim_nums[i]) 
                  for i in range(n)]
        
        # 生成多维网格
        mesh = np.meshgrid(*samples)
        candidates = np.column_stack([m.flatten() for m in mesh])
        
        # 筛选总需求约束
        demand_satisfied = np.abs(np.sum(candidates, axis=1) - 10000) < 1e-6
        candidates = candidates[demand_satisfied]
        
        # 评估候选解
        ff, err = self.evaluate_solution(candidates, relation_matrix, zeta)
        
        # 筛选有效解
        valid = err == 0
        valid_solutions = []
        
        if np.any(valid):
            valid_ff = ff[valid]
            valid_flows = candidates[valid]
            
            for i in range(len(valid_ff)):
                solution = {
                    'flows': {j+1: flow for j, flow in enumerate(valid_flows[i])},
                    'objective': {
                        'time': valid_ff[i, 0],
                        'money': valid_ff[i, 1]
                    }
                }
                valid_solutions.append(solution)
                
            # 更新搜索范围
            self.console.print(f"找到的解的范围:")
            self.console.print(f"最小值: {np.min(valid_flows, axis=0)}")
            self.console.print(f"最大值: {np.max(valid_flows, axis=0)}")
        
        return valid_solutions
    
    def visualize_solutions(self, solutions: List[Dict], zeta: float):
        """可视化结果
        Args:
            solutions: 解集合
            zeta: 容忍度参数
        """
        if not solutions:
            self.console.print("[red]没有找到有效解")
            return
        
        # 提取数据
        flows = np.array([list(sol['flows'].values()) for sol in solutions])
        objectives = np.array([[sol['objective']['time'], sol['objective']['money']] 
                             for sol in solutions])
        
        # 1. 绘制流量空间
        plt.figure(figsize=(10, 8))
        plt.scatter(flows[:, 0], flows[:, 1], alpha=0.6, s=5)
        plt.title(f'Flow Space (ζ = {zeta})')
        plt.xlabel('Flow 1')
        plt.ylabel('Flow 2')
        plt.grid(True, alpha=0.3)
        plt.show()
        
        # 2. 绘制目标函数空间
        plt.figure(figsize=(10, 8))
        plt.scatter(objectives[:, 0], objectives[:, 1], alpha=0.6, s=5)
        plt.title(f'Objective Space (ζ = {zeta})')
        plt.xlabel('Total Travel Time')
        plt.ylabel('Total Money Cost')
        plt.grid(True, alpha=0.3)
        plt.show()

def main():
    solver = BRUEMatlabSolver()
    
    # 测试用例1
    ordered_pair = np.array([
        [1,1], [2,2], [3,3], [3,7], [4,4], [4,8], [5,3], [5,5], [5,8]
    ])
    
    # 设置搜索范围
    range_min = np.array([4200, 3000, 0, 0, 0])
    range_max = np.array([4900, 4200, 1610, 1610, 2000])
    
    # 使用不同的zeta值测试
    zeta_values = [15, 31]
    for zeta in zeta_values:
        print(f"\n测试 ζ = {zeta}")
        solutions = solver.find_brue_set(zeta, ordered_pair, range_min, range_max, num_points=20)
        solver.visualize_solutions(solutions, zeta)

if __name__ == "__main__":
    main() 