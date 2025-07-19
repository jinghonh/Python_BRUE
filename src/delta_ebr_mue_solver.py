import numpy as np
import matplotlib.pyplot as plt
from pyomo.environ import *
from pyomo.opt import SolverFactory
from rich.console import Console
from rich.table import Table
import copy
import time

from .traffic_network_config import TrafficNetworkConfig
from .ue_solver import UESolver


class DeltaEBRMUESolver:
    """
    δ-EBR-MUE（Delta-Epsilon Bounded Rational Multi-Objective User Equilibrium）求解器
    
    实现了参考文献中的迭代算法，求解有界理性多目标用户均衡问题。
    考虑了：
    1. δ-有界理性：用户可能选择时间成本不超过最短路径时间 + δ 的路径
    2. ε-非严格支配：用户选择的路径不被其他路径ε-严格支配
    """
    
    def __init__(self, config: TrafficNetworkConfig, delta=2.0, epsilon=(0.5, 0.5)):
        """
        初始化求解器
        
        Args:
            config: 交通网络配置
            delta: 时间成本容忍度（对最短路径的容忍值）
            epsilon: 支配容忍度向量，形式为(ε_time, ε_money)
        """
        self.config = config
        self.delta = delta
        self.epsilon = epsilon  # (时间容忍度, 金钱容忍度)
        self.console = Console()
        self.results = []
        self.path_costs = []  # 存储每次迭代的路径成本
        self.effective_paths = {}  # 存储有效路径
        
    def initialize_flow(self, init_strategy='uniform'):
        """
        初始化路径流量
        
        Args:
            init_strategy: 初始化策略，可选 'uniform'（均匀分布）或 'shortest_path'（最短路径）
            
        Returns:
            初始流量向量
        """
        flows = {}
        
        if init_strategy == 'uniform':
            # 对每个OD组均匀分配流量
            for group_name, group_pairs in self.config.od_groups.items():
                demand = self.config.od_demands[group_name]
                paths_count = len(group_pairs)
                flow_per_path = demand / paths_count
                
                for path_idx in group_pairs:
                    flows[path_idx] = flow_per_path
        
        elif init_strategy == 'shortest_path':
            # 创建临时UE求解器来获取初始解
            ue_solver = UESolver(self.config)
            ue_solver.initialize_sets()
            ue_solver.initialize_parameters()
            ue_solver.initialize_variables()
            ue_solver.add_constraints()
            ue_solver.set_objective()
            ue_solver.solve(tee=False)
            
            # 从UE解中获取流量
            for i in ue_solver.model.od_pairs:
                flows[i] = ue_solver.model.flow[i].value
        
        return flows
        
    def calculate_path_costs(self, flows):
        """
        计算给定流量下的路径成本
        
        Args:
            flows: 路径流量字典
            
        Returns:
            时间成本字典和金钱成本字典
        """
        time_costs = {}
        money_costs = {}
        
        # 首先计算链接流量
        link_flows = {j: 0.0 for j in range(1, self.config.num_links + 1)}
        for i, flow in flows.items():
            for j in range(1, self.config.num_links + 1):
                if (i, j) in self.config.path_link_matrix and self.config.path_link_matrix[(i, j)] > 0:
                    link_flows[j] += flow
        
        # 计算链接时间成本（使用BPR函数）
        link_times = {}
        for j in range(1, self.config.num_links + 1):
            link_times[j] = self.config.free_flow_time[j] * (
                1 + 0.15 * (link_flows[j] / self.config.link_capacity[j]) ** 4
            )
        
        # 计算路径时间成本和金钱成本
        for i in range(1, self.config.num_paths + 1):
            # 时间成本
            time_costs[i] = sum(
                self.config.path_link_matrix.get((i, j), 0) * link_times[j]
                for j in range(1, self.config.num_links + 1)
            )
            
            # 金钱成本（固定的）
            money_costs[i] = sum(
                self.config.path_link_matrix.get((i, j), 0) * self.config.link_money_cost[j]
                for j in range(1, self.config.num_links + 1)
            )
        
        return time_costs, money_costs
    
    def find_minimum_time_costs(self, time_costs):
        """
        找出每个OD对的最小时间成本
        
        Args:
            time_costs: 路径时间成本字典
            
        Returns:
            每个OD组的最小时间成本字典
        """
        min_time_costs = {}
        
        for group_name, group_pairs in self.config.od_groups.items():
            min_time = float('inf')
            for path_idx in group_pairs:
                if time_costs[path_idx] < min_time:
                    min_time = time_costs[path_idx]
            min_time_costs[group_name] = min_time
        
        return min_time_costs
    
    def is_strictly_dominated(self, path_idx, compared_path_idx, time_costs, money_costs, use_epsilon=False):
        """
        判断一个路径是否被另一个路径严格支配
        
        如果compared_path_idx在时间和金钱两个维度上都严格优于path_idx，则path_idx被支配。
        当use_epsilon=True时，使用ε-严格支配概念。
        
        Args:
            path_idx: 被检查的路径索引
            compared_path_idx: 比较的路径索引
            time_costs: 时间成本字典
            money_costs: 金钱成本字典
            use_epsilon: 是否使用ε-严格支配概念
            
        Returns:
            如果被支配，则返回True；否则返回False
        """
        if use_epsilon:
            eps_time, eps_money = self.epsilon
            return (time_costs[compared_path_idx] + eps_time < time_costs[path_idx] and
                    money_costs[compared_path_idx] + eps_money < money_costs[path_idx])
        else:
            return (time_costs[compared_path_idx] < time_costs[path_idx] and
                    money_costs[compared_path_idx] < money_costs[path_idx])
    
    def filter_feasible_paths(self, time_costs, money_costs, min_time_costs, use_epsilon=False):
        """
        过滤满足δ-有界且非严格支配的路径
        
        Args:
            time_costs: 时间成本字典
            money_costs: 金钱成本字典
            min_time_costs: 每个OD组的最小时间成本字典
            use_epsilon: 是否使用ε-严格支配概念
            
        Returns:
            每个OD组的可行路径集合字典
        """
        feasible_paths = {}
        
        for group_name, group_pairs in self.config.od_groups.items():
            feasible_paths[group_name] = []
            min_time = min_time_costs[group_name]
            
            for path_idx in group_pairs:
                # 检查时间约束
                if time_costs[path_idx] <= min_time + self.delta:
                    # 检查是否被其他路径严格支配
                    is_dominated = False
                    for other_path_idx in group_pairs:
                        if other_path_idx != path_idx and self.is_strictly_dominated(
                            path_idx, other_path_idx, time_costs, money_costs, use_epsilon
                        ):
                            is_dominated = True
                            break
                    
                    if not is_dominated:
                        feasible_paths[group_name].append(path_idx)
        
        return feasible_paths
    
    def solve_flow_assignment(self, feasible_paths):
        """
        求解流量分配子问题
        
        在给定可行路径集合的情况下，使用凸优化求解流量分配问题。
        
        Args:
            feasible_paths: 每个OD组的可行路径集合字典
            
        Returns:
            优化后的路径流量字典
        """
        # 创建子问题的模型
        model = ConcreteModel()
        
        # 创建变量集合
        all_paths = set()
        for paths in feasible_paths.values():
            all_paths.update(paths)
        model.paths = Set(initialize=all_paths)
        model.links = RangeSet(1, self.config.num_links)
        
        # 创建变量
        model.flow = Var(model.paths, domain=NonNegativeReals)
        model.link_flow = Var(model.links, domain=NonNegativeReals)
        
        # 需求约束
        model.demand_constraints = ConstraintList()
        for group_name, paths in feasible_paths.items():
            if paths:  # 确保有可行路径
                model.demand_constraints.add(
                    sum(model.flow[p] for p in paths) == self.config.od_demands[group_name]
                )
        
        # 链接流量约束
        def link_flow_rule(model, j):
            return model.link_flow[j] == sum(
                model.flow[i] * self.config.path_link_matrix.get((i, j), 0)
                for i in model.paths if (i, j) in self.config.path_link_matrix
            )
        
        model.link_flow_constraints = Constraint(model.links, rule=link_flow_rule)
        
        # 目标函数（BPR函数积分）
        def obj_rule(model):
            result = 0
            for j in model.links:
                # BPR函数积分
                result += self.config.free_flow_time[j] * (
                    model.link_flow[j] + 0.15 * (model.link_flow[j]**5) / (5 * self.config.link_capacity[j]**4)
                )
            return result
        
        model.objective = Objective(rule=obj_rule, sense=minimize)
        
        # 求解
        solver = SolverFactory('ipopt')
        solver.options['max_iter'] = 5000
        solver.solve(model, tee=False)
        
        # 提取结果
        flows = {i: model.flow[i].value for i in model.paths}
        
        # 确保所有路径都有流量值（包括不在可行路径集中的）
        for i in range(1, self.config.num_paths + 1):
            if i not in flows:
                flows[i] = 0.0
        
        return flows
    
    def solve_flow_assignment_vi(self, feasible_paths, time_costs, money_costs):
        """
        使用变分不等式(VI)方法求解流量分配子问题
        
        基于论文中的算法实现，使用投影梯度法求解VI问题：
        找到f^(k+1) ∈ F^(δ,ε)，使得 ⟨F(f^(k+1)), f - f^(k+1)⟩ ≥ 0，∀f ∈ F^(δ,ε)
        
        Args:
            feasible_paths: 每个OD组的可行路径集合字典
            time_costs: 当前流量下的路径时间成本
            money_costs: 当前流量下的路径金钱成本
            
        Returns:
            优化后的路径流量字典
        """
        # 初始化路径流量（均匀分配给可行路径）
        flows = {}
        for group_name, paths in feasible_paths.items():
            if not paths:  # 如果没有可行路径，跳过
                continue
            demand = self.config.od_demands[group_name]
            path_count = len(paths)
            if path_count > 0:
                flow_per_path = demand / path_count
                for path_idx in paths:
                    flows[path_idx] = flow_per_path
        
        # 确保所有路径都有流量值（未在可行集中的路径流量为0）
        for i in range(1, self.config.num_paths + 1):
            if i not in flows:
                flows[i] = 0.0
        
        # 投影梯度法参数
        max_vi_iter = 100  # 最大迭代次数
        step_size = 0.01   # 步长
        vi_tolerance = 1e-5  # 收敛阈值
        
        # 投影梯度法迭代
        for vi_iter in range(max_vi_iter):
            # 计算当前流量下的路径成本
            current_time_costs, current_money_costs = self.calculate_path_costs(flows)
            
            # 构造多目标成本运算符 F(f)
            # 在这里我们使用时间成本作为主要成本，并考虑金钱成本的影响
            operator_f = {}
            for i in range(1, self.config.num_paths + 1):
                # F(f) = (t_p(f), c_p(f))，简化为加权和
                operator_f[i] = current_time_costs.get(i, 0) + 0.5 * current_money_costs.get(i, 0)
            
            # 更新流量（梯度步骤）
            new_flows = {}
            for i in flows.keys():
                if flows[i] > 0:  # 只更新正流量
                    new_flows[i] = max(0, flows[i] - step_size * operator_f[i])
                else:
                    new_flows[i] = 0
            
            # 投影回可行集（确保需求约束）
            for group_name, paths in feasible_paths.items():
                feasible_indices = [p for p in paths if p in new_flows]
                if not feasible_indices:
                    continue
                
                # 计算当前组的总流量
                total_flow = sum(new_flows[p] for p in feasible_indices)
                
                if total_flow > 0:
                    # 缩放以满足需求约束
                    scale_factor = self.config.od_demands[group_name] / total_flow
                    for p in feasible_indices:
                        new_flows[p] *= scale_factor
            
            # 检查收敛性
            flow_diff = sum(abs(new_flows.get(i, 0) - flows.get(i, 0)) for i in flows.keys())
            if flow_diff < vi_tolerance:
                break
            
            # 更新流量
            flows = new_flows.copy()
        
        # 确保所有路径都有流量值
        for i in range(1, self.config.num_paths + 1):
            if i not in flows:
                flows[i] = 0.0
        
        return flows
    
    def run(self, max_iter=50, convergence_threshold=0.01, init_strategy='uniform', use_epsilon=False, solver_method='convex'):
        """
        运行完整的δ-EBR-MUE求解算法
        
        Args:
            max_iter: 最大迭代次数
            convergence_threshold: 收敛阈值
            init_strategy: 初始流量策略
            use_epsilon: 是否使用ε-严格支配概念
            solver_method: 子问题求解方法，'convex'使用凸优化，'vi'使用变分不等式
            
        Returns:
            最终的流量分配结果
        """
        # 初始化
        flows = self.initialize_flow(init_strategy)
        
        self.results = []
        self.path_costs = []
        self.effective_paths = {}
        
        # 开始计时
        start_time = time.time()
        
        # 迭代
        for iter_count in range(max_iter):
            self.console.print(f"\n[bold]迭代 {iter_count + 1}/{max_iter}[/bold]")
            
            # 计算当前流量下的路径成本
            time_costs, money_costs = self.calculate_path_costs(flows)
            
            # 存储当前迭代的路径成本
            self.path_costs.append((copy.deepcopy(time_costs), copy.deepcopy(money_costs)))
            
            # 找出每个OD对的最小时间成本
            min_time_costs = self.find_minimum_time_costs(time_costs)
            
            # 过滤满足δ-有界且非严格支配的路径
            feasible_paths = self.filter_feasible_paths(
                time_costs, money_costs, min_time_costs, use_epsilon
            )
            
            # 存储最终迭代的有效路径
            if iter_count == max_iter - 1:
                self.effective_paths = feasible_paths
            
            # 根据选择的方法求解流量分配子问题
            if solver_method == 'vi':
                self.console.print(f"[blue]使用变分不等式(VI)方法求解子问题[/blue]")
                new_flows = self.solve_flow_assignment_vi(feasible_paths, time_costs, money_costs)
            else:
                self.console.print(f"[blue]使用凸优化方法求解子问题[/blue]")
                new_flows = self.solve_flow_assignment(feasible_paths)
            
            # 检查收敛性
            flow_diff = sum(abs(new_flows.get(i, 0) - flows.get(i, 0)) for i in range(1, self.config.num_paths + 1))
            relative_gap = flow_diff / sum(flows.values()) if sum(flows.values()) > 0 else float('inf')
            
            self.console.print(f"相对流量差: {relative_gap:.6f}")
            
            # 存储每次迭代的结果
            iteration_result = {
                'iteration': iter_count + 1,
                'flows': copy.deepcopy(new_flows),
                'time_costs': copy.deepcopy(time_costs),
                'money_costs': copy.deepcopy(money_costs),
                'feasible_paths': copy.deepcopy(feasible_paths),
                'relative_gap': relative_gap,
                'elapsed_time': time.time() - start_time,
                'solver_method': solver_method
            }
            self.results.append(iteration_result)
            
            # 判断是否收敛
            if relative_gap <= convergence_threshold:
                self.console.print(f"[green]算法已收敛，在第 {iter_count + 1} 次迭代[/green]")
                break
            
            # 更新流量
            flows = new_flows
            
            # 如果是最后一次迭代且尚未收敛，打印警告
            if iter_count == max_iter - 1:
                self.console.print("[yellow]警告：达到最大迭代次数，但算法未收敛[/yellow]")
        
        # 计算最终结果的总用时
        total_time = time.time() - start_time
        self.console.print(f"总用时: {total_time:.2f} 秒")
        
        # 返回最终的流量分配结果，同时包含迭代历史
        final_result = self.results[-1]
        final_result['results'] = self.results  # 添加完整的迭代历史
        return final_result
    
    def display_results(self, result=None):
        """显示计算结果"""
        if result is None:
            if not self.results:
                self.console.print("[red]没有可用的结果[/red]")
                return
            result = self.results[-1]
        
        flows = result['flows']
        time_costs = result['time_costs']
        money_costs = result['money_costs']
        feasible_paths = result['feasible_paths']
        
        table = Table(title="δ-EBR-MUE 计算结果")
        
        columns = ["路径", "流量", "时间成本", "金钱成本", "OD组", "是否可行"]
        for col in columns:
            table.add_column(col, justify="right")
        
        for group_name, group_pairs in self.config.od_groups.items():
            for i in group_pairs:
                is_feasible = i in feasible_paths.get(group_name, [])
                table.add_row(
                    f"{i}",
                    f"{flows.get(i, 0):.2f}",
                    f"{time_costs.get(i, 0):.2f}",
                    f"{money_costs.get(i, 0):.2f}",
                    f"{group_name}",
                    "✓" if is_feasible else "✗"
                )
        
        self.console.print(table)
        
        # 显示每个OD组的总流量和平均成本
        table_summary = Table(title="OD组汇总")
        table_summary.add_column("OD组", justify="left")
        table_summary.add_column("总流量", justify="right")
        table_summary.add_column("平均时间成本", justify="right")
        table_summary.add_column("平均金钱成本", justify="right")
        table_summary.add_column("可行路径数", justify="right")
        table_summary.add_column("总路径数", justify="right")
        
        for group_name, group_pairs in self.config.od_groups.items():
            paths_count = len(group_pairs)
            feasible_paths_count = len(feasible_paths.get(group_name, []))
            
            total_flow = sum(flows.get(i, 0) for i in group_pairs)
            avg_time = sum(flows.get(i, 0) * time_costs.get(i, 0) 
                         for i in group_pairs) / total_flow if total_flow > 0 else 0
            avg_money = sum(flows.get(i, 0) * money_costs.get(i, 0) 
                          for i in group_pairs) / total_flow if total_flow > 0 else 0
            
            table_summary.add_row(
                group_name,
                f"{total_flow:.2f}",
                f"{avg_time:.2f}",
                f"{avg_money:.2f}",
                f"{feasible_paths_count}",
                f"{paths_count}"
            )
        
        self.console.print(table_summary)
    
    def display_convergence(self):
        """显示收敛过程"""
        if not self.results:
            self.console.print("[red]没有可用的结果[/red]")
            return
        
        iterations = [r['iteration'] for r in self.results]
        gaps = [r['relative_gap'] for r in self.results]
        times = [r['elapsed_time'] for r in self.results]
        
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
        
        # 相对间隙随迭代次数变化
        ax1.plot(iterations, gaps, 'o-', linewidth=2)
        ax1.set_xlabel('迭代次数')
        ax1.set_ylabel('相对流量差')
        ax1.set_title('收敛过程 - 相对流量差')
        ax1.grid(True, linestyle='--', alpha=0.6)
        ax1.set_yscale('log')
        
        # 累计用时随迭代次数变化
        ax2.plot(iterations, times, 's-', linewidth=2, color='orange')
        ax2.set_xlabel('迭代次数')
        ax2.set_ylabel('累计用时 (秒)')
        ax2.set_title('收敛过程 - 计算时间')
        ax2.grid(True, linestyle='--', alpha=0.6)
        
        plt.tight_layout()
        plt.savefig('output/delta_ebr_mue_convergence.png')
        plt.close()
        
        self.console.print(f"收敛过程图已保存至: output/delta_ebr_mue_convergence.png")
    
    def plot_cost_analysis(self, result=None, fig_name='delta_ebr_mue_cost_analysis.png'):
        """绘制路径成本分析图"""
        if result is None:
            if not self.results:
                self.console.print("[red]没有可用的结果[/red]")
                return
            result = self.results[-1]
        
        flows = result['flows']
        time_costs = result['time_costs']
        money_costs = result['money_costs']
        feasible_paths = result['feasible_paths']
        
        fig, ax = plt.subplots(figsize=(12, 8))
        
        # 添加验证，确保有可绘制的路径
        has_valid_paths = False
        
        for group_name, group_pairs in self.config.od_groups.items():
            # 确保该OD组有有效路径
            valid_group_pairs = [i for i in group_pairs if i in time_costs and i in money_costs]
            if not valid_group_pairs:
                continue
                
            has_valid_paths = True
            min_time = min(time_costs.get(i, float('inf')) for i in valid_group_pairs)
            min_money = min(money_costs.get(i, float('inf')) for i in valid_group_pairs)
            
            # 绘制各路径的成本点
            for i in valid_group_pairs:
                is_feasible = i in feasible_paths.get(group_name, [])
                flow = flows.get(i, 0)
                time_cost = time_costs.get(i, 0)
                money_cost = money_costs.get(i, 0)
                
                marker = 'o' if is_feasible else 'x'
                color = 'green' if flow > 0.01 else 'red'
                
                # 设置大小，但确保有最小值以保证可见性
                size = max(20, flow)  # 确保即使流量为0也能显示
                
                ax.scatter(time_cost, money_cost, s=size, c=color, marker=marker, alpha=0.7)
                ax.annotate(f"{i}", (time_cost, money_cost), xytext=(5, 5), 
                          textcoords='offset points')
            
            # 绘制最短路径到δ边界的线
            ax.axvline(x=min_time + self.delta, color='blue', linestyle='--', alpha=0.5)
            ax.axvline(x=min_time, color='green', linestyle='-', alpha=0.5)
            
            # 如果使用了ε-非支配，则绘制ε-边界
            eps_time, eps_money = self.epsilon
            used_feasible_paths = [i for i in valid_group_pairs 
                                  if i in feasible_paths.get(group_name, []) and flows.get(i, 0) > 0.01]
            
            for i in used_feasible_paths:
                time_cost = time_costs.get(i, 0)
                money_cost = money_costs.get(i, 0)
                
                # 绘制ε-支配区域边界
                ax.plot([time_cost - eps_time, time_cost - eps_time], 
                      [money_cost - eps_money, money_cost], 
                      'k--', alpha=0.3)
                ax.plot([time_cost - eps_time, time_cost], 
                      [money_cost - eps_money, money_cost - eps_money], 
                      'k--', alpha=0.3)
        
        # 如果没有有效路径，显示提示信息
        if not has_valid_paths:
            self.console.print("[red]没有有效的路径数据可绘制[/red]")
            return
        
        ax.set_xlabel('时间成本')
        ax.set_ylabel('金钱成本')
        ax.set_title(f'δ-EBR-MUE 路径成本分析 (δ={self.delta}, ε={self.epsilon})')
        ax.grid(True, linestyle='--', alpha=0.6)
        
        # 添加图例
        from matplotlib.lines import Line2D
        legend_elements = [
            Line2D([0], [0], marker='o', color='w', markerfacecolor='green', markersize=10, label='使用的可行路径'),
            Line2D([0], [0], marker='x', color='red', markersize=10, label='未使用的路径'),
            Line2D([0], [0], color='green', label='最短路径时间'),
            Line2D([0], [0], color='blue', linestyle='--', label=f'δ边界 ({self.delta})'),
            Line2D([0], [0], color='k', linestyle='--', label=f'ε边界 ({self.epsilon})')
        ]
        ax.legend(handles=legend_elements, loc='best')
        
        plt.tight_layout()
        plt.savefig(f'output/{fig_name}')
        plt.close()
        
        self.console.print(f"路径成本分析图已保存至: output/{fig_name}")
    
    def compare_delta_sensitivity(self, deltas=[0.5, 1.0, 2.0, 5.0], max_iter=20, convergence_threshold=0.01):
        """比较不同δ值对结果的影响"""
        results = {}
        
        for delta in deltas:
            self.delta = delta
            self.console.print(f"\n[bold]求解 δ = {delta}[/bold]")
            result = self.run(max_iter=max_iter, convergence_threshold=convergence_threshold)
            results[delta] = result
        
        # 分析不同δ值下的路径使用情况
        fig, ax = plt.subplots(figsize=(12, 8))
        
        bar_width = 0.15
        index = np.arange(self.config.num_paths)
        
        for i, delta in enumerate(deltas):
            flows = results[delta]['flows']
            flow_values = [flows.get(p+1, 0) for p in index]
            ax.bar(index + i * bar_width, flow_values, bar_width,
                 label=f'δ = {delta}')
        
        ax.set_xlabel('路径ID')
        ax.set_ylabel('流量')
        ax.set_title('不同δ值下的路径流量分布')
        ax.set_xticks(index + bar_width * (len(deltas) - 1) / 2)
        ax.set_xticklabels([f'{i+1}' for i in index])
        ax.legend()
        ax.grid(True, linestyle='--', alpha=0.6)
        
        plt.tight_layout()
        plt.savefig('output/delta_sensitivity_paths.png')
        plt.close()
        
        self.console.print(f"δ敏感性分析图已保存至: output/delta_sensitivity_paths.png")
        
        # 恢复初始δ值
        self.delta = deltas[-1]
        
        return results


def main():
    """主函数，用于测试δ-EBR-MUE求解器"""
    # 创建网络配置
    config = TrafficNetworkConfig.create_multi_od_network()
    
    # 初始化并运行求解器
    solver = DeltaEBRMUESolver(config, delta=2.0, epsilon=(0.5, 0.5))
    result = solver.run(max_iter=20, convergence_threshold=0.01)
    
    # 显示结果
    solver.display_results(result)
    solver.display_convergence()
    solver.plot_cost_analysis(result)
    
    # 比较不同δ值的影响
    solver.compare_delta_sensitivity([0.5, 1.0, 2.0, 5.0])
    
    return solver, result


if __name__ == "__main__":
    main() 