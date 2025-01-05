from brue_base import BRUEBase
from pyomo.environ import *
from rich.console import Console
from rich.table import Table
import matplotlib.pyplot as plt


class SimpleBRUE(BRUEBase):
    def initialize_sets(self):
        self.model.od_pairs = RangeSet(1, 8)
        self.model.paths = RangeSet(1, 19)
        self.model.od_demand = Set(initialize=range(1, 9))
    def initialize_parameters(self):
        # 初始化自由流时间
        free_flow_time = {1: 10, 2: 10, 3: 10, 4: 40, 5: 10, 6: 10, 7: 10, 8: 20, 9: 10,
      10: 25, 11: 10, 12: 10, 13: 40, 14: 10, 15: 10, 16: 10, 17: 10, 18: 80, 19: 10}
        link_money_cost = {1: 10, 2: 10, 3: 10, 4: 20, 5: 10, 6: 10, 7: 10, 8: 20, 9: 10,
      10: 25, 11: 10, 12: 10, 13: 20, 14: 10, 15: 10, 16: 10, 17: 10, 18: 30, 19: 10}
        self.model.free_flow_time = Param(self.model.paths, initialize=free_flow_time)
        self.model.link_money_cost = Param(self.model.paths, initialize=link_money_cost)

A_data = {
    (1, 1): 1, (1, 5): 1, (1, 7): 1, (1, 9): 1, (1, 11): 1,
    (2, 1): 1, (2, 5): 1, (2, 7): 1, (2, 10): 1, (2, 15): 1,
    (3, 1): 1, (3, 5): 1, (3, 8): 1, (3, 14): 1, (3, 15): 1,
    (4, 1): 1, (4, 6): 1, (4, 12): 1, (4, 14): 1, (4, 15): 1,
    (5, 2): 1, (5, 7): 1, (5, 9): 1, (5, 11): 1, (5, 17): 1,
    (6, 2): 1, (6, 7): 1, (6, 10): 1, (6, 15): 1, (6, 17): 1,
    (7, 2): 1, (7, 8): 1, (7, 14): 1, (7, 15): 1, (7, 17): 1,
    (8, 2): 1, (8, 11): 1, (8, 18): 1,
    (9, 1): 1, (9, 5): 1, (9, 7): 1, (9, 10): 1, (9, 16): 1,
    (10, 1): 1, (10, 5): 1, (10, 8): 1, (10, 14): 1, (10, 16): 1,
    (11, 1): 1, (11, 6): 1, (11, 12): 1, (11, 14): 1, (11, 16): 1,
    (12, 1): 1, (12, 13): 1, (12, 19): 1,
    (13, 2): 1, (13, 7): 1, (13, 10): 1, (13, 16): 1, (13, 17): 1,
    (14, 2): 1, (14, 8): 1, (14, 14): 1, (14, 16): 1, (14, 17): 1
}

from brue_base import BRUEBase
from pyomo.environ import *
from rich.console import Console
from rich.table import Table
import matplotlib.pyplot as plt


class SimpleBRUE(BRUEBase):
    def initialize_sets(self):
        self.model.od_pairs = RangeSet(1, 8)
        self.model.paths = RangeSet(1, 19)
        self.model.od_demand = Set(initialize=range(1, 9))
    def initialize_parameters(self):
        # 初始化自由流时间
        free_flow_time = {1: 10, 2: 10, 3: 10, 4: 40, 5: 10, 6: 10, 7: 10, 8: 20, 9: 10,
      10: 25, 11: 10, 12: 10, 13: 40, 14: 10, 15: 10, 16: 10, 17: 10, 18: 80, 19: 10}
        link_money_cost = {1: 10, 2: 10, 3: 10, 4: 20, 5: 10, 6: 10, 7: 10, 8: 20, 9: 10,
      10: 25, 11: 10, 12: 10, 13: 20, 14: 10, 15: 10, 16: 10, 17: 10, 18: 30, 19: 10}
        self.model.free_flow_time = Param(self.model.paths, initialize=free_flow_time)
        self.model.link_money_cost = Param(self.model.paths, initialize=link_money_cost)

        # 初始化路径-链接矩阵
        path_link_matrix = {
            (1, 1): 1, (1, 5): 1, (1, 7): 1, (1, 9): 1, (1, 11): 1,
            (2, 1): 1, (2, 5): 1, (2, 7): 1, (2, 10): 1, (2, 15): 1,
            (3, 1): 1, (3, 5): 1, (3, 8): 1, (3, 14): 1, (3, 15): 1,
            (4, 1): 1, (4, 6): 1, (4, 12): 1, (4, 14): 1, (4, 15): 1,
            (5, 2): 1, (5, 7): 1, (5, 9): 1, (5, 11): 1, (5, 17): 1,
            (6, 2): 1, (6, 7): 1, (6, 10): 1, (6, 15): 1, (6, 17): 1,
            (7, 2): 1, (7, 8): 1, (7, 14): 1, (7, 15): 1, (7, 17): 1,
            (8, 2): 1, (8, 11): 1, (8, 18): 1,
        }
        self.model.path_link_matrix = Param(self.model.od_pairs, self.model.paths,
                                            initialize=path_link_matrix, default=0)

        # 初始化链接容量
        link_capacity = {1: 2500, 2: 2500, 3: 2500, 4: 2500,
                         5: 2500, 6: 2500, 7: 2500, 8: 2500,
                         9: 2500, 10: 2500, 11: 2500, 12: 2500,
                         13: 2500, 14: 2500, 15: 2500, 16: 2500,
                         17: 2500, 18: 2500, 19: 2500}
        self.model.link_capacity = Param(self.model.paths, initialize=link_capacity)

    def initialize_variables(self):
        m = self.model
        m.flow = Var(m.od_pairs, domain=NonNegativeReals)
        m.travel_time = Var(m.paths, domain=NonNegativeReals)
        m.path_cost = Var(m.od_pairs, domain=NonNegativeReals)
        m.residual = Var(m.od_pairs, domain=NonNegativeReals)
        m.perception = Var(domain=NonNegativeReals)
        m.epsilon = Var(domain=NonNegativeReals)

    def add_constraints(self):
        m = self.model

        # 总需求约束
        m.demand_constraint = ConstraintList()
        m.demand_constraint.add(sum(m.flow[i] for i in m.od_demand) == 3000)

        # 误差约束
        m.epsilon_constraints = Constraint(
            m.od_pairs,
            rule=lambda m, i: m.epsilon >= m.residual[i]
        )

        # 平衡约束
        m.balance_constraints = Constraint(
            m.od_pairs,
            rule=lambda m, i: m.flow[i] + m.path_cost[i] + m.residual[i] - m.perception >= 0.0
        )

        # 路径成本非负约束
        m.path_cost_nonnegativity = Constraint(
            m.od_pairs,
            rule=lambda m, i: m.path_cost[i] + m.residual[i] - m.perception >= 0
        )

        # 路径成本约束
        m.path_cost_constraints = Constraint(
            m.paths,
            rule=lambda m, j: m.travel_time[j] == m.free_flow_time[j] *
                              (1 + 0.15 * (sum(m.path_link_matrix[i, j] * m.flow[i]
                                               for i in m.od_pairs) / m.link_capacity[j]) ** 4)
        )

        # 路径成本计算约束
        m.path_cost_calculation = Constraint(
            m.od_pairs,
            rule=lambda m, i: m.path_cost[i] == sum(m.path_link_matrix[i, j] * m.travel_time[j]
                                                    for j in m.paths)
        )

        # 总成本约束
        m.total_cost_constraint = ConstraintList()
        m.total_cost_constraint.add(
            sum(m.flow[i] * (m.path_cost[i] + m.residual[i] - m.perception)
                for i in m.od_pairs) == 0
        )

    def set_objective(self):
        self.model.objective = Objective(expr=self.model.epsilon, sense=minimize)

    def display_results(self):
        table = Table(title="交通分配结果")

        columns = ["流量", "成本", "残差", "感知值", "目标函数值"]
        for col in columns:
            table.add_column(col, justify="right", style="cyan")

        for i in self.model.od_pairs:
            table.add_row(
                f"{self.model.flow[i].value:.3f}",
                f"{self.model.path_cost[i].value:.3f}",
                f"{self.model.residual[i].value:.3f}",
                f"{self.model.perception.value:.3f}",
                f"{self.model.objective():.3f}"
            )

        self.console.print(table)

    def calculate_money_cost(self):
        """计算每条路径的金钱成本"""
        m = self.model
        money_costs = {}
        for i in m.od_pairs:
            money_costs[i] = sum(m.path_link_matrix[i, j] * m.link_money_cost[j]
                                 for j in m.paths)
        return money_costs

    def is_dominated(self, path_i, path_costs, money_costs, all_paths):
        """
        判断一个路径是否被其他路径支配
        路径i被路径j支配的条件：
        1. j的时间成本 <= i的时间成本
        2. j的金钱成本 <= i的金钱成本
        3. 至少有一个严格小于
        """
        cost_i = path_costs[path_i]
        money_i = money_costs[path_i]

        for path_j in all_paths:
            if path_j != path_i:
                cost_j = path_costs[path_j]
                money_j = money_costs[path_j]

                # 检查是否被支配
                if (cost_j <= cost_i and money_j <= money_i and
                        (cost_j < cost_i or money_j < money_i)):
                    return True
        return False

    def plot_cost_analysis(self, is_effective_path, iteration_data=None):
        """
        绘制路径成本与金钱成本的散点图
        Args:
            is_effective_path: 有效路径列表
            iteration_data: 包含当前迭代信息的字典，包括epsilon和restricted_paths
        """
        m = self.model

        # 计算金钱成本
        money_costs = self.calculate_money_cost()

        # 获取路径成本
        path_costs = {i: m.path_cost[i].value for i in m.od_pairs}

        # 获取最小成本和上界
        min_cost = min(path_costs.values())
        epsilon = m.epsilon.value
        upper_bound = min_cost + epsilon

        # 首先确定可行路径（在有效区间内的路径）
        feasible_paths = [i for i in m.od_pairs
                          if min_cost <= path_costs[i] <= upper_bound + 5e-2]
        infeasible_paths = [i for i in m.od_pairs if i not in feasible_paths]

        # 在可行路径中分析支配关系
        dominated_paths = []
        non_dominated_paths = []

        # 首先找出非支配路径
        for path in feasible_paths:
            if not self.is_dominated(path, path_costs, money_costs, feasible_paths):
                non_dominated_paths.append(path)
            else:
                dominated_paths.append(path)

        # 创建图形
        plt.figure(figsize=(10, 6))

        # 绘制不可行路径点（灰色）
        if infeasible_paths:
            plt.scatter([path_costs[i] for i in infeasible_paths],
                        [money_costs[i] for i in infeasible_paths],
                        color='grey', alpha=0.5, label='Infeasible paths')

        # 绘制被支配的可行路径点（红色）
        if dominated_paths:
            plt.scatter([path_costs[i] for i in dominated_paths],
                        [money_costs[i] for i in dominated_paths],
                        color='red', alpha=0.7, label='Dominated feasible paths')

        # 绘制非支配路径点（绿色）
        if non_dominated_paths:
            plt.scatter([path_costs[i] for i in non_dominated_paths],
                        [money_costs[i] for i in non_dominated_paths],
                        color='green', alpha=1.0, label='Non-dominated paths')

        # 添加路径标签和成本值
        for i in m.od_pairs:
            plt.annotate(f'P{i}\n({path_costs[i]:.1f}, {money_costs[i]:.1f})',
                         (path_costs[i], money_costs[i]),
                         xytext=(5, 5), textcoords='offset points',
                         fontsize=8)

        # 绘制有效区间的竖线
        plt.axvline(x=min_cost, color='blue', linestyle='--',
                    alpha=0.5, label='Min cost')
        plt.axvline(x=upper_bound, color='blue', linestyle='--',
                    alpha=0.5, label='Max cost')

        # 设置图形属性
        plt.xlabel('Travel Time Cost')
        plt.ylabel('Money Cost')

        # 添加迭代信息到标题
        if iteration_data:
            plt.title(f'Iteration {iteration_data["iteration"]}\n'
                      f'Restricted Paths: {iteration_data["restricted_paths"]}\n'
                      f'Epsilon: {epsilon:.2f}')
        else:
            plt.title('Path Cost Analysis with Dominance Relations')

        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.grid(True, alpha=0.3)

        # 添加说明文本
        info_text = (
            f'Feasible paths: {sorted(feasible_paths)}\n'
            f'Non-dominated paths: {sorted(non_dominated_paths)}\n'
            f'Dominated paths: {sorted(dominated_paths)}\n'
            f'Infeasible paths: {sorted(infeasible_paths)}\n'
            f'Current epsilon: {epsilon:.2f}'
        )
        plt.text(0.6, 1, info_text,
                 transform=plt.gca().transAxes,
                 verticalalignment='top',
                 bbox=dict(boxstyle='round', facecolor='white', alpha=0.8),
                 fontsize=8)

        # 调整布局以适应文本
        plt.tight_layout()

        # 显示图形
        plt.show()

        # 返回分析结果
        return {
            'non_dominated_paths': non_dominated_paths,
            'dominated_paths': dominated_paths,
            'infeasible_paths': infeasible_paths,
            'feasible_paths': feasible_paths,
            'current_epsilon': epsilon,
            'path_costs': path_costs,
            'money_costs': money_costs
        }

    def analyze_path_costs(self):
        """分析路径成本并找出有效路径"""
        m = self.model
        is_effective_path = []
        # 获取最小成本
        min_cost = min(m.path_cost[i].value for i in m.od_pairs)
        epsilon = m.epsilon.value
        upper_bound = min_cost + epsilon

        # 创建结果表格
        table = Table(title="路径成本分析")
        table.add_column("路径编号", justify="right", style="cyan")
        table.add_column("路径成本", justify="right", style="magenta")
        table.add_column("流量", justify="right", style="green")
        table.add_column("是否在有效区间", justify="center", style="yellow")

        self.console.print(f"\n成本分析区间: [{min_cost:.3f}, {upper_bound:.3f}]")

        # 分析每条路径
        for i in m.od_pairs:
            cost = m.path_cost[i].value
            flow = m.flow[i].value
            is_effective = min_cost <= cost <= upper_bound + 5e-2

            table.add_row(
                str(i),
                f"{cost:.3f}",
                f"{flow:.3f}",
                "✓" if is_effective else "✗"
            )
            if is_effective:
                is_effective_path.append(i)

        self.console.print(table)

        # 输出实际使用的路径（有流量的路径）
        used_paths = [i for i in m.od_pairs if m.flow[i].value > 1e-6]
        self.console.print(f"\n实际使用的路径: {used_paths}")

        # return 有效路径 is_effective
        self.console.print("\n有效路径: ", is_effective_path)

        # 在函数末尾添加绘图
        # self.plot_cost_analysis(is_effective_path)
        return is_effective_path

    def run(self):
        """运行完整求解过程"""
        super().run()
        is_effective_path = self.analyze_path_costs()
        self.plot_cost_analysis(is_effective_path)

    def add_feasible_path_constraints(self, restricted_paths, min_epsilon=None):
        """
        添加可行路径约束
        Args:
            restricted_paths: 限制的路径列表
            min_epsilon: 最小epsilon值
        """
        m = self.model
        m.feasible_path_constraint = ConstraintList()
        m.feasible_path_constraint.add(
            3000 - sum(m.flow[i] for i in restricted_paths) >= 0.01
        )
        if min_epsilon is not None:
            m.feasible_path_constraint.add(
                m.epsilon >= min_epsilon
            )
        # m.feasible_path_constraint.add(
        #     sum(m.flow[i] for i in [1,2]) == 0
        # )

    def iterate_path_analysis(self, path_groups):
        """
        迭代分析不同路径组合的epsilon和有效路径
        Args:
            path_groups: 路径组列表的列表，例如 [[1,2,5], [1,2,5,3], [1,2,5,3,4]]
        Returns:
            list of dict: 每次迭代的结果，包含epsilon和有效路径
        """
        results = []
        prev_epsilon = None

        for paths in path_groups:
            # 重新初始化模型
            self.model = ConcreteModel()
            self.initialize_sets()
            self.initialize_parameters()
            self.initialize_variables()

            # 添加基本约束
            self.add_constraints()

            # 添加可行路径约束
            self.add_feasible_path_constraints(paths, prev_epsilon)

            # 设置目标函数和求解
            self.set_objective()
            self.solve()

            # 分析结果
            is_effective_path = self.analyze_path_costs()
            current_epsilon = self.model.epsilon.value

            # 存储结果
            iteration_result = {
                'restricted_paths': paths,
                'epsilon': current_epsilon,
                'effective_paths': is_effective_path
            }
            results.append(iteration_result)

            # 更新前一次的epsilon
            prev_epsilon = current_epsilon

        return results

    def display_iteration_results(self, results):
        """
        显示迭代结果
        Args:
            results: iterate_path_analysis的返回结果
        """
        table = Table(title="路径约束迭代分析结果")
        table.add_column("迭代", style="cyan")
        table.add_column("限制路径", style="magenta")
        table.add_column("Epsilon", style="yellow")
        table.add_column("有效路径", style="green")

        for i, result in enumerate(results, 1):
            table.add_row(
                str(i),
                str(result['restricted_paths']),
                f"{result['epsilon']:.3f}",
                str(result['effective_paths'])
            )

        self.console.print(table)

    def run_path_iteration_analysis(self):
        """运行路径迭代分析"""
        results = []
        current_paths = [1]  # 从路径1开始
        prev_epsilon = None
        iteration = 1

        while True:
            # 重新初始化模型
            self.model = ConcreteModel()
            self.initialize_sets()
            self.initialize_parameters()
            self.initialize_variables()
            self.add_constraints()
            self.add_feasible_path_constraints(current_paths, prev_epsilon)
            self.set_objective()
            self.solve()

            # 分析结果
            is_effective_path = self.analyze_path_costs()
            current_epsilon = self.model.epsilon.value

            # 创建当前迭代的数据字典
            iteration_data = {
                'iteration': iteration,
                'restricted_paths': current_paths,
                'epsilon': current_epsilon,
                'effective_paths': is_effective_path
            }

            # 存储结果并创建图表
            results.append(iteration_data)
            self.plot_cost_analysis(is_effective_path, iteration_data)

            # 检查是否所有有效路径都已包含在限制路径中
            new_paths = [p for p in is_effective_path if p not in current_paths]

            # 如果没有新的有效路径，或者所有路径都已包含，则终止迭代
            if not new_paths or set(current_paths) == set(range(1, 7)):
                break

            # 更新路径集合和epsilon
            current_paths.extend(new_paths)
            current_paths.sort()
            prev_epsilon = current_epsilon
            iteration += 1

            self.console.print(f"\n迭代 {iteration}:")
            self.console.print(f"当前限制路径: {current_paths}")
            self.console.print(f"新增有效路径: {new_paths}")

        # 显示迭代结果汇总
        self.display_iteration_results(results)
        return results

    def display_iteration_results(self, results):
        """
        显示迭代结果
        Args:
            results: iterate_path_analysis的返回结果
        """
        table = Table(title="路径约束迭代分析结果")
        table.add_column("迭代", style="cyan")
        table.add_column("限制路径", style="magenta")
        table.add_column("新增有效路径", style="blue")
        table.add_column("Epsilon", style="yellow")
        table.add_column("所有有效路径", style="green")

        for i, result in enumerate(results):
            # 计算新增的有效路径
            prev_paths = set() if i == 0 else set(results[i - 1]['effective_paths'])
            new_effective = set(result['effective_paths']) - prev_paths

            table.add_row(
                str(result['iteration']),
                str(result['restricted_paths']),
                str(sorted(list(new_effective))),
                f"{result['epsilon']:.3f}",
                str(result['effective_paths'])
            )

        self.console.print(table)

    def plot_initial_costs(self):
        """
        绘制三种情况下的时间和金钱成本的散点图：
        1. 初始状态 (t0)
        2. 无流量状态
        3. 满容量状态
        """
        m = self.model
        fig, axes = plt.subplots(1, 3, figsize=(18, 6))

        # 计算三种情况下的时间成本
        scenarios = {
            'Initial State': {
                i: sum(m.path_link_matrix[i, j] * m.free_flow_time[j]
                       for j in m.paths) for i in m.od_pairs
            },
            'No Flow': {
                i: sum(m.path_link_matrix[i, j] * m.free_flow_time[j]
                       for j in m.paths) for i in m.od_pairs
            },
            'Full Capacity': {
                i: sum(m.path_link_matrix[i, j] * m.free_flow_time[j] * 
                       (1 + 0.15 * (1.0) ** 4)  # 假设流量/容量比为1
                       for j in m.paths) for i in m.od_pairs
            }
        }

        # 获取金钱成本（对所有情况都相同）
        money_costs = self.calculate_money_cost()

        # 为每种情况分析和绘图
        for idx, (scenario_name, time_costs) in enumerate(scenarios.items()):
            # 分析支配关系
            dominated_paths = []
            non_dominated_paths = []
            
            # 检查每条路径是否被支配
            for path_i in m.od_pairs:
                is_dominated = False
                time_i = time_costs[path_i]
                money_i = money_costs[path_i]
                
                for path_j in m.od_pairs:
                    if path_j != path_i:
                        time_j = time_costs[path_j]
                        money_j = money_costs[path_j]
                        
                        if (time_j <= time_i and money_j <= money_i and 
                            (time_j < time_i or money_j < money_i)):
                            is_dominated = True
                            break
                
                if is_dominated:
                    dominated_paths.append(path_i)
                else:
                    non_dominated_paths.append(path_i)

            # 绘制散点图
            ax = axes[idx]
            
            # 绘制被支配的路径点（红色）
            if dominated_paths:
                ax.scatter([time_costs[i] for i in dominated_paths],
                          [money_costs[i] for i in dominated_paths],
                          color='red', alpha=0.7, label='Dominated paths')
            
            # 绘制非支配路径点（绿色）
            if non_dominated_paths:
                ax.scatter([time_costs[i] for i in non_dominated_paths],
                          [money_costs[i] for i in non_dominated_paths],
                          color='green', alpha=1.0, label='Non-dominated paths')
            
            # 添加路径标签和成本值
            for i in m.od_pairs:
                ax.annotate(f'P{i}\n({time_costs[i]:.1f}, {money_costs[i]:.1f})',
                           (time_costs[i], money_costs[i]),
                           xytext=(5, 5), textcoords='offset points',
                           fontsize=8)
            
            # 设置图形属性
            ax.set_xlabel('Travel Time Cost')
            ax.set_ylabel('Money Cost')
            ax.set_title(f'{scenario_name}\n'
                        f'Non-dominated: {sorted(non_dominated_paths)}\n'
                        f'Dominated: {sorted(dominated_paths)}')
            ax.legend()
            ax.grid(True, alpha=0.3)

        # 调整布局
        plt.tight_layout()
        plt.show()

        # 返回分析结果
        return {
            'scenarios': {
                name: {
                    'time_costs': costs,
                    'money_costs': money_costs,
                    'non_dominated_paths': [p for p in m.od_pairs 
                                          if not self.is_dominated(p, costs, money_costs, m.od_pairs)],
                    'dominated_paths': [p for p in m.od_pairs 
                                      if self.is_dominated(p, costs, money_costs, m.od_pairs)]
                }
                for name, costs in scenarios.items()
            }
        }

if __name__ == "__main__":
    solver = SimpleBRUE()
    solver.run_path_iteration_analysis()

    solver = SimpleBRUE()
    # 初始化模型参数
    solver.initialize_sets()
    solver.initialize_parameters()
    # 绘制初始成本分析图
    initial_analysis = solver.plot_initial_costs()
    # 打印详细结果
    print("\n初始路径分析结果:")
    try:
        print(f"非支配路径: {sorted(initial_analysis['non_dominated_paths'])}")
    except KeyError:
        print("无非支配路径")
    try:
        print(f"被支配路径: {sorted(initial_analysis['dominated_paths'])}")
    except KeyError:
        print("无被支配路径") 
