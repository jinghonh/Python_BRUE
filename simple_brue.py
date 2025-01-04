from brue_base import BRUEBase
from pyomo.environ import *
from rich.console import Console
from rich.table import Table


class SimpleBRUE(BRUEBase):
    def initialize_sets(self):
        self.model.od_pairs = RangeSet(1, 6)
        self.model.paths = RangeSet(1, 8)
        self.model.od_demand = Set(initialize=range(1, 7))

    def initialize_parameters(self):
        # 初始化自由流时间
        free_flow_time = {1: 18, 2: 22.5, 3: 12, 4: 24, 5: 2.4, 6: 6, 7: 24, 8: 12}
        self.model.free_flow_time = Param(self.model.paths, initialize=free_flow_time)

        # 初始化路径-链接矩阵
        path_link_matrix = {
            (1, 1): 1, (1, 2): 0, (1, 3): 0, (1, 4): 0, (1, 5): 0, (1, 6): 0, (1, 7): 0, (1, 8): 0,
            (2, 1): 0, (2, 2): 1, (2, 3): 0, (2, 4): 0, (2, 5): 0, (2, 6): 0, (2, 7): 0, (2, 8): 0,
            (3, 1): 0, (3, 2): 0, (3, 3): 1, (3, 4): 0, (3, 5): 0, (3, 6): 0, (3, 7): 1, (3, 8): 0,
            (4, 1): 0, (4, 2): 0, (4, 3): 0, (4, 4): 1, (4, 5): 0, (4, 6): 0, (4, 7): 0, (4, 8): 1,
            (5, 1): 0, (5, 2): 0, (5, 3): 1, (5, 4): 0, (5, 5): 1, (5, 6): 0, (5, 7): 0, (5, 8): 1,
            (6, 1): 0, (6, 2): 0, (6, 3): 0, (6, 4): 1, (6, 5): 0, (6, 6): 1, (6, 7): 1, (6, 8): 0
        }
        self.model.path_link_matrix = Param(self.model.od_pairs, self.model.paths,
                                            initialize=path_link_matrix, default=0)

        # 初始化链接容量
        link_capacity = {1: 3600, 2: 3600, 3: 1800, 4: 1800,
                         5: 1800, 6: 1800, 7: 1800, 8: 1800}
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
        m.demand_constraint.add(sum(m.flow[i] for i in m.od_demand) == 10000)

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

        # 可行路径约束 - 原来的版本
        m.feasible_path_constraint = ConstraintList()
        m.feasible_path_constraint.add(
            10000 - sum(m.flow[i] for i in [1, 2, 3, 4, 5]) >= 0.01  # 原来包含了1-5的路径
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
        return is_effective_path


if __name__ == "__main__":
    solver = SimpleBRUE()
    solver.run()
