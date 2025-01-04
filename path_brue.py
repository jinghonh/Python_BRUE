from brue_base import BRUEBase
from pyomo.environ import *
from rich.console import Console
from rich.table import Table

class PathBRUE(BRUEBase):
    def initialize_sets(self):
        self.model.i = RangeSet(1, 14)  # OD set
        self.model.j = RangeSet(1, 19)  # Path set
        self.model.OD1 = Set(initialize=range(1, 9))
        self.model.OD2 = Set(initialize=range(9, 15))

    def initialize_parameters(self):
        # 初始化自由流时间
        t0 = {1: 18, 2: 22.5, 3: 12, 4: 24, 5: 2.4, 6: 6, 7: 24, 8: 12}
        self.model.t0 = Param(self.model.j, initialize=t0)

        # 初始化路径-链接矩阵
        A_data = {
        (1, 1): 1, (1, 2): 0, (1, 3): 0, (1, 4): 0, (1, 5): 0, (1, 6): 0, (1, 7): 0, (1, 8): 0,
        (2, 1): 0, (2, 2): 1, (2, 3): 0, (2, 4): 0, (2, 5): 0, (2, 6): 0, (2, 7): 0, (2, 8): 0,
        (3, 1): 0, (3, 2): 0, (3, 3): 1, (3, 4): 0, (3, 5): 0, (3, 6): 0, (3, 7): 1, (3, 8): 0,
        (4, 1): 0, (4, 2): 0, (4, 3): 0, (4, 4): 1, (4, 5): 0, (4, 6): 0, (4, 7): 0, (4, 8): 1,
        (5, 1): 0, (5, 2): 0, (5, 3): 1, (5, 4): 0, (5, 5): 1, (5, 6): 0, (5, 7): 0, (5, 8): 1,
        (6, 1): 0, (6, 2): 0, (6, 3): 0, (6, 4): 1, (6, 5): 0, (6, 6): 1, (6, 7): 1, (6, 8): 0
    }
        self.model.A = Param(self.model.i, self.model.j, initialize=A_data, default=0)

    def initialize_variables(self):
        m = self.model
        m.f = Var(m.i, domain=NonNegativeReals)
        m.t = Var(m.j, domain=NonNegativeReals)
        m.c = Var(m.i, domain=NonNegativeReals)
        m.r = Var(m.i, domain=NonNegativeReals)
        m.p1 = Var(domain=NonNegativeReals)
        m.p2 = Var(domain=NonNegativeReals)
        m.ep1 = Var(domain=NonNegativeReals)
        m.ep2 = Var(domain=NonNegativeReals)

    def add_constraints(self):
        # ... 实现所有约束 ...
        pass

    def set_objective(self):
        self.model.obj = Objective(
            expr=self.model.ep1 * 4 + self.model.ep2 * 7,
            sense=minimize
        )

    def display_results(self):
        # ... 实现结果显示 ...
        pass

    def analyze_path_costs(self):
        """分析路径成本并找出有效路径"""
        m = self.model
        
        # 分别分析OD1和OD2的路径
        for od_group, epsilon in [("OD1", m.ep1.value), ("OD2", m.ep2.value)]:
            paths = m.OD1 if od_group == "OD1" else m.OD2
            
            # 获取该组的最小成本
            min_cost = min(m.c[i].value for i in paths)
            upper_bound = min_cost + epsilon
            
            # 创建结果表格
            table = Table(title=f"{od_group}路径成本分析")
            table.add_column("路径编号", justify="right", style="cyan")
            table.add_column("路径成本", justify="right", style="magenta")
            table.add_column("流量", justify="right", style="green")
            table.add_column("是否在有效区间", justify="center", style="yellow")
            
            self.console.print(f"\n{od_group}成本分析区间: [{min_cost:.3f}, {upper_bound:.3f}]")
            
            # 分析每条路径
            for i in paths:
                cost = m.c[i].value
                flow = m.f[i].value
                is_effective = min_cost <= cost <= upper_bound
                
                table.add_row(
                    str(i),
                    f"{cost:.3f}",
                    f"{flow:.3f}",
                    "✓" if is_effective else "✗"
                )
            
            self.console.print(table)
            
            # 输出实际使用的路径
            used_paths = [i for i in paths if m.f[i].value > 1e-6]
            self.console.print(f"\n{od_group}实际使用的路径: {used_paths}")

if __name__ == "__main__":
    solver = PathBRUE()
    solver.run()
