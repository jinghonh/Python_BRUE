from pyomo.environ import *
from pyomo.opt import SolverFactory
from rich.console import Console
from rich.table import Table

def initialize_model():
    model = ConcreteModel()
    
    # 定义集合
    model.od_pairs = RangeSet(1, 6)  # OD对集合
    model.paths = RangeSet(1, 8)     # 路径集合
    model.od_demand = Set(initialize=range(1, 7))  # OD需求集合

    # 初始化参数
    free_flow_time = {1: 18, 2: 22.5, 3: 12, 4: 24, 5: 2.4, 6: 6, 7: 24, 8: 12}
    path_link_matrix = {
        (1, 1): 1, (1, 2): 0, (1, 3): 0, (1, 4): 0, (1, 5): 0, (1, 6): 0, (1, 7): 0, (1, 8): 0,
        (2, 1): 0, (2, 2): 1, (2, 3): 0, (2, 4): 0, (2, 5): 0, (2, 6): 0, (2, 7): 0, (2, 8): 0,
        (3, 1): 0, (3, 2): 0, (3, 3): 1, (3, 4): 0, (3, 5): 0, (3, 6): 0, (3, 7): 1, (3, 8): 0,
        (4, 1): 0, (4, 2): 0, (4, 3): 0, (4, 4): 1, (4, 5): 0, (4, 6): 0, (4, 7): 0, (4, 8): 1,
        (5, 1): 0, (5, 2): 0, (5, 3): 1, (5, 4): 0, (5, 5): 1, (5, 6): 0, (5, 7): 0, (5, 8): 1,
        (6, 1): 0, (6, 2): 0, (6, 3): 0, (6, 4): 1, (6, 5): 0, (6, 6): 1, (6, 7): 1, (6, 8): 0
    }
    link_capacity = {
        1: 3600, 2: 3600, 3: 1800, 4: 1800, 5: 1800, 6: 1800, 7: 1800, 8: 1800
    }

    # 设置模型参数
    model.free_flow_time = Param(model.paths, initialize=free_flow_time)
    model.path_link_matrix = Param(model.od_pairs, model.paths, initialize=path_link_matrix, default=0)
    model.link_capacity = Param(model.paths, initialize=link_capacity)

    # 定义变量
    model.flow = Var(model.od_pairs, domain=NonNegativeReals)       # 流量
    model.travel_time = Var(model.paths, domain=NonNegativeReals)   # 旅行时间
    model.path_cost = Var(model.od_pairs, domain=NonNegativeReals)  # 路径成本
    model.residual = Var(model.od_pairs, domain=NonNegativeReals)   # 残差
    model.perception = Var(domain=NonNegativeReals)                 # 感知值
    model.epsilon = Var(domain=NonNegativeReals)                    # 误差项

    return model

def add_constraints(model):
    # 总需求约束
    model.demand_constraint = ConstraintList()
    model.demand_constraint.add(sum(model.flow[i] for i in model.od_demand) == 10000)

    # 误差约束
    model.epsilon_constraints = Constraint(
        model.od_pairs,
        rule=lambda m, i: m.epsilon >= m.residual[i]
    )

    # 平衡约束
    model.balance_constraints = Constraint(
        model.od_pairs,
        rule=lambda m, i: m.flow[i] + m.path_cost[i] + m.residual[i] - m.perception >= 0.0
    )

    # 路径成本非负约束
    model.path_cost_nonnegativity = Constraint(
        model.od_pairs,
        rule=lambda m, i: m.path_cost[i] + m.residual[i] - m.perception >= 0
    )

    # 路径成本约束
    model.path_cost_constraints = Constraint(
        model.paths,
        rule=lambda m, j: m.travel_time[j] == m.free_flow_time[j] * 
             (1 + 0.15 * (sum(m.path_link_matrix[i, j] * m.flow[i] 
              for i in m.od_pairs) / m.link_capacity[j]) ** 4)
    )

    # 路径成本计算约束
    model.path_cost_calculation = Constraint(
        model.od_pairs,
        rule=lambda m, i: m.path_cost[i] == sum(m.path_link_matrix[i, j] * m.travel_time[j] 
                                               for j in m.paths)
    )

    # 总成本约束
    model.total_cost_constraint = ConstraintList()
    model.total_cost_constraint.add(
        sum(model.flow[i] * (model.path_cost[i] + model.residual[i] - model.perception) 
            for i in model.od_pairs) == 0
    )

    # 可行路径约束 - 原来的版本
    model.feasible_path_constraint = ConstraintList()
    model.feasible_path_constraint.add(
        10000 - sum(model.flow[i] for i in [1,2,5,3,4]) >= 0.01  # 原来包含了1-5的路径
    )

def solve_model(model):
    # 设置目标函数
    model.objective = Objective(expr=model.epsilon, sense=minimize)
    
    # 求解
    solver = SolverFactory('ipopt')
    results = solver.solve(model, tee=True)
    return results

def display_results(model):
    console = Console()
    table = Table(title="交通分配结果")
    
    # 设置表格列
    columns = ["流量", "成本", "残差", "感知值", "目标函数值"]
    for col in columns:
        table.add_column(col, justify="right", style="cyan")
    
    # 添加结果行
    for i in model.od_pairs:
        table.add_row(
            f"{model.flow[i].value:.3f}",
            f"{model.path_cost[i].value:.3f}",
            f"{model.residual[i].value:.3f}",
            f"{model.perception.value:.3f}",
            f"{model.objective():.3f}"
        )
    
    console.print(table)

def main():
    model = initialize_model()
    add_constraints(model)
    results = solve_model(model)
    display_results(model)

if __name__ == "__main__":
    main()