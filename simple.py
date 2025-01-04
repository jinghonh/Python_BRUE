from pyomo.environ import *
from pyomo.opt import SolverFactory
from rich.console import Console
from rich.table import Table

# 创建模型
model = ConcreteModel()
model.i = RangeSet(1, 6)  # OD set
model.j = RangeSet(1, 8)  # Path set
model.OD1 = Set(initialize=range(1, 7))

# Parameters
t0 = {1: 18, 2: 22.5, 3: 12, 4: 24, 5: 2.4, 6: 6, 7: 24, 8: 12}
A_data = {
    (1, 1): 1, (1, 2): 0, (1, 3): 0, (1, 4): 0, (1, 5): 0, (1, 6): 0,(1, 7): 0,(1, 8): 0,
    (2, 1): 0, (2, 2): 1, (2, 3): 0, (2, 4): 0, (2, 5): 0, (2, 6): 0,(2, 7): 0,(2, 8): 0,
    (3, 1): 0, (3, 2): 0, (3, 3): 1, (3, 4): 0, (3, 5): 0, (3, 6): 0,(3, 7): 1,(3, 8): 0,
    (4, 1): 0, (4, 2): 0, (4, 3): 0, (4, 4): 1, (4, 5): 0, (4, 6): 0,(4, 7): 0,(4, 8): 1,
    (5, 1): 0, (5, 2): 0, (5, 3): 1, (5, 4): 0, (5, 5): 1, (5, 6): 0,(5, 7): 0,(5, 8): 1,
    (6, 1): 0, (6, 2): 0, (6, 3): 0, (6, 4): 1, (6, 5): 0, (6, 6): 1,(6, 7): 1,(6, 8): 0
}
link_capacities = {
    1: 3600, 2: 3600, 3: 1800, 4: 1800, 5: 1800, 6: 1800, 7: 1800, 8: 1800
}
model.t0 = Param(model.j, initialize=t0)
model.A = Param(model.i, model.j, initialize=A_data, default=0)
model.link_capacities = Param(model.j, initialize=link_capacities)

# variables
model.f = Var(model.i, domain=NonNegativeReals)  # Flow
model.t = Var(model.j, domain=NonNegativeReals)  # Travel time
model.c = Var(model.i, domain=NonNegativeReals)  # Cost
model.r = Var(model.i, domain=NonNegativeReals)  # Residual

model.p = Var(domain=NonNegativeReals)
model.ep = Var(domain=NonNegativeReals)

# 定义约束
model.e1 = Constraint(expr=model.x1 + model.x2 + model.x3 + model.x4 + model.x5 + model.x6 == 10000)

# 定义路径约束
def ep_constraints_rule(model, i):
    return model.ep >= model.r[i]

def balance_constraints_rule(model, i):
    return model.f[i] + model.c[i] + model.r[i] - model.p >= 0.0

def path_cost_constraints_rule(model, i):
    return (model.c[i] + model.r[i] - model.p) >= 0

def link_cost_rule(model, j):
    return model.t[j] == model.t0[j] * (1 + 0.15 * (sum(model.A[i, j] * model.f[i] for i in model.i) / model.link_capacities[j]) ** 4)

def path_cost_rule(model, i):
    return model.c[i] == sum(model.A[i, j] * model.t[j] for j in model.j)

def objective_rule(model):
    return model.ep

model.demand_constraints = ConstraintList()
model.demand_constraints.add(sum(model.f[i] for i in model.OD1) == 10000)
model.ep_constraints = Constraint(model.i, rule=ep_constraints_rule)
model.balance_constraints = Constraint(model.i, rule=balance_constraints_rule)
model.path_cost_constraints = Constraint(model.i, rule=path_cost_constraints_rule)
model.link_cost = Constraint(model.j, rule=link_cost_rule)
model.path_cost = Constraint(model.i, rule=path_cost_rule)
model.total_constraints = ConstraintList()
model.total_constraints.add(sum(model.f[i] * (model.c[i] + model.r[i] - model.p) for i in model.i) == 0)
# 添加路径约束条件
model.path_constraints = ConstraintList()
model.path_constraints.add(10000 - model.f[1]-model.f[2]-model.f[5]-model.f[3]-model.f[4]>= 0.01)
# 定义目标函数
model.obj = Objective(rule=objective_rule, sense=minimize)
# 求解
solver = SolverFactory('ipopt')
results = solver.solve(model, tee=True)

# 显示结果
console = Console()
table = Table(title="Results")
table.add_column("Flow", justify="right", style="cyan", no_wrap=True)
table.add_column("Cost", justify="right", style="green", no_wrap=True)
table.add_column("Residual", justify="right", style="yellow", no_wrap=True)
table.add_column("p", justify="right", style="blue", no_wrap=True)
table.add_column("Objective", justify="right", style="cyan", no_wrap=True)
for i in model.i:
    table.add_row(
        str(model.f[i].value),
        str(model.c[i].value),
        str(model.r[i].value),
        str(model.p.value),
        str(model.obj())
    )
console.print(table)