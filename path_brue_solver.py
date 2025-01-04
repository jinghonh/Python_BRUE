from pyomo.environ import *
from pyomo.opt import SolverFactory
import rich
from rich.console import Console
from rich.table import Table

# Define model
model = ConcreteModel()

# Sets and indices
model.i = RangeSet(1, 14)  # OD set
model.j = RangeSet(1, 19)  # Path set
model.OD1 = Set(initialize=range(1, 9))
model.OD2 = Set(initialize=range(9, 15))

# Parameters
t0 = {1: 10, 2: 10, 3: 10, 4: 40, 5: 10, 6: 10, 7: 10, 8: 20, 9: 10,
      10: 25, 11: 10, 12: 10, 13: 40, 14: 10, 15: 10, 16: 10, 17: 10, 18: 80, 19: 10}
model.t0 = Param(model.j, initialize=t0)

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
model.A = Param(model.i, model.j, initialize=A_data, default=0)

# Variables
model.f = Var(model.i, domain=NonNegativeReals)  # Flow
model.t = Var(model.j, domain=NonNegativeReals)  # Travel time
model.c = Var(model.i, domain=NonNegativeReals)  # Cost
model.r = Var(model.i, domain=NonNegativeReals)  # Residual

model.p1 = Var(domain=NonNegativeReals)
model.p2 = Var(domain=NonNegativeReals)
model.ep1 = Var(domain=NonNegativeReals)
model.ep2 = Var(domain=NonNegativeReals)


# Constraints
def ep_constraints_rule(model, i):
    if i in model.OD1:
        return model.ep1 >= model.r[i]
    elif i in model.OD2:
        return model.ep2 >= model.r[i]
    else:
        return Constraint.Skip


# Flow balance constraints
def balance_constraints_rule(model, i):
    if i in model.OD1:
        return model.f[i] + model.c[i] + model.r[i] - model.p1 >= 0
    elif i in model.OD2:
        return model.f[i] + model.c[i] + model.r[i] - model.p2 >= 0
    else:
        return Constraint.Skip


def path_cost_constraints_rule(model, i):
    if i in model.OD1:
        return model.f[i] * (model.c[i] + model.r[i] - model.p1) == 0
    elif i in model.OD2:
        return model.f[i] * (model.c[i] + model.r[i] - model.p2) == 0
    else:
        return Constraint.Skip


# Travel time constraints
def link_cost_rule(model, j):
    return model.t[j] == model.t0[j] * (1 + 0.15 * (sum(model.A[i, j] * model.f[i] for i in model.i) / 2500) ** 4)


# Path cost constraints
def path_cost_rule(model, i):
    return model.c[i] == sum(model.A[i, j] * model.t[j] for j in model.j)


# Demand constraints
model.demand_constraints = ConstraintList()
model.demand_constraints.add(sum(model.f[i] for i in model.OD1) == 3000)
model.demand_constraints.add(sum(model.f[i] for i in model.OD2) == 3000)
model.ep_constraints = Constraint(model.i, rule=ep_constraints_rule)
model.balance_constraints = Constraint(model.i, rule=balance_constraints_rule)
model.path_cost_constraints = Constraint(model.i, rule=path_cost_constraints_rule)
model.link_cost = Constraint(model.j, rule=link_cost_rule)
model.path_cost = Constraint(model.i, rule=path_cost_rule)

# 添加路径约束条件
model.path_constraints = ConstraintList()
model.path_constraints.add(3000 - model.f[7] >= 0.001)


# 修改目标函数
def objective_rule(model):
    return model.ep1 * 4 + model.ep2 * 7


model.obj = Objective(rule=objective_rule, sense=minimize)

# Solve the model
solver = SolverFactory('ipopt')  # Using IPOPT as the solver
solver.solve(model, tee=True)

# Display results
print("Optimal flow (f):", [model.f[i].value for i in model.i])
print("Optimal travel times (t):", [model.t[j].value for j in model.j])
print("Optimal costs (c):", [model.c[i].value for i in model.i])
print("Objective value (z):", model.obj())
print("r=", [model.r[i].value for i in model.i])
# 显示ep1与ep2的值
print("ep1 =", model.ep1.value)
print("ep2 =", model.ep2.value)
# p1与p2的值
print("p1 =", model.p1.value)
print("p2 =", model.p2.value)
# 将结果显示在表格中
console = Console()
table = Table(title="Results")
table.add_column("OD", style="cyan")
table.add_column("Flow", style="magenta")
table.add_column("Travel time", style="green")
table.add_column("Cost", style="yellow")
table.add_column("Cost+ep", style="blue")
table.add_column("Residual", style="red")
for i in model.i:
    if i in model.OD1:
        cost_plus_ep = model.c[i].value + model.ep1.value
    else:
        cost_plus_ep = model.c[i].value + model.ep2.value
    table.add_row(
        str(i),
        str(model.f[i].value),
        str(model.t[i].value if i in model.j else ""),
        str(model.c[i].value),
        str(cost_plus_ep),
        str(model.r[i].value)
    )

console.print(table)
# 显示目标函数的值
console.print(f"Objective function value: {model.obj()}")
