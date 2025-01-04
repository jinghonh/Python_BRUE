from pyomo.environ import *
from pyomo.opt import SolverFactory

# 创建模型
model = ConcreteModel()

# 定义变量
model.x1 = Var(domain=NonNegativeReals)
model.x2 = Var(domain=NonNegativeReals)
model.x3 = Var(domain=NonNegativeReals)
model.x4 = Var(domain=NonNegativeReals)
model.x5 = Var(domain=NonNegativeReals)
model.x6 = Var(domain=NonNegativeReals)
model.x7 = Var(domain=NonNegativeReals)
model.x8 = Var(domain=NonNegativeReals)
model.x9 = Var(domain=NonNegativeReals)
model.x10 = Var(domain=NonNegativeReals)
model.x11 = Var(domain=NonNegativeReals)
model.x12 = Var(domain=NonNegativeReals)
model.x13 = Var(domain=NonNegativeReals)
model.x14 = Var(domain=Reals)
model.c1 = Var(domain=NonNegativeReals)
model.c2 = Var(domain=NonNegativeReals)
model.c3 = Var(domain=NonNegativeReals)
model.c4 = Var(domain=NonNegativeReals)
model.c5 = Var(domain=NonNegativeReals)
model.c6 = Var(domain=NonNegativeReals)

# 定义约束
model.e1 = Constraint(expr=model.x1 + model.x2 + model.x3 + model.x4 + model.x5 + model.x6 == 10000)
model.e2 = Constraint(expr=10000 - model.x1 - model.x2 - model.x5 >= 0.01)

# x14 约束
model.e3 = Constraint(expr=model.x14 - model.x7 >= 0)
model.e4 = Constraint(expr=model.x14 - model.x8 >= 0)
model.e5 = Constraint(expr=model.x14 - model.x9 >= 0)
model.e6 = Constraint(expr=model.x14 - model.x10 >= 0)
model.e7 = Constraint(expr=model.x14 - model.x11 >= 0)
model.e8 = Constraint(expr=model.x14 - model.x12 >= 0)

# 流量平衡约束
model.e9 = Constraint(expr=model.x1 + model.c1 + model.x7 - model.x13 >= 0.0001)
model.e10 = Constraint(expr=model.x2 + model.c2 + model.x8 - model.x13 >= 0.0001)
model.e11 = Constraint(expr=model.x3 + model.c3 + model.x9 - model.x13 >= 0.0001)
model.e12 = Constraint(expr=model.x4 + model.c4 + model.x10 - model.x13 >= 0.0001)
model.e13 = Constraint(expr=model.x5 + model.c5 + model.x11 - model.x13 >= 0.0001)
model.e14 = Constraint(expr=model.x6 + model.c6 + model.x12 - model.x13 >= 0.0001)

# 路径成本约束
model.e15 = Constraint(expr=model.c1 + model.x7 - model.x13 >= 0)
model.e16 = Constraint(expr=model.c2 + model.x8 - model.x13 >= 0)
model.e17 = Constraint(expr=model.c3 + model.x9 - model.x13 >= 0)
model.e18 = Constraint(expr=model.c4 + model.x10 - model.x13 >= 0)
model.e19 = Constraint(expr=model.c5 + model.x11 - model.x13 >= 0)
model.e20 = Constraint(expr=model.c6 + model.x12 - model.x13 >= 0)

# 定义目标函数
model.objective = Objective(expr=model.x14, sense=minimize)
# e21：定义目标函数的权重平衡约束
model.e21 = Constraint(expr=(
    model.x1 * (model.c1 + model.x7 - model.x13) +
    model.x2 * (model.c2 + model.x8 - model.x13) +
    model.x3 * (model.c3 + model.x9 - model.x13) +
    model.x4 * (model.c4 + model.x10 - model.x13) +
    model.x5 * (model.c5 + model.x11 - model.x13) +
    model.x6 * (model.c6 + model.x12 - model.x13)
) == 0)
# 路径成本的非线性约束
model.e22 = Constraint(expr=model.c1 == 18 * 1 + 18 * 0.15 * (model.x1 / 3600)**4 + 
                        15 * (1 - exp(-0.02 * (18 * 1 + 18 * 0.15 * (model.x1 / 3600)**4))))
model.e23 = Constraint(expr=model.c2 == 22.5 * 1 + 22.5 * 0.15 * (model.x2 / 3600)**4 + 
                        15 * (1 - exp(-0.02 * (22.5 * 1 + 22.5 * 0.15 * (model.x2 / 3600)**4))))
model.e24 = Constraint(expr=model.c3 == 12 * 1 + 12 * 0.15 * ((model.x3 + model.x5) / 1800)**4 + 24 +
                        24 * 0.15 * ((model.x3 + model.x6) / 1800)**4 + 
                        15 * (1 - exp(-0.02 * (12 * 1 + 12 * 0.15 * ((model.x3 + model.x5) / 1800)**4 + 
                                              24 * 1 + 24 * 0.15 * ((model.x3 + model.x6) / 1800)**4))))
model.e25 = Constraint(expr=model.c4 == 24 * 1 + 24 * 0.15 * ((model.x4 + model.x6) / 1800)**4 + 12 +
                        12 * 0.15 * ((model.x4 + model.x5) / 1800)**4 + 
                        15 * (1 - exp(-0.02 * (24 * 1 + 24 * 0.15 * ((model.x4 + model.x6) / 1800)**4 + 
                                              12 * 1 + 12 * 0.15 * ((model.x4 + model.x5) / 1800)**4))))
model.e26 = Constraint(expr=model.c5 == 12 * 1 + 12 * 0.15 * ((model.x3 + model.x5) / 1800)**4 + 
                        2.4 * 1 + 2.4 * 0.15 * (model.x5 / 1800)**4 + 
                        12 * 1 + 12 * 0.15 * ((model.x4 + model.x5) / 1800)**4 + 
                        15 * (1 - exp(-0.02 * (12 * 1 + 12 * 0.15 * ((model.x3 + model.x5) / 1800)**4 + 
                                              2.4 * 1 + 2.4 * 0.15 * (model.x5 / 1800)**4 + 
                                              12 * 1 + 12 * 0.15 * ((model.x4 + model.x5) / 1800)**4))))
model.e27 = Constraint(expr=model.c6 == 24 * 1 + 24 * 0.15 * ((model.x4 + model.x6) / 1800)**4 + 
                        6 * 1 + 6 * 0.15 * (model.x6 / 1800)**4 + 
                        24 * 1 + 24 * 0.15 * ((model.x3 + model.x6) / 1800)**4 + 
                        15 * (1 - exp(-0.02 * (24 * 1 + 24 * 0.15 * ((model.x4 + model.x6) / 1800)**4 + 
                                              6 * 1 + 6 * 0.15 * (model.x6 / 1800)**4 + 
                                              24 * 1 + 24 * 0.15 * ((model.x3 + model.x6) / 1800)**4))))

# 求解
solver = SolverFactory('ipopt')
results = solver.solve(model, tee=True)

# 显示结果
print("x1 =", model.x1.value)
print("x2 =", model.x2.value)
print("x3 =", model.x3.value)
print("x4 =", model.x4.value)
print("x5 =", model.x5.value)
print("x6 =", model.x6.value)
print("x14 =", model.x14.value)
print("c1 =", model.c1.value)
print("c2 =", model.c2.value)
print("c3 =", model.c3.value)
print("c4 =", model.c4.value)
print("c5 =", model.c5.value)
print("c6 =", model.c6.value)
