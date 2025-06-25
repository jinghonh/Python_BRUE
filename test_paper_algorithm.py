# test_paper_algorithm.py - 严格按照论文算法实现的多OD对BRUE求解器

from traffic_network_config import TrafficNetworkConfig
from pyomo.environ import *
from rich.console import Console
from rich.table import Table
import copy

class PaperAlgorithmBRUESolver:
    def __init__(self, config: TrafficNetworkConfig):
        self.config = config
        self.console = Console()

    def _solve_ue(self):
        """求解用户均衡作为起点"""
        model = ConcreteModel()

        # 集合和参数
        model.od_pairs = Set(initialize=list(range(1, self.config.num_od_pairs + 1)))
        model.links = Set(initialize=list(range(1, self.config.num_paths + 1)))
        model.group_names = Set(initialize=self.config.od_groups.keys())

        model.free_flow_time = Param(model.links, initialize=self.config.free_flow_time)
        model.path_link_matrix = Param(model.od_pairs, model.links, 
                                     initialize=self.config.path_link_matrix, default=0)
        model.od_demands = Param(model.group_names, initialize=self.config.od_demands)

        # 变量
        model.flow = Var(model.od_pairs, domain=NonNegativeReals, initialize=0)
        model.link_flow = Var(model.links, domain=NonNegativeReals)
        model.travel_time = Var(model.links, domain=NonNegativeReals)
        model.path_cost = Var(model.od_pairs, domain=NonNegativeReals)

        # 约束
        def link_flow_rule(m, l):
            return m.link_flow[l] == sum(m.path_link_matrix[p, l] * m.flow[p] for p in m.od_pairs)
        model.link_flow_constraint = Constraint(model.links, rule=link_flow_rule)
        
        def travel_time_rule(m, l):
            # return m.travel_time[l] == m.free_flow_time[l] + m.link_flow[l]
            return m.travel_time[l] == m.free_flow_time[l] * (
                1 + 0.15 * (m.link_flow[l] / self.config.link_capacity[l]) ** 4
            )
        model.travel_time_constraint = Constraint(model.links, rule=travel_time_rule)

        def path_cost_rule(m, p):
            return m.path_cost[p] == sum(m.path_link_matrix[p, l] * m.travel_time[l] for l in m.links)
        model.path_cost_constraint = Constraint(model.od_pairs, rule=path_cost_rule)

        # 需求约束
        model.demand_constraint = ConstraintList()
        for g, paths in self.config.od_groups.items():
            model.demand_constraint.add(sum(model.flow[p] for p in paths) == model.od_demands[g])

        # UE目标函数
        model.objective = Objective(
            expr=sum(model.path_cost[p] * model.flow[p] for p in model.od_pairs), 
            sense=minimize
        )

        # 求解
        solver = SolverFactory('ipopt')
        status = solver.solve(model, tee=False)
        
        if status.solver.termination_condition != TerminationCondition.optimal:
            raise Exception("UE求解失败")

        # 提取流量
        flows = {}
        for p in model.od_pairs:
            flows[p] = model.flow[p].value if model.flow[p].value else 0
        
        past_cost = {}
        for p in model.od_pairs:
            past_cost[p] = model.path_cost[p].value if model.path_cost[p].value else 0
        print(past_cost)
        return flows

    def _solve_equation_19(self, target_od_group: str, other_flows: dict, 
                           acceptable_paths: dict):
        """
        实现论文公式(19) - 为特定OD对计算临界点
        """
        model = ConcreteModel()
        
        target_paths = self.config.od_groups[target_od_group]
        current_acceptable = acceptable_paths.get(target_od_group, [])
        potential_new_paths = [p for p in target_paths if p not in current_acceptable]
        
        if not potential_new_paths:
            return None, None, None
        
        self.console.print(f"OD组 '{target_od_group}': 已接受路径 {current_acceptable}, 潜在新路径 {potential_new_paths}")
        
        # 集合
        model.target_paths = Set(initialize=target_paths)
        model.links = Set(initialize=list(range(1, self.config.num_paths + 1)))
        model.od_pairs = Set(initialize=list(range(1, self.config.num_od_pairs + 1)))
        
        # 参数
        model.free_flow_time = Param(model.links, initialize=self.config.free_flow_time)
        model.path_link_matrix = Param(model.od_pairs, model.links, 
                                     initialize=self.config.path_link_matrix, default=0)
        model.target_demand = Param(initialize=self.config.od_demands[target_od_group])
        
        # 变量
        model.target_flow = Var(model.target_paths, domain=NonNegativeReals)
        model.link_flow = Var(model.links, domain=NonNegativeReals)
        model.travel_time = Var(model.links, domain=NonNegativeReals)
        model.path_cost = Var(model.od_pairs, domain=NonNegativeReals)
        
        # 论文公式(19)的关键变量
        model.epsilon = Var(domain=NonNegativeReals)  # 目标：最小化这个ε
        model.pi = Var(domain=NonNegativeReals)       # π^ν
        model.rho = Var(model.target_paths, domain=NonNegativeReals)  # ρ_i^ν
        
        # 目标函数：min ε (论文公式19第一行)
        model.objective = Objective(expr=model.epsilon, sense=minimize)
        
        # 链路流量：目标OD组流量 + 其他固定流量
        def link_flow_rule(m, l):
            target_contrib = sum(m.path_link_matrix[p, l] * m.target_flow[p] 
                               for p in m.target_paths if p in m.od_pairs)
            other_contrib = sum(m.path_link_matrix[p, l] * other_flows.get(p, 0) 
                              for p in m.od_pairs if p not in m.target_paths)
            return m.link_flow[l] == target_contrib + other_contrib
        model.link_flow_constraint = Constraint(model.links, rule=link_flow_rule)
        
        # 旅行时间 cost = to + f
        # def travel_time_rule(m, l):
        #     return m.travel_time[l] == m.free_flow_time[l] + m.link_flow[l]
        # model.travel_time_constraint = Constraint(model.links, rule=travel_time_rule)

        # 另一种非线性旅行时间 cost = t0（1 + 0.15(f/c)^4）
        def travel_time_rule_nonlinear(m, l):
            return m.travel_time[l] == m.free_flow_time[l] * (
                1 + 0.15 * (m.link_flow[l] / self.config.link_capacity[l]) ** 4
            )
        model.travel_time_constraint_nonlinear = Constraint(model.links, rule=travel_time_rule_nonlinear)
        
        # 路径成本
        def path_cost_rule(m, p):
            return m.path_cost[p] == sum(m.path_link_matrix[p, l] * m.travel_time[l] 
                                       for l in m.links)
        model.path_cost_constraint = Constraint(model.od_pairs, rule=path_cost_rule)
        
        # 需求约束：∑f_j^ν = d^ν (论文公式19e第三行)
        model.demand_constraint = Constraint(
            expr=sum(model.target_flow[p] for p in model.target_paths) == model.target_demand
        )
        
        # BRUE条件：所有路径成本 + ρ ≤ π
        model.brue_constraints = ConstraintList()
        for p in model.target_paths:
            model.brue_constraints.add(model.path_cost[p] + model.rho[p] >= model.pi)
        
        model.total_cost_constraint = Constraint(
                expr=sum(model.target_flow[i] * (model.path_cost[i] + 
                    model.rho[i] - model.pi) for i in model.target_paths) == 0
            )
        
        # ρ的边界：0 ≤ ρ_i^ν ≤ ε (论文公式19第四行)
        model.rho_bounds = ConstraintList()
        for p in model.target_paths:
            model.rho_bounds.add(model.rho[p] <= model.epsilon)
        
        # 强制使用新路径：d^ν - ∑f_j^{j-1} > 0 (论文公式19e第三行变形)
        # 即强制至少有δ流量分配给潜在新路径
        delta = 1e-4
        model.force_new_path = Constraint(
            expr=sum(model.target_flow[p] for p in potential_new_paths) >= delta
        )
        
        # # 为了确保π的准确性，添加最小成本约束
        # model.min_cost_var = Var(domain=NonNegativeReals)
        # model.min_cost_constraints = ConstraintList()
        # for p in model.target_paths:
        #     model.min_cost_constraints.add(model.min_cost_var <= model.path_cost[p])
        # # π应该等于最小成本 + ε（简化的BRUE条件）
        # model.pi_definition = Constraint(expr=model.pi == model.min_cost_var + model.epsilon)
        
        # 求解
        solver = SolverFactory('ipopt')
        status = solver.solve(model, tee=False)
        
        if status.solver.termination_condition != TerminationCondition.optimal:
            return None, None, None
        
        # 提取结果
        critical_epsilon = model.epsilon.value
        new_paths = []
        target_flows = {}
        
        for p in model.target_paths:
            flow_val = model.target_flow[p].value if model.target_flow[p].value else 0
            target_flows[p] = flow_val
            if p in potential_new_paths and flow_val > delta / 2:
                new_paths.append(p)
        
        return critical_epsilon, new_paths, target_flows

    def solve_multi_od_brue(self):
        """
        按照论文算法求解多OD对BRUE问题
        """
        self.console.print("[bold green]========== 论文算法：多OD对BRUE求解 ===========[/bold green]")
        
        # 步骤0：求解UE获得初始流量
        self.console.print("[bold]步骤0: 求解用户均衡(UE)...[/bold]")
        initial_flows = self._solve_ue()
        
        # 初始化可接受路径集合（UE中有流量的路径）
        acceptable_paths = {g: [] for g in self.config.od_groups.keys()}
        critical_points = {g: [0.0] for g in self.config.od_groups.keys()}
        # 用于跟踪每个临界点的路径集合
        path_evolution = {g: [] for g in self.config.od_groups.keys()}
        
        for g, paths in self.config.od_groups.items():
            for p in paths:
                if initial_flows.get(p, 0) > 1e-5:
                    acceptable_paths[g].append(p)
            # 记录初始可接受路径
            path_evolution[g].append(acceptable_paths[g].copy())
        
        self.console.print(f"UE解下的可接受路径: {acceptable_paths}")
        
        # 主迭代循环
        current_flows = copy.deepcopy(initial_flows)
        max_iterations = 10
        
        for iteration in range(1, max_iterations + 1):
            self.console.print(f"\n[bold green]========== 迭代 {iteration} ===========[/bold green]")
            found_improvement = False
            
            # 论文算法步骤1：对每个OD对计算临界点
            for target_od_group in self.config.od_groups.keys():
                self.console.print(f"\n--- 处理OD组 '{target_od_group}' ---")
                
                # 准备其他OD组的固定流量（论文：f^w, w≠ν作为参数）
                other_flows = {}
                for p in range(1, self.config.num_od_pairs + 1):
                    if p not in self.config.od_groups[target_od_group]:
                        other_flows[p] = current_flows.get(p, 0)
                
                # 使用公式(19)计算临界点
                critical_epsilon, new_paths, target_flows = self._solve_equation_19(
                    target_od_group, other_flows, acceptable_paths
                )
                self.console.print(f"critical_epsilon: {critical_epsilon}")
                
                if critical_epsilon is not None and new_paths:
                    last_epsilon = critical_points[target_od_group][-1]
                    
                    if abs(critical_epsilon - last_epsilon) > 1e-6:
                        self.console.print(f"[green]找到新临界点! ε* = {critical_epsilon:.6f}[/green]")
                        self.console.print(f"新增路径: {new_paths}")
                        
                        # 更新临界点和可接受路径
                        critical_points[target_od_group].append(critical_epsilon)
                        acceptable_paths[target_od_group].extend(new_paths)
                        acceptable_paths[target_od_group] = sorted(list(set(acceptable_paths[target_od_group])))
                        
                        # 记录当前临界点的路径集合
                        path_evolution[target_od_group].append(acceptable_paths[target_od_group].copy())
                        
                        # 更新流量（论文算法步骤2的副产品）
                        for p, flow in target_flows.items():
                            current_flows[p] = flow
                        
                        found_improvement = True
                    else:
                        self.console.print(f"[yellow]临界点无改进: {critical_epsilon:.6f} ≤ {last_epsilon:.6f}[/yellow]")
                else:
                    self.console.print(f"[yellow]OD组 '{target_od_group}' 无新路径可加入[/yellow]")
            
            if not found_improvement:
                self.console.print("\n[bold green]算法收敛：所有OD组都无法找到更大的临界点[/bold green]")
                break
        
        # 显示结果
        self._display_results(critical_points, path_evolution, initial_flows)

    def _display_results(self, critical_points, path_evolution, initial_flows):
        """显示最终结果"""
        self.console.print("\n[bold underline]========== 最终结果 ==========[/bold underline]")
        
        for g in self.config.od_groups.keys():
            table = Table(title=f"OD组 '{g}' 的临界点演化")
            table.add_column("阶段", style="cyan")
            table.add_column("临界点 ε*", style="yellow")
            table.add_column("可接受路径集合", style="green")
            
            # UE阶段 - 使用最初的可接受路径
            table.add_row("UE", f"{critical_points[g][0]:.6f}", str(sorted(path_evolution[g][0])))
            
            # 临界点阶段 - 循环显示每个临界点及其对应的可接受路径
            for i, eps in enumerate(critical_points[g][1:], 1):
                table.add_row(f"临界点{i}", f"{eps:.6f}", str(sorted(path_evolution[g][i])))
            
            self.console.print(table)

def main():
    # 测试两OD网络
    config = TrafficNetworkConfig.create_multi_od_network()
    solver = PaperAlgorithmBRUESolver(config)
    solver.solve_multi_od_brue()

if __name__ == "__main__":
    main()
