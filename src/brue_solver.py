from .brue_base import BRUEBase
from .traffic_network_config import TrafficNetworkConfig
from pyomo.environ import *
import matplotlib.pyplot as plt
from rich.console import Console
from rich.table import Table


class BRUESolver(BRUEBase):
    def __init__(self, config: TrafficNetworkConfig):
        super().__init__()
        self.config = config

    def initialize_sets(self):
        self.model.od_pairs = RangeSet(1, self.config.num_paths)
        self.model.paths = RangeSet(1, self.config.num_links)
        self.model.od_demand = Set(initialize=range(1, self.config.num_paths + 1))
        # 初始化OD组
        for group_name, group_pairs in self.config.od_groups.items():
            setattr(self.model, f'OD_{group_name}', Set(initialize=group_pairs))

    def initialize_parameters(self):
        # 初始化所有参数
        self.model.free_flow_time = Param(
            self.model.paths,
            initialize=self.config.free_flow_time
        )
        self.model.link_money_cost = Param(
            self.model.paths,
            initialize=self.config.link_money_cost
        )
        self.model.link_capacity = Param(
            self.model.paths,
            initialize=self.config.link_capacity
        )
        self.model.path_link_matrix = Param(
            self.model.od_pairs,
            self.model.paths,
            initialize=self.config.path_link_matrix,
            default=0
        )

    def initialize_variables(self):
        m = self.model
        m.flow = Var(m.od_pairs, domain=NonNegativeReals)
        m.travel_time = Var(m.paths, domain=NonNegativeReals)
        m.path_cost = Var(m.od_pairs, domain=NonNegativeReals)
        m.residual = Var(m.od_pairs, domain=NonNegativeReals)

        # 修改perception的定义，为每个OD组创建独立的perception
        if len(self.config.od_groups) > 1:
            m.perceptions = Var(self.config.od_groups.keys(), domain=NonNegativeReals)
            m.epsilons = Var(self.config.od_groups.keys(), domain=NonNegativeReals)
        else:
            m.perception = Var(domain=NonNegativeReals)
            m.epsilon = Var(domain=NonNegativeReals)

    def add_constraints(self):
        m = self.model
        # 总需求约束 - 修改为使用每个OD组独立的需求
        m.demand_constraint = ConstraintList()
        for group_name, group_pairs in self.config.od_groups.items():
            m.demand_constraint.add(
                sum(m.flow[i] for i in group_pairs) == self.config.od_demands[group_name]
            )

        # 误差约束
        def epsilon_rule(m, i):
            for group_name, group_pairs in self.config.od_groups.items():
                if i in group_pairs:
                    if len(self.config.od_groups) > 1:
                        return m.epsilons[group_name] >= m.residual[i]
                    else:
                        return m.epsilon >= m.residual[i]
            return Constraint.Skip

        m.epsilon_constraints = Constraint(m.od_pairs, rule=epsilon_rule)

        # 平衡约束 - 修改为使用对应组的perception
        def balance_rule(m, i):
            for group_name, group_pairs in self.config.od_groups.items():
                if i in group_pairs:
                    if len(self.config.od_groups) > 1:
                        return m.flow[i] + m.path_cost[i] + m.residual[i] - m.perceptions[group_name] >= 0.0
                    else:
                        return m.flow[i] + m.path_cost[i] + m.residual[i] - m.perception >= 0.0
            return Constraint.Skip

        m.balance_constraints = Constraint(m.od_pairs, rule=balance_rule)

        # 路径成本非负约束 - 修改为使用对应组的perception
        def cost_nonnegativity_rule(m, i):
            for group_name, group_pairs in self.config.od_groups.items():
                if i in group_pairs:
                    if len(self.config.od_groups) > 1:
                        return m.path_cost[i] + m.residual[i] - m.perceptions[group_name] >= 0
                    else:
                        return m.path_cost[i] + m.residual[i] - m.perception >= 0
            return Constraint.Skip

        m.path_cost_nonnegativity = Constraint(m.od_pairs, rule=cost_nonnegativity_rule)

        # 总成本约束 - 修改为使用对应组的perception
        m.total_cost_constraint = ConstraintList()
        if len(self.config.od_groups) > 1:
            for group_name, group_pairs in self.config.od_groups.items():
                m.total_cost_constraint.add(
                    sum(m.flow[i] * (m.path_cost[i] + m.residual[i] - m.perceptions[group_name])
                        for i in group_pairs) == 0
                )
        else:
            m.total_cost_constraint.add(
                sum(m.flow[i] * (m.path_cost[i] + m.residual[i] - m.perception)
                    for i in m.od_pairs) == 0
            )

        # 路径成本约束
        def travel_time_rule(m, j):
            return m.travel_time[j] == m.free_flow_time[j] * (
                    1 + 0.15 * (sum(m.path_link_matrix[i, j] * m.flow[i]
                                    for i in m.od_pairs) / m.link_capacity[j]) ** 4
            )

        m.travel_time_constraints = Constraint(m.paths, rule=travel_time_rule)

        # 路径成本计算约束
        def path_cost_rule(m, i):
            return m.path_cost[i] == sum(
                m.path_link_matrix[i, j] * m.travel_time[j] for j in m.paths
            )

        m.path_cost_constraints = Constraint(m.od_pairs, rule=path_cost_rule)

    def set_objective(self):
        if len(self.config.od_groups) > 1:
            # 多OD组情况下的目标函数
            weights = {'OD1': 4, 'OD2': 7}  # 可以根据需要调整权重
            self.model.objective = Objective(
                expr=sum(self.model.epsilons[g] * weights.get(g, 1)
                         for g in self.config.od_groups.keys()),
                sense=minimize
            )
        else:
            # 单一OD组情况下的目标函数
            self.model.objective = Objective(
                expr=self.model.epsilon,
                sense=minimize
            )

    def calculate_money_cost(self):
        """计算每条路径的金钱成本"""
        m = self.model
        money_costs = {}
        for i in m.od_pairs:
            money_costs[i] = sum(m.path_link_matrix[i, j] * m.link_money_cost[j]
                                 for j in m.paths)
        return money_costs

    @staticmethod
    def is_dominated(path_i, path_costs, money_costs, all_paths):
        """判断一个路径是否被其他路径支配"""
        cost_i = path_costs[path_i]
        money_i = money_costs[path_i]

        for path_j in all_paths:
            if path_j != path_i:
                cost_j = path_costs[path_j]
                money_j = money_costs[path_j]

                if (cost_j <= cost_i and money_j <= money_i and
                        (cost_j < cost_i or money_j < money_i)):
                    return True
        return False

    def display_results(self):
        """显示计算结果"""
        table = Table(title="计算结果")

        columns = ["路径", "流量", "时间成本", "金钱成本", "残差", "感知值"]
        for col in columns:
            table.add_column(col, justify="right")

        for group_name, group_pairs in self.config.od_groups.items():
            perception = (self.model.perceptions[group_name].value
                          if len(self.config.od_groups) > 1
                          else self.model.perception.value)

            for i in group_pairs:
                money_cost = sum(self.model.path_link_matrix[i, j] *
                                 self.model.link_money_cost[j] for j in self.model.paths)
                table.add_row(
                    f"{i} ({group_name})",
                    f"{self.model.flow[i].value:.2f}",
                    f"{self.model.path_cost[i].value:.2f}",
                    f"{money_cost:.2f}",
                    f"{self.model.residual[i].value:.2f}",
                    f"{perception:.2f}"
                )

        self.console.print(table)

    def analyze_path_costs(self):
        """分析路径成本并找出有效路径（向后兼容）"""
        _, all_effective_paths = self.analyze_path_costs_by_group()
        return all_effective_paths

    def plot_cost_analysis(self, is_effective_path, iteration_data=None):
        """Plot path cost analysis with money cost"""
        m = self.model
        money_costs = self.calculate_money_cost()
        path_costs = {i: m.path_cost[i].value for i in m.od_pairs}

        # Analyze each OD group
        for group_name, group_pairs in self.config.od_groups.items():
            if len(self.config.od_groups) > 1:
                epsilon = m.epsilons[group_name].value
            else:
                epsilon = m.epsilon.value
            min_cost = min(path_costs[i] for i in group_pairs)
            upper_bound = min_cost + epsilon + 2e-2

            # 修改判断可行路径的逻辑
            feasible_paths = []
            infeasible_paths = []
            for i in group_pairs:
                cost = path_costs[i]
                # 判断是否在可行范围内
                if abs(cost - min_cost) <= 1e-2:  # 最小成本路径
                    feasible_paths.append(i)
                elif cost <= upper_bound + 1e-2:  # 在上界范围内的路径
                    feasible_paths.append(i)
                else:  # 超出上界的路径
                    infeasible_paths.append(i)

            # 在可行路径中判断支配关系
            dominated_paths = []
            non_dominated_paths = []
            for path in feasible_paths:
                if not self.is_dominated(path, path_costs, money_costs, feasible_paths):
                    non_dominated_paths.append(path)
                else:
                    dominated_paths.append(path)

            # Create plot
            plt.figure(figsize=(10, 6))

            # Plot points
            if infeasible_paths:
                plt.scatter([path_costs[i] for i in infeasible_paths],
                            [money_costs[i] for i in infeasible_paths],
                            color='grey', alpha=0.5, label='Infeasible Paths')
            if dominated_paths:
                plt.scatter([path_costs[i] for i in dominated_paths],
                            [money_costs[i] for i in dominated_paths],
                            color='red', alpha=0.7, label='Dominated Paths')
            if non_dominated_paths:
                plt.scatter([path_costs[i] for i in non_dominated_paths],
                            [money_costs[i] for i in non_dominated_paths],
                            color='green', alpha=1.0, label='Non-dominated Paths')

            # 添加带有调整位置的标签
            label_positions = {}  # 用于存储已使用的标签位置
            for i in group_pairs:
                x, y = path_costs[i], money_costs[i]

                # 计算标签位置，避免重叠
                offset_x = 5
                offset_y = 5
                label_x = x
                label_y = y

                # 检查附近是否有其他标签
                while any(abs(label_x - lx) < 10 and abs(label_y - ly) < 10
                          for lx, ly in label_positions.values()):
                    offset_y += 5  # 逐渐增加垂直偏移
                    if offset_y > 30:  # 如果垂直偏移太大，尝试水平偏移
                        offset_y = 5
                        offset_x += 10
                    label_x = x + offset_x
                    label_y = y + offset_y

                # 存储这个标签位置
                label_positions[i] = (label_x, label_y)

                # 绘制带有连接线的标签
                plt.annotate(
                    f'P{i}\n(t:{path_costs[i]:.1f}, m:{money_costs[i]:.1f})',
                    xy=(x, y),  # 点的位置
                    xytext=(label_x, label_y),  # 标签的位置
                    textcoords='offset points',
                    fontsize=8,
                    bbox=dict(facecolor='white', edgecolor='none', alpha=0.7),
                    arrowprops=dict(arrowstyle='-', color='gray', alpha=0.5)
                )

            # 设置图表属性
            if iteration_data:
                # 处理迭代数据中的epsilon可能是字典的情况
                if isinstance(iteration_data["epsilon"], dict):
                    group_epsilon = iteration_data["epsilon"].get(group_name, epsilon)
                    epsilon_display = f"{group_epsilon:.2f}"
                else:
                    epsilon_display = f"{epsilon:.2f}"
                    
                plt.title(f'Path Cost Analysis - {group_name}\n'
                          f'Iteration {iteration_data["iteration"]}, '
                          f'Restricted Paths: {iteration_data["restricted_paths"]}\n'
                          f'ε = {epsilon_display}')
            else:
                plt.title(f'Path Cost Analysis - {group_name}\nε = {epsilon:.2f}')

            plt.xlabel('Travel Time Cost')
            plt.ylabel('Money Cost')
            plt.grid(True, alpha=0.3)

            # Add bound lines
            plt.axvline(x=min_cost, color='blue', linestyle='--', alpha=0.5,
                        label=f'Min Cost ({min_cost:.2f})')
            plt.axvline(x=upper_bound, color='red', linestyle='--', alpha=0.5,
                        label=f'Upper Bound ({upper_bound:.2f})')

            # Add info text
            info_text = (
                f'Feasible Paths: {sorted(feasible_paths)}\n'
                f'Non-dominated Paths: {sorted(non_dominated_paths)}\n'
                f'Dominated Paths: {sorted(dominated_paths)}\n'
                f'Infeasible Paths: {sorted(infeasible_paths)}'
            )
            plt.text(0.98, 0.98, info_text,
                     transform=plt.gca().transAxes,
                     verticalalignment='top',
                     horizontalalignment='right',
                     bbox=dict(boxstyle='round', facecolor='white', alpha=0.8),
                     fontsize=8)

            plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
            plt.tight_layout()
            plt.show()

        # 返回适当的epsilon值，对于多OD组情况应返回字典
        if len(self.config.od_groups) > 1:
            epsilon_return = {g: m.epsilons[g].value for g in self.config.od_groups.keys()}
        else:
            epsilon_return = m.epsilon.value
            
        return {
            'effective_paths': is_effective_path,
            'path_costs': path_costs,
            'money_costs': money_costs,
            'epsilon': epsilon_return
        }

    def plot_initial_costs(self):
        """Plot three scenarios of path costs"""
        m = self.model
        fig, axes = plt.subplots(1, 3, figsize=(18, 6))

        # Calculate costs for three scenarios
        scenarios = {
            'Initial State (t0)': {
                i: sum(m.path_link_matrix[i, j] * m.free_flow_time[j]
                       for j in m.paths) for i in m.od_pairs
            },
            'No Flow State': {
                i: sum(m.path_link_matrix[i, j] * m.free_flow_time[j]
                       for j in m.paths) for i in m.od_pairs
            },
            'Full Capacity State': {
                i: sum(m.path_link_matrix[i, j] * m.free_flow_time[j] *
                       (1 + 0.15 * (1.0) ** 4)
                       for j in m.paths) for i in m.od_pairs
            }
        }

        money_costs = self.calculate_money_cost()

        # Plot each scenario
        for idx, (scenario_name, time_costs) in enumerate(scenarios.items()):
            ax = axes[idx]

            # Analyze dominance
            dominated_paths = []
            non_dominated_paths = []

            for path_i in m.od_pairs:
                if self.is_dominated(path_i, time_costs, money_costs, m.od_pairs):
                    dominated_paths.append(path_i)
                else:
                    non_dominated_paths.append(path_i)

            # Plot points
            if dominated_paths:
                ax.scatter([time_costs[i] for i in dominated_paths],
                           [money_costs[i] for i in dominated_paths],
                           color='red', alpha=0.7, label='Dominated')
            if non_dominated_paths:
                ax.scatter([time_costs[i] for i in non_dominated_paths],
                           [money_costs[i] for i in non_dominated_paths],
                           color='green', alpha=1.0, label='Non-dominated')

            # Add labels
            for i in m.od_pairs:
                ax.annotate(f'P{i}\n(t:{time_costs[i]:.1f}, m:{money_costs[i]:.1f})',
                            (time_costs[i], money_costs[i]),
                            xytext=(5, 5), textcoords='offset points',
                            fontsize=8)

            # Set properties
            ax.set_xlabel('Travel Time Cost')
            ax.set_ylabel('Money Cost')
            ax.set_title(f'{scenario_name}\n'
                         f'Non-dominated: {sorted(non_dominated_paths)}\n'
                         f'Dominated: {sorted(dominated_paths)}')
            ax.legend()
            ax.grid(True, alpha=0.3)

        plt.suptitle('Path Cost Analysis Under Different Scenarios', fontsize=14)
        plt.tight_layout()
        plt.show()

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

    def add_path_constraints(self, restricted_paths=None, prev_epsilon=None):
        """添加路径约束"""
        m = self.model

        # 添加路径流量约束
        m.path_constraint = ConstraintList()
        for group_name, group_pairs in self.config.od_groups.items():
            group_restricted_paths = [p for p in restricted_paths if p in group_pairs] if restricted_paths else []
            # 修改约束条件，使其更加宽松
            if group_restricted_paths:  # 只在有受限路径时添加约束
                # 修改约束：从严格限制改为更宽松的上限
                m.path_constraint.add(
                    sum(m.flow[i] for i in group_restricted_paths) <=
                    self.config.od_demands[group_name] * 0.99
                )

        # 如果有前一次迭代的epsilon，添加epsilon下界约束
        if prev_epsilon is not None:
            if len(self.config.od_groups) > 1:
                for group_name in self.config.od_groups.keys():
                    m.path_constraint.add(m.epsilons[group_name] >= prev_epsilon[group_name])
            else:
                m.path_constraint.add(m.epsilon >= prev_epsilon)

    def run_with_iterations(self, initial_paths=None):
        """
        使用迭代方法求解，每个OD组有独立的路径和epsilon
        Args:
            initial_paths: 初始路径字典，格式为{od_group_name: [paths]}，默认为None
        Returns:
            list of dict: 迭代结果列表
        """
        results = []
        # 初始化每个OD组的路径集合
        current_paths = {}
        if initial_paths is not None:
            current_paths = initial_paths
        else:
            for group_name in self.config.od_groups.keys():
                current_paths[group_name] = []
        
        # 每个OD组独立的epsilon
        prev_epsilon = {group_name: None for group_name in self.config.od_groups.keys()}
        iteration = 1
        max_iterations = 10  # 添加最大迭代次数限制
        relax_factor = 0.99  # 初始放松因子
        min_relax_factor = 0.5  # 最小放松因子

        while iteration <= max_iterations:
            # 初始化新的模型
            self.model = ConcreteModel()
            self.initialize_sets()
            self.initialize_parameters()
            self.initialize_variables()
            self.add_constraints()

            # 添加当前迭代的路径约束 - 需要修改add_path_constraints以支持按OD组约束
            self.add_path_constraints_by_group(current_paths, prev_epsilon, relax_factor)

            # 设置目标函数并求解
            self.set_objective()

            # tee 为 False 时，不显示求解日志
            solve_status = self.solve(tee=False)

            if solve_status.solver.status != SolverStatus.ok:
                self.console.print(f"[red]迭代 {iteration} 求解失败，放松因子为 {relax_factor}[/red]")
                
                # 如果有结果，先尝试放松约束重新求解
                if relax_factor > min_relax_factor:
                    relax_factor -= 0.1  # 每次减少放松因子
                    self.console.print(f"[yellow]降低放松因子至 {relax_factor:.2f} 重新尝试...[/yellow]")
                    continue  # 重新尝试当前迭代，不增加迭代计数
                
                # 如果放松到最低限度还是不行，考虑独立求解各OD组
                self.console.print("[yellow]尝试为各OD组独立求解...[/yellow]")
                
                # 尝试独立求解每个OD组
                group_results = self.solve_groups_independently(current_paths, prev_epsilon)
                if group_results:
                    # 如果独立求解成功，合并结果
                    effective_paths_by_group = {}
                    all_effective_paths = []
                    current_epsilon = {}
                    
                    for group_name, result in group_results.items():
                        effective_paths_by_group[group_name] = result['effective_paths']
                        all_effective_paths.extend(result['effective_paths'])
                        current_epsilon[group_name] = result['epsilon']
                    
                    # 记录本次迭代结果
                    iteration_data = {
                        'iteration': iteration,
                        'restricted_paths': {g: paths.copy() for g, paths in current_paths.items()},
                        'effective_paths': effective_paths_by_group,
                        'epsilon': current_epsilon,
                        'path_costs': {}, # 独立求解时无法获取统一的路径成本
                        'money_costs': {},
                        'flows': {}
                    }
                    results.append(iteration_data)
                    
                    # 更新路径和epsilon
                    for group_name, result in group_results.items():
                        new_paths = [p for p in result['effective_paths'] if p not in current_paths[group_name]]
                        if new_paths:
                            current_paths[group_name].extend(new_paths)
                            current_paths[group_name].sort()
                        prev_epsilon[group_name] = result['epsilon']
                        
                        # 输出该组的迭代信息
                        self.console.print(f"\n{group_name} 迭代 {iteration} (独立求解):")
                        self.console.print(f"当前限制路径: {current_paths[group_name]}")
                        self.console.print(f"新增有效路径: {new_paths}")
                        self.console.print(f"当前epsilon: {result['epsilon']:.3f}")
                    
                    iteration += 1
                    relax_factor = 0.99  # 重置放松因子
                    continue
                
                # 如果独立求解也失败，回退到使用已有结果
                if results:
                    self.console.print("[yellow]无法继续求解，使用最后一次有效结果[/yellow]")
                    break
                else:
                    # 如果没有任何结果，进行无约束求解
                    self.console.print("[yellow]尝试无约束求解...[/yellow]")
                    self.model = ConcreteModel()
                    self.initialize_sets()
                    self.initialize_parameters()
                    self.initialize_variables()
                    self.add_constraints()
                    self.set_objective()
                    solve_status = self.solve(tee=False)
                    if solve_status.solver.status != SolverStatus.ok:
                        self.console.print("[red]无约束求解也失败，退出[/red]")
                        break

            # 分析结果，修改为按OD组返回有效路径
            effective_paths_by_group, all_effective_paths = self.analyze_path_costs_by_group()

            # 获取当前epsilon - 已经是按OD组的字典
            if len(self.config.od_groups) > 1:
                current_epsilon = {g: self.model.epsilons[g].value for g in self.config.od_groups.keys()}
            else:
                group_name = list(self.config.od_groups.keys())[0]
                current_epsilon = {group_name: self.model.epsilon.value}

            # 记录本次迭代结果
            iteration_data = {
                'iteration': iteration,
                'restricted_paths': {g: paths.copy() for g, paths in current_paths.items()},
                'effective_paths': effective_paths_by_group,
                'epsilon': current_epsilon,
                'path_costs': {i: self.model.path_cost[i].value for i in self.model.od_pairs},
                'money_costs': self.calculate_money_cost(),
                'flows': {i: self.model.flow[i].value for i in self.model.od_pairs}
            }
            results.append(iteration_data)

            # 生成结果图表
            self.plot_cost_analysis(all_effective_paths, iteration_data)

            # 检查每个OD组的新路径并更新
            has_new_paths = False
            all_paths_complete = True
            
            for group_name, group_pairs in self.config.od_groups.items():
                # 计算该组的新路径
                group_effective_paths = effective_paths_by_group[group_name]
                new_paths = [p for p in group_effective_paths if p not in current_paths[group_name]]
                
                # 更新状态变量
                has_new_paths = has_new_paths or bool(new_paths)
                
                # 检查该OD组是否所有路径都已考虑
                all_group_paths = [p for p in group_pairs]
                is_group_complete = set(current_paths[group_name] + new_paths) == set(all_group_paths)
                all_paths_complete = all_paths_complete and is_group_complete
                
                # 更新该组的路径集合
                if new_paths:
                    current_paths[group_name].extend(new_paths)
                    current_paths[group_name].sort()
                
                # 更新该组的epsilon
                prev_epsilon[group_name] = current_epsilon[group_name]
                
                # 输出该组的迭代信息
                self.console.print(f"\n{group_name} 迭代 {iteration}:")
                self.console.print(f"当前限制路径: {current_paths[group_name]}")
                self.console.print(f"新增有效路径: {new_paths}")
                self.console.print(f"当前epsilon: {current_epsilon[group_name]:.3f}")

            # 终止条件：所有OD组都没有新路径，或所有OD组的所有路径都已考虑
            if not has_new_paths or all_paths_complete:
                self.console.print(f"[green]迭代完成，共 {iteration} 次迭代[/green]")
                break

            iteration += 1
            # 重置放松因子
            relax_factor = 0.99

        # 如果达到最大迭代次数
        if iteration > max_iterations:
            self.console.print(f"[yellow]已达最大迭代次数 {max_iterations}，停止迭代[/yellow]")

        # 显示最终结果汇总
        self.display_iteration_results_by_group(results)
        return results

    def add_path_constraints_by_group(self, restricted_paths_by_group=None, prev_epsilon=None, relax_factor=0.99):
        """添加按OD组的路径约束
        Args:
            restricted_paths_by_group: 字典，键为OD组名称，值为该组的受限路径列表
            prev_epsilon: 字典，键为OD组名称，值为该组的上一次迭代epsilon值
            relax_factor: 流量约束的松弛因子，值越小约束越松
        """
        m = self.model

        # 添加路径流量约束
        m.path_constraint = ConstraintList()
        for group_name, group_pairs in self.config.od_groups.items():
            if restricted_paths_by_group and group_name in restricted_paths_by_group:
                group_restricted_paths = [p for p in restricted_paths_by_group[group_name] if p in group_pairs]
                # 只在有受限路径时添加约束
                if group_restricted_paths:
                    m.path_constraint.add(
                        sum(m.flow[i] for i in group_restricted_paths) <=
                        self.config.od_demands[group_name] * relax_factor
                    )

        # 如果有前一次迭代的epsilon，添加epsilon下界约束，但允许一定的松弛
        if prev_epsilon is not None:
            if len(self.config.od_groups) > 1:
                for group_name in self.config.od_groups.keys():
                    if prev_epsilon[group_name] is not None:
                        # 松弛epsilon的下界约束，允许一定的波动
                        m.path_constraint.add(m.epsilons[group_name] >= prev_epsilon[group_name] * 0.9)
            else:
                group_name = list(self.config.od_groups.keys())[0]
                if prev_epsilon[group_name] is not None:
                    m.path_constraint.add(m.epsilon >= prev_epsilon[group_name] * 0.9)

    def solve_groups_independently(self, current_paths, prev_epsilon):
        """单独为每个OD组求解
        
        Args:
            current_paths: 当前限制路径
            prev_epsilon: 上一次迭代的epsilon值
            
        Returns:
            Dict: 各个OD组的求解结果
        """
        results = {}
        
        for group_name, group_pairs in self.config.od_groups.items():
            self.console.print(f"[blue]尝试单独求解 {group_name}...[/blue]")
            
            # 创建单独的模型
            group_model = ConcreteModel()
            
            # 创建组变量
            group_model.od_pairs = Set(initialize=group_pairs)
            group_model.paths = RangeSet(1, self.config.num_links)
            
            # 创建路径流量变量
            group_model.flow = Var(group_model.od_pairs, domain=NonNegativeReals)
            group_model.travel_time = Var(group_model.paths, domain=NonNegativeReals)
            group_model.path_cost = Var(group_model.od_pairs, domain=NonNegativeReals)
            group_model.residual = Var(group_model.od_pairs, domain=NonNegativeReals)
            group_model.epsilon = Var(domain=NonNegativeReals)
            group_model.perception = Var(domain=NonNegativeReals)
            
            # 添加约束
            # 需求约束
            group_model.demand_constraint = Constraint(
                expr=sum(group_model.flow[i] for i in group_pairs) == self.config.od_demands[group_name]
            )
            
            # 限制路径约束
            if group_name in current_paths and current_paths[group_name]:
                group_model.path_constraint = Constraint(
                    expr=sum(group_model.flow[i] for i in current_paths[group_name]) <= 
                    self.config.od_demands[group_name] * 0.9  # 使用较低的约束值
                )
            
            # 误差约束
            def epsilon_rule(m, i):
                return m.epsilon >= m.residual[i]
            group_model.epsilon_constraints = Constraint(group_model.od_pairs, rule=epsilon_rule)
            
            # 平衡约束
            def balance_rule(m, i):
                return m.flow[i] + m.path_cost[i] + m.residual[i] - m.perception >= 0.0
            group_model.balance_constraints = Constraint(group_model.od_pairs, rule=balance_rule)
            
            # 路径成本非负约束
            def cost_nonnegativity_rule(m, i):
                return m.path_cost[i] + m.residual[i] - m.perception >= 0
            group_model.cost_nonnegativity = Constraint(group_model.od_pairs, rule=cost_nonnegativity_rule)
            
            # 总成本约束
            group_model.total_cost_constraint = Constraint(
                expr=sum(group_model.flow[i] * (group_model.path_cost[i] + 
                    group_model.residual[i] - group_model.perception) for i in group_pairs) == 0
            )
            
            # 旅行时间约束 cost = t0(1 + 0.15(f/c)^4)
            def travel_time_rule(m, j):
                # 只考虑当前组中路径使用的链路
                return m.travel_time[j] == self.config.free_flow_time[j] * (
                    1 + 0.15 * (sum(self.config.path_link_matrix.get((i, j), 0) * m.flow[i] 
                                for i in group_pairs) / self.config.link_capacity[j]) ** 4
                )
            group_model.travel_time_constraints = Constraint(group_model.paths, rule=travel_time_rule)
            
            # 路径成本计算约束
            def path_cost_rule(m, i):
                return m.path_cost[i] == sum(
                    self.config.path_link_matrix.get((i, j), 0) * m.travel_time[j] for j in group_model.paths
                )
            group_model.path_cost_constraints = Constraint(group_model.od_pairs, rule=path_cost_rule)
            
            # 目标函数
            group_model.objective = Objective(expr=group_model.epsilon, sense=minimize)
            
            # 求解
            solver = SolverFactory('ipopt')
            solver.options['tol'] = 1e-6
            solver.options['max_iter'] = 3000
            try:
                solve_result = solver.solve(group_model, tee=False)
                
                if solve_result.solver.status == SolverStatus.ok:
                    self.console.print(f"[green]{group_name} 单独求解成功[/green]")
                    
                    # 分析有效路径
                    epsilon_value = group_model.epsilon.value
                    min_cost = min(group_model.path_cost[i].value for i in group_pairs)
                    upper_bound = min_cost + epsilon_value
                    effective_paths = []
                    
                    for i in group_pairs:
                        cost = group_model.path_cost[i].value
                        is_effective = abs(cost - min_cost) <= 1e-2 or cost <= upper_bound + 1e-2
                        if is_effective:
                            effective_paths.append(i)
                    
                    results[group_name] = {
                        'effective_paths': effective_paths,
                        'epsilon': epsilon_value,
                        'flows': {i: group_model.flow[i].value for i in group_pairs}
                    }
                    
                    self.console.print(f"{group_name} 有效路径: {sorted(effective_paths)}")
                    self.console.print(f"{group_name} epsilon: {epsilon_value:.3f}")
                else:
                    self.console.print(f"[red]{group_name} 单独求解失败[/red]")
            except Exception as e:
                self.console.print(f"[red]{group_name} 单独求解出错: {str(e)}[/red]")
                
        return results if results else None

    def analyze_path_costs_by_group(self):
        """分析路径成本并按OD组返回有效路径"""
        m = self.model
        effective_paths_by_group = {group_name: [] for group_name in self.config.od_groups.keys()}
        all_effective_paths = []  # 用于向后兼容

        # 分析每个OD组
        for group_name, group_pairs in self.config.od_groups.items():
            if len(self.config.od_groups) > 1:
                epsilon = m.epsilons[group_name].value
            else:
                epsilon = m.epsilon.value
            min_cost = min(m.path_cost[i].value for i in group_pairs)
            upper_bound = min_cost + epsilon

            # 创建结果表格
            table = Table(title=f"{group_name} 路径成本分析")
            table.add_column("路径", style="cyan")
            table.add_column("成本", style="magenta")
            table.add_column("流量", style="green")
            table.add_column("有效性", style="yellow")
            table.add_column("相对差异", style="blue")  # 添加与最小成本的差异

            # 分析每条路径，修改判断有效路径的逻辑
            for i in group_pairs:
                cost = m.path_cost[i].value
                flow = m.flow[i].value
                # 与plot_cost_analysis保持一致的判断标准
                is_effective = abs(cost - min_cost) <= 1e-2 or cost <= upper_bound + 1e-2
                relative_diff = cost - min_cost  # 相对于最小成本的差异

                table.add_row(
                    str(i),
                    f"{cost:.3f}",
                    f"{flow:.3f}",
                    "✓" if is_effective else "✗",
                    f"{relative_diff:.3f}"
                )

                if is_effective:
                    effective_paths_by_group[group_name].append(i)
                    all_effective_paths.append(i)  # 添加到总的有效路径列表

            self.console.print(table)
            
            # 添加组路径总结
            self.console.print(f"{group_name} 有效路径: {sorted(effective_paths_by_group[group_name])}\n")

        # 为了保持向后兼容，返回全局的有效路径列表
        return effective_paths_by_group, all_effective_paths

    def display_iteration_results_by_group(self, results):
        """按OD组显示迭代结果汇总"""
        # 为每个OD组创建一个表格
        for group_name in self.config.od_groups.keys():
            table = Table(title=f"{group_name} 迭代分析结果")
            table.add_column("迭代", style="cyan")
            table.add_column("限制路径", style="magenta")
            table.add_column("有效路径", style="green")
            table.add_column("新增路径", style="blue")
            table.add_column("Epsilon", style="yellow")

            for i, result in enumerate(results):
                # 确保数据存在
                if group_name not in result['effective_paths']:
                    continue
                
                # 计算新增路径
                prev_paths = set()
                if i > 0 and group_name in results[i - 1]['effective_paths']:
                    prev_paths = set(results[i - 1]['effective_paths'][group_name])
                
                new_paths = sorted(list(set(result['effective_paths'][group_name]) - prev_paths))
                
                # 获取当前epsilon
                epsilon_value = result['epsilon'].get(group_name, 0.0)
                
                table.add_row(
                    str(result['iteration']),
                    str(sorted(result['restricted_paths'].get(group_name, []))),
                    str(sorted(result['effective_paths'][group_name])),
                    str(new_paths),
                    f"{epsilon_value:.3f}"
                )

            self.console.print(table)
            self.console.print("\n")


def main():
    base_config = TrafficNetworkConfig.create_basic_network()
    base_solver = BRUESolver(base_config)
    base_solver.run_with_iterations()
    base_solver.plot_initial_costs()

    # # 测试简单网络
    # simple_config = TrafficNetworkConfig.create_single_od_network()
    # simple_solver = BRUESolver(simple_config)
    # simple_solver.run_with_iterations()
    # simple_solver.plot_initial_costs()

    # # 测试路径网络
    # path_config = TrafficNetworkConfig.create_multi_od_network()
    # path_solver = BRUESolver(path_config)
    # path_solver.run_with_iterations()
    # path_solver.plot_initial_costs()

    # 测试两起终点对网络
    # two_od_config = TrafficNetworkConfig.create_two_od_network()
    # two_od_solver = BRUESolver(two_od_config)
    # two_od_solver.run_with_iterations()
    # two_od_solver.plot_initial_costs()

    # 测试反例网络
    # anti_example_config = TrafficNetworkConfig.create_anti_example_network()
    # anti_example_solver = BRUESolver(anti_example_config)
    # anti_example_solver.run_with_iterations()
    # anti_example_solver.plot_initial_costs()


if __name__ == "__main__":
    main()
