from pyomo.environ import *
from pyomo.opt import SolverFactory
import matplotlib.pyplot as plt
import numpy as np
from rich.console import Console
from rich.table import Table

from .traffic_network_config import TrafficNetworkConfig


class UESolver:
    """
    用户均衡(User Equilibrium)求解器
    
    实现了基于Beckmann公式的经典交通分配模型，用于求解网络中的流量均衡。
    """
    
    def __init__(self, config: TrafficNetworkConfig):
        """
        初始化求解器
        
        Args:
            config: 交通网络配置
        """
        self.config = config
        self.model = ConcreteModel()
        self.console = Console()
        self.results = None
        
    def initialize_sets(self):
        """初始化模型集合"""
        self.model.od_pairs = RangeSet(1, self.config.num_paths)  # 路径
        self.model.links = RangeSet(1, self.config.num_links)  # 链接（道路段）
        
        # 初始化OD组
        for group_name, group_pairs in self.config.od_groups.items():
            setattr(self.model, f'OD_{group_name}', Set(initialize=group_pairs))
    
    def initialize_parameters(self):
        """初始化模型参数"""
        m = self.model
        
        # 初始化链接参数
        m.free_flow_time = Param(
            m.links,
            initialize=self.config.free_flow_time
        )
        m.link_money_cost = Param(
            m.links,
            initialize=self.config.link_money_cost
        )
        m.link_capacity = Param(
            m.links,
            initialize=self.config.link_capacity
        )
        
        # 初始化路径-链接关系矩阵
        m.path_link_matrix = Param(
            m.od_pairs,
            m.links,
            initialize=self.config.path_link_matrix,
            default=0
        )
    
    def initialize_variables(self):
        """初始化模型变量"""
        m = self.model
        
        # 路径流量变量
        m.flow = Var(m.od_pairs, domain=NonNegativeReals)
        
        # 链接流量变量
        m.link_flow = Var(m.links, domain=NonNegativeReals)
        
        # 链接时间成本变量
        m.link_time = Var(m.links, domain=NonNegativeReals)
        
        # 路径时间成本变量
        m.path_time = Var(m.od_pairs, domain=NonNegativeReals)
    
    def add_constraints(self):
        """添加模型约束"""
        m = self.model
        
        # OD需求约束
        m.demand_constraint = ConstraintList()
        for group_name, group_pairs in self.config.od_groups.items():
            m.demand_constraint.add(
                sum(m.flow[i] for i in group_pairs) == self.config.od_demands[group_name]
            )
        
        # 链接流量计算约束
        def link_flow_rule(m, j):
            return m.link_flow[j] == sum(m.path_link_matrix[i, j] * m.flow[i] 
                                        for i in m.od_pairs)
        
        m.link_flow_constraints = Constraint(m.links, rule=link_flow_rule)
        
        # 链接时间计算约束（BPR函数）
        def link_time_rule(m, j):
            return m.link_time[j] == m.free_flow_time[j] * (
                    1 + 0.15 * (m.link_flow[j] / m.link_capacity[j]) ** 4
            )
        
        m.link_time_constraints = Constraint(m.links, rule=link_time_rule)
        
        # 路径时间计算约束
        def path_time_rule(m, i):
            return m.path_time[i] == sum(
                m.path_link_matrix[i, j] * m.link_time[j] for j in m.links
            )
        
        m.path_time_constraints = Constraint(m.od_pairs, rule=path_time_rule)
    
    def set_objective(self):
        """设置目标函数（Beckmann公式）"""
        m = self.model
        
        # 使用链接积分作为目标函数
        def link_integral_rule(m, j):
            # BPR函数的积分
            return m.free_flow_time[j] * (
                    m.link_flow[j] + 0.15 * (m.link_flow[j]**5) / (5 * m.link_capacity[j]**4)
            )
        
        m.objective = Objective(
            expr=sum(link_integral_rule(m, j) for j in m.links),
            sense=minimize
        )
    
    def solve(self, solver_name='ipopt', tee=True):
        """求解模型"""
        solver = SolverFactory(solver_name)
        solver.options['max_iter'] = 5000
        self.results = solver.solve(self.model, tee=tee)
        return self.results
    
    def calculate_link_flows(self):
        """计算最终链接流量"""
        link_flows = {}
        for j in self.model.links:
            link_flows[j] = self.model.link_flow[j].value
        return link_flows
    
    def calculate_link_times(self):
        """计算最终链接时间"""
        link_times = {}
        for j in self.model.links:
            link_times[j] = self.model.link_time[j].value
        return link_times
    
    def calculate_money_costs(self):
        """计算每条路径的金钱成本"""
        m = self.model
        money_costs = {}
        for i in m.od_pairs:
            money_costs[i] = sum(m.path_link_matrix[i, j] * m.link_money_cost[j]
                               for j in m.links)
        return money_costs
    
    def display_results(self):
        """显示计算结果"""
        table = Table(title="用户均衡解算结果")
        
        columns = ["路径", "流量", "时间成本", "金钱成本", "OD组"]
        for col in columns:
            table.add_column(col, justify="right")
        
        money_costs = self.calculate_money_costs()
        
        for group_name, group_pairs in self.config.od_groups.items():
            for i in group_pairs:
                table.add_row(
                    f"{i}",
                    f"{self.model.flow[i].value:.2f}",
                    f"{self.model.path_time[i].value:.2f}",
                    f"{money_costs[i]:.2f}",
                    f"{group_name}"
                )
        
        self.console.print(table)
        
        # 显示每个OD组的总流量和平均成本
        table_summary = Table(title="OD组汇总")
        table_summary.add_column("OD组", justify="left")
        table_summary.add_column("总流量", justify="right")
        table_summary.add_column("平均时间成本", justify="right")
        table_summary.add_column("平均金钱成本", justify="right")
        
        for group_name, group_pairs in self.config.od_groups.items():
            total_flow = sum(self.model.flow[i].value for i in group_pairs)
            avg_time = sum(self.model.flow[i].value * self.model.path_time[i].value 
                         for i in group_pairs) / total_flow if total_flow > 0 else 0
            avg_money = sum(self.model.flow[i].value * money_costs[i] 
                          for i in group_pairs) / total_flow if total_flow > 0 else 0
            
            table_summary.add_row(
                group_name,
                f"{total_flow:.2f}",
                f"{avg_time:.2f}",
                f"{avg_money:.2f}"
            )
        
        self.console.print(table_summary)
    
    def plot_link_flows(self):
        """绘制链接流量图"""
        link_flows = self.calculate_link_flows()
        link_ids = list(link_flows.keys())
        flow_values = [link_flows[j] for j in link_ids]
        
        plt.figure(figsize=(10, 6))
        plt.bar(link_ids, flow_values)
        plt.xlabel('链接ID')
        plt.ylabel('流量')
        plt.title('链接流量分布')
        plt.grid(True, linestyle='--', alpha=0.6)
        plt.tight_layout()
        plt.savefig('output/link_flows.png')
        plt.close()
        
        return 'output/link_flows.png'
    
    def plot_path_distribution(self):
        """绘制路径分布图"""
        plt.figure(figsize=(12, 8))
        
        for group_name, group_pairs in self.config.od_groups.items():
            path_ids = list(group_pairs)
            flow_values = [self.model.flow[i].value for i in path_ids]
            
            plt.subplot(len(self.config.od_groups), 1, 
                      list(self.config.od_groups.keys()).index(group_name) + 1)
            plt.bar(path_ids, flow_values, label=group_name)
            plt.xlabel('路径ID')
            plt.ylabel('流量')
            plt.title(f'{group_name}组路径流量分布')
            plt.grid(True, linestyle='--', alpha=0.6)
            plt.legend()
        
        plt.tight_layout()
        plt.savefig('output/path_distribution.png')
        plt.close()
        
        return 'output/path_distribution.png'
    
    def plot_cost_flow_relationship(self):
        """绘制成本-流量关系图"""
        money_costs = self.calculate_money_costs()
        
        plt.figure(figsize=(10, 8))
        
        for group_name, group_pairs in self.config.od_groups.items():
            path_times = [self.model.path_time[i].value for i in group_pairs]
            path_money = [money_costs[i] for i in group_pairs]
            path_flows = [self.model.flow[i].value for i in group_pairs]
            
            # 散点大小与流量成比例
            sizes = [50 * flow + 10 for flow in path_flows]
            
            plt.scatter(path_times, path_money, s=sizes, alpha=0.7, label=group_name)
            
            # 标记路径ID
            for i, txt in enumerate(group_pairs):
                plt.annotate(txt, (path_times[i], path_money[i]),
                           xytext=(5, 5), textcoords='offset points')
        
        plt.xlabel('时间成本')
        plt.ylabel('金钱成本')
        plt.title('路径成本与流量分布')
        plt.grid(True, linestyle='--', alpha=0.6)
        plt.legend()
        plt.tight_layout()
        plt.savefig('output/cost_flow_relationship.png')
        plt.close()
        
        return 'output/cost_flow_relationship.png'
    
    def run(self):
        """运行完整求解过程"""
        self.initialize_sets()
        self.initialize_parameters()
        self.initialize_variables()
        self.add_constraints()
        self.set_objective()
        self.solve()
        self.display_results()
        
        # 生成可视化图表
        link_flow_plot = self.plot_link_flows()
        path_dist_plot = self.plot_path_distribution()
        cost_flow_plot = self.plot_cost_flow_relationship()
        
        return {
            "link_flow_plot": link_flow_plot,
            "path_distribution_plot": path_dist_plot,
            "cost_flow_plot": cost_flow_plot
        }

def main():
    """主函数，用于测试UE求解器"""
    # 创建网络配置
    config = TrafficNetworkConfig.create_multi_od_network()
    
    # 初始化并运行求解器
    solver = UESolver(config)
    plots = solver.run()
    
    print(f"链接流量图保存至: {plots['link_flow_plot']}")
    print(f"路径分布图保存至: {plots['path_distribution_plot']}")
    print(f"成本-流量关系图保存至: {plots['cost_flow_plot']}")

if __name__ == "__main__":
    main() 