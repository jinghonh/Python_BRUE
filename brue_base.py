from pyomo.environ import *
from pyomo.opt import SolverFactory
from abc import ABC, abstractmethod
from rich.console import Console
from rich.table import Table

class BRUEBase(ABC):
    def __init__(self):
        self.model = ConcreteModel()
        self.console = Console()
        
    @abstractmethod
    def initialize_sets(self):
        """初始化集合"""
        pass

    @abstractmethod
    def initialize_parameters(self):
        """初始化参数"""
        pass

    @abstractmethod
    def initialize_variables(self):
        """初始化变量"""
        pass

    @abstractmethod
    def add_constraints(self):
        """添加约束"""
        pass

    @abstractmethod
    def set_objective(self):
        """设置目标函数"""
        pass

    @abstractmethod
    def analyze_path_costs(self):
        """分析路径成本并找出有效路径"""
        pass

    @abstractmethod
    def calculate_money_cost(self):
        """计算每条路径的金钱成本"""
        pass

    @abstractmethod
    def plot_cost_analysis(self, is_effective_path):
        """绘制路径成本与金钱成本的散点图"""
        pass

    def solve(self, solver_name='ipopt', tee=True):
        """求解模型"""
        solver = SolverFactory(solver_name)
        self.results = solver.solve(self.model, tee=tee)
        return self.results

    def display_results(self):
        """显示结果"""
        pass

    def run(self):
        """运行完整求解过程"""
        self.initialize_sets()
        self.initialize_parameters()
        self.initialize_variables()
        self.add_constraints()
        self.set_objective()
        self.solve()
        self.display_results()
        self.analyze_path_costs()  # 添加新的分析步骤
