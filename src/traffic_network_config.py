from dataclasses import dataclass
from typing import Dict, Set, List
import numpy as np


@dataclass
class TrafficNetworkConfig:
    """交通网络配置类，用于存储和管理交通网络的所有参数"""
    # 基础参数
    num_links: int  # 路径数量
    num_paths: int  # 起终点对数量
    total_demand: float  # 总交通需求量
    od_groups: Dict[str, List[int]]  # 起终点组
    od_demands: Dict[str, float]  # 每个起终点组的需求量

    # 成本参数
    free_flow_time: Dict[int, float]  # 自由流时间（无拥堵时的行程时间）
    link_money_cost: Dict[int, float]  # 道路费用（如过路费）
    link_capacity: Dict[int, float]  # 道路通行能力

    # 路径-链接关系
    path_link_matrix: Dict[tuple, float]  # 路径与道路段的关联矩阵

    def get_path_link_array(self) -> "np.ndarray":
        """以矩阵形式返回路径-道路段关系"""
        import numpy as np

        num_paths = self.num_paths
        num_links = self.num_links
        matrix = np.zeros((num_paths, num_links))

        for (i, j), v in self.path_link_matrix.items():
            if i <= num_paths and j <= num_links:
                matrix[i - 1, j - 1] = v

        return matrix

    @classmethod
    def create_basic_network(cls):
        """创建基础示例网络配置"""
        return cls(
            num_links=8,
            num_paths=6,
            total_demand=10000,
            od_groups={'ALL': list(range(1, 7))},
            od_demands={'ALL': 10000},
            free_flow_time={1: 18, 2: 22.5, 3: 12, 4: 24, 5: 2.4, 6: 6, 7: 24, 8: 12},
            link_money_cost={1: 20, 2: 15, 3: 1, 4: 0, 5: 0, 6: 0, 7: 0, 8: 1},
            link_capacity={1: 3600, 2: 3600, 3: 1800, 4: 1800, 5: 1800, 6: 1800, 7: 1800, 8: 1800},
            path_link_matrix={
                (1, 1): 1, (1, 2): 0, (1, 3): 0, (1, 4): 0, (1, 5): 0, (1, 6): 0, (1, 7): 0, (1, 8): 0,
                (2, 1): 0, (2, 2): 1, (2, 3): 0, (2, 4): 0, (2, 5): 0, (2, 6): 0, (2, 7): 0, (2, 8): 0,
                (3, 1): 0, (3, 2): 0, (3, 3): 1, (3, 4): 0, (3, 5): 0, (3, 6): 0, (3, 7): 1, (3, 8): 0,
                (4, 1): 0, (4, 2): 0, (4, 3): 0, (4, 4): 1, (4, 5): 0, (4, 6): 0, (4, 7): 0, (4, 8): 1,
                (5, 1): 0, (5, 2): 0, (5, 3): 1, (5, 4): 0, (5, 5): 1, (5, 6): 0, (5, 7): 0, (5, 8): 1,
                (6, 1): 0, (6, 2): 0, (6, 3): 0, (6, 4): 1, (6, 5): 0, (6, 6): 1, (6, 7): 1, (6, 8): 0
            }
        )

    @classmethod
    def create_multi_od_network(cls):
        """创建多起终点对网络配置"""
        return cls(
            num_links=19,
            num_paths=14,
            total_demand=6000,  # 总需求（OD1+OD2）
            od_groups={
                'OD1': list(range(1, 9)),
                'OD2': list(range(9, 15))
            },
            od_demands={
                'OD1': 3000,  # OD1组的需求
                'OD2': 3000  # OD2组的需求
            },
            free_flow_time={
                1: 10, 2: 10, 3: 10, 4: 40, 5: 10, 6: 10, 7: 10, 8: 20, 9: 10,
                10: 25, 11: 10, 12: 10, 13: 40, 14: 10, 15: 10, 16: 10, 17: 10, 18: 80, 19: 10
            },
            link_money_cost={
                1: 10, 2: 10, 3: 10, 4: 20, 5: 10, 6: 10, 7: 10, 8: 20, 9: 10,
                10: 25, 11: 10, 12: 10, 13: 20, 14: 10, 15: 10, 16: 10, 17: 10, 18: 30, 19: 10
            },
            link_capacity={i: 2500 for i in range(1, 20)},
            path_link_matrix={
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
        )

    @classmethod
    def create_single_od_network(cls):
        """创建单一起终点对网络配置"""
        return cls(
            num_links=19,
            num_paths=8,  # 只有OD1组的8个OD对
            total_demand=3000,
            od_groups={'OD1': list(range(1, 9))},  # 只保留OD1组
            od_demands={'OD1': 3000},  # OD1组的需求
            free_flow_time={
                1: 10, 2: 10, 3: 10, 4: 40, 5: 10, 6: 10, 7: 10, 8: 20, 9: 10,
                10: 25, 11: 10, 12: 10, 13: 40, 14: 10, 15: 10, 16: 10, 17: 10, 18: 80, 19: 10
            },
            link_money_cost={
                1: 10, 2: 10, 3: 10, 4: 20, 5: 10, 6: 10, 7: 10, 8: 20, 9: 10,
                10: 25, 11: 10, 12: 10, 13: 20, 14: 10, 15: 10, 16: 10, 17: 10, 18: 30, 19: 10
            },
            link_capacity={i: 2500 for i in range(1, 20)},
            path_link_matrix={
                (1, 1): 1, (1, 5): 1, (1, 7): 1, (1, 9): 1, (1, 11): 1,
                (2, 1): 1, (2, 5): 1, (2, 7): 1, (2, 10): 1, (2, 15): 1,
                (3, 1): 1, (3, 5): 1, (3, 8): 1, (3, 14): 1, (3, 15): 1,
                (4, 1): 1, (4, 6): 1, (4, 12): 1, (4, 14): 1, (4, 15): 1,
                (5, 2): 1, (5, 7): 1, (5, 9): 1, (5, 11): 1, (5, 17): 1,
                (6, 2): 1, (6, 7): 1, (6, 10): 1, (6, 15): 1, (6, 17): 1,
                (7, 2): 1, (7, 8): 1, (7, 14): 1, (7, 15): 1, (7, 17): 1,
                (8, 2): 1, (8, 11): 1, (8, 18): 1
            }
        )
    
    @classmethod
    def create_two_od_network(cls):
        """创建两起终点对网络配置"""
        return cls(
            num_links=4,
            num_paths=6,
            total_demand=3000,
            od_groups={
                'OD1': list(range(5,7)),
                'OD2': list(range(1,5))
                },
            od_demands={
                'OD1': 1, 
                'OD2': 1
                },
            free_flow_time={1: 0, 2: 10, 3: 5, 4: 0},
            link_money_cost={1: 10, 2: 20, 3: 15, 4: 10},
            link_capacity={1: 0.5, 2: 1, 3: 0.5, 4: 0.5},
            path_link_matrix={
                (1, 1): 1, (1, 3): 1,
                (2, 1): 1, (2, 4): 1,
                (3, 2): 1, (3, 3): 1,
                (4, 2): 1, (4, 4): 1,
                (5, 1): 1, 
                (6, 2): 1,
            }
        )
    
    #创建一个反例网络配置
    @classmethod
    def create_anti_example_network(cls):
        """创建一个反例网络配置"""
        return cls(
            num_links=12,
            num_paths=10,
            total_demand=10000,
            od_groups={
                'OD1': list(range(1, 6)),
                'OD2': list(range(6, 11))
            },
            od_demands={
                'OD1': 5000,
                'OD2': 5000
            },
            free_flow_time={1:4,2:5,3:3,4:2,5:3,6:1,7:1,8:5,9:3,10:4,11:1,12:5},#FFT = [4, 5, 3, 2, 3, 1, 1, 5, 3, 4, 1, 5]
            link_money_cost={1:2,2:1,3:3,4:4,5:3,6:5,7:5,8:1,9:3,10:2,11:5,12:1},#M = [2, 1, 3, 4, 3, 5, 5, 1, 3, 2, 5, 1] 
            link_capacity={1:2000,2:2000,3:2000,4:2000,5:2000,6:2000,7:2000,8:2000,9:2000,10:2000,11:2000,12:2000},
            path_link_matrix={
                (1, 1): 1, (1, 7): 1,# OD1-Path1
                (2, 1): 1, (2, 8): 1, (2, 2): 1,# OD1-Path2
                (3, 2): 1, (3, 9): 1, (3, 12): 1,# OD1-Path3
                (4, 3): 1, (4, 11): 1,# OD1-Path4
                (5, 3): 1, (5, 10): 1, (5, 12): 1,# OD1-Path5
                (6, 6): 1, (6, 11): 1,# OD2-Path1
                (7, 6): 1, (7, 10): 1, (7, 12): 1,# OD2-Path2
                (8, 5): 1, (8, 9): 1, (8, 12): 1,# OD2-Path3
                (9, 4): 1, (9, 7): 1,# OD2-Path4
                (10, 4): 1, (10, 8): 1, (10, 12): 1,# OD2-Path5
            }
        )

    @classmethod
    def create_large_single_od_network(cls):
        """创建单一OD对多路径网络配置（至少15条路径）"""
        num_paths = 15
        num_links = 30
        
        # 初始化空的path_link_matrix字典
        path_link_matrix = {}
        
        # 创建每条路径的link连接
        for i in range(1, num_paths + 1):
            # 为每条路径分配2-4条link
            num_links_for_path = min(2 + (i % 3), 4)  # 使用2到4条link
            
            # 选择连续的link作为这条路径的组成部分
            start_link = (i % (num_links - num_links_for_path)) + 1
            
            # 为每条路径创建link连接
            for j in range(start_link, start_link + num_links_for_path):
                path_link_matrix[(i, j)] = 1
        
        # 创建link参数
        free_flow_time = {}
        link_money_cost = {}
        link_capacity = {}
        
        for j in range(1, num_links + 1):
            # 随机分配参数，但保持一定的模式
            free_flow_time[j] = 10 + (j % 10) * 2  # 10-28的自由流时间
            link_money_cost[j] = 5 + j % 15  # 5-19的金钱成本
            link_capacity[j] = 2000 + (j % 5) * 200  # 2000-2800的通行能力
        
        return cls(
            num_links=num_links,
            num_paths=num_paths,
            total_demand=5000,
            od_groups={'OD1': list(range(1, num_paths + 1))},  # 单OD对，15条路径
            od_demands={'OD1': 5000},  # OD1组的需求
            free_flow_time=free_flow_time,
            link_money_cost=link_money_cost,
            link_capacity=link_capacity,
            path_link_matrix=path_link_matrix
        )

