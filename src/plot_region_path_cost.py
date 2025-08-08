"""demo.py – Python 版的 BRUE 路径成本区域绘制示例
====================================================
该脚本根据给出的 Mathematica 代码 (见 README) 重写而成，
完成如下功能：
1. 定义链路/路径成本函数；
2. 构造三类可行域 (reg, reg2, reg3)；
3. 计算各路径成本的最值及中点；
4. 采样并绘制区域及给定的方案点；
5. 导出多个图表展示不同区域组合。

与 Mathematica 的 RegionPlot 不同，Python 这里采用网格采样 +
布尔掩码的方式近似绘制区域。若需更高精度，可适当缩小步长。
"""

from __future__ import annotations

import math
import json
import os
import time
from pathlib import Path
from typing import Callable, Tuple, Dict, Any, List

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
# 绘图相关
from matplotlib.patches import Patch
import matplotlib.patheffects as patheffects  # 用于给文字添加描边效果
from scipy.interpolate import interp1d

# 从统一配置读取区域样式
from plot_styles import REGION_STYLES



def link1(f1: np.ndarray | float) -> np.ndarray | float:
    return 18.0 * (1 + 0.15 * (f1 / 3600.0) ** 4)


def link2(f2: np.ndarray | float) -> np.ndarray | float:
    return 22.5 * (1 + 0.15 * (f2 / 3600.0) ** 4)


def _shared_denom(f1: np.ndarray | float, f2: np.ndarray | float) -> np.ndarray | float:
    """Denominator  (10000 - f1 - f2) / 1800  多处复用。"""
    return (10000.0 - f1 - f2) / 1800.0


def link3(f1: np.ndarray | float, f2: np.ndarray | float) -> np.ndarray | float:
    return 12.0 * (1 + 0.15 * _shared_denom(f1, f2) ** 4)


def link5(f1: np.ndarray | float, f2: np.ndarray | float) -> np.ndarray | float:
    return 2.4 * (1 + 0.15 * _shared_denom(f1, f2) ** 4)


def link8(f1: np.ndarray | float, f2: np.ndarray | float) -> np.ndarray | float:
    return 12.0 * (1 + 0.15 * _shared_denom(f1, f2) ** 4)


def _rg_cost(raw_time: np.ndarray | float) -> np.ndarray | float:
    """为原始旅行时间添加感知成本。"""
    return raw_time + RHO * (1 - np.exp(-SIGMA * raw_time))


def path1(f1: np.ndarray | float, f2: np.ndarray | float) -> np.ndarray | float:
    return _rg_cost(link1(f1))


def path2(f1: np.ndarray | float, f2: np.ndarray | float) -> np.ndarray | float:
    return _rg_cost(link2(f2))


def path5(f1: np.ndarray | float, f2: np.ndarray | float) -> np.ndarray | float:
    raw = link3(f1, f2) + link5(f1, f2) + link8(f1, f2)
    return _rg_cost(raw)


def max_path(f1: np.ndarray | float, f2: np.ndarray | float) -> np.ndarray | float:
    return np.maximum.reduce([path1(f1, f2), path2(f1, f2), path5(f1, f2)])


def min_path(f1: np.ndarray | float, f2: np.ndarray | float) -> np.ndarray | float:
    return np.minimum.reduce([path1(f1, f2), path2(f1, f2), path5(f1, f2)])

# ============================ 约束工具函数 ===========================

def _order_constraints(f1: np.ndarray, f2: np.ndarray) -> np.ndarray:
    """path1 <= path2 - EPS, path1 <= path5 - EPS, path2 <= path5 - EPS"""
    c1 = path1(f1, f2) - path2(f1, f2) + EPS <= 0
    c2 = path1(f1, f2) - path5(f1, f2) + EPS <= 0
    c3 = path2(f1, f2) - path5(f1, f2) + EPS <= 0
    return c1 & c2 & c3


def _pair_constraints(f1: np.ndarray, f2: np.ndarray) -> np.ndarray:
    diff12 = np.abs(path1(f1, f2) - path2(f1, f2)) <= E
    diff15 = np.abs(path1(f1, f2) - path5(f1, f2)) <= E
    diff25 = np.abs(path2(f1, f2) - path5(f1, f2)) <= E
    return diff12 & diff15 & diff25


def _box_constraints(f1: np.ndarray, f2: np.ndarray) -> np.ndarray:
    range_ok = (F1_MIN <= f1) & (f1 <= F1_MAX) & (F2_MIN <= f2) & (f2 <= F2_MAX)
    total_ok = (f1 + f2) <= TOTAL_MAX
    return range_ok & total_ok


def feasible_mask(f1: np.ndarray, f2: np.ndarray) -> np.ndarray:
    """整体可行域 reg 掩码。"""
    return _order_constraints(f1, f2) & _pair_constraints(f1, f2) & _box_constraints(
        f1, f2
    )

# ============================ 平衡线计算函数 ===========================

def calculate_equilibrium_line(left_boundary_x, upper_limit_x, upper_limit_y):
    """计算平衡线位置
    
    参数:
        left_boundary_x: 左边界X值
        upper_limit_x: 上限X值 (mid_path1, mid_path2, mid_path5)
        upper_limit_y: 上限Y值 (货币成本值)
        
    返回:
        eqm_limit_x: 平衡线X值
    """
    # 获取归一化的货币成本
    money_max = max(upper_limit_y)
    money_min = min(upper_limit_y)
    money_range = money_max - money_min
    
    if money_range > 0:
        # 计算权重: 较低的货币成本获得较高权重
        weights = 0.7 + 0.1 * ((upper_limit_y - money_min) / money_range)
    else:
        weights = 0.5 * np.ones_like(upper_limit_y)
    
    # 计算平衡线位置
    # 较低的货币成本（较大的权重）时，eqm_limit_x更接近左边界
    # 较高的货币成本（较小的权重）时，eqm_limit_x更接近上限
    eqm_limit_x = left_boundary_x + (upper_limit_x - left_boundary_x) * weights
    
    return eqm_limit_x

# ========================= JSON加载/保存功能 =========================

def get_or_compute_path_limits(left_boundary_x, upper_limit_x, upper_limit_y, e_value, results_dir='results', zeta=16):
    """
    获取TTB_max_delta和TTB_max值，从plot_path_costs.py生成的JSON文件中读取
    
    参数:
        left_boundary_x: 左边界X值
        upper_limit_x: 上限X值 (用于默认值，如果文件不存在)
        upper_limit_y: 上限Y值 (货币成本值)
        e_value: 误差上限值(用于日志显示)
        results_dir: 结果存储目录
        zeta: zeta值，对应plot_path_costs.py生成的文件名
        
    返回:
        Tuple[np.ndarray, np.ndarray]: mid_path(TTB_max_delta)和eqm_limit_x(TTB_max)数组
    """
    # 创建存储目录
    os.makedirs(results_dir, exist_ok=True)
    # 使用与plot_path_costs.py相同的文件名格式
    cache_file = os.path.join(results_dir, f"tmax_teqm_zeta{int(zeta)}.json")
    
    # 检查是否存在缓存文件
    if os.path.exists(cache_file):
        print(f"正在从 {cache_file} 加载 TTB_max_delta(mid_path) 和 TTB_max(eqm_limit_x) 数据")
        try:
            with open(cache_file, 'r') as f:
                data = json.load(f)
            
            # 从JSON加载数据 - plot_path_costs.py中的TTB_max_delta对应mid_path，TTB_max对应eqm_limit_x
            money_values = np.array(data['money_values'])
            TTB_max_delta_values = np.array(data['TTB_max_delta'])  # 对应mid_path
            TTB_max_values = np.array(data['TTB_max'])  # 对应eqm_limit_x
            
            # 如果货币值与缓存不同，需要进行插值
            if len(money_values) != len(upper_limit_y) or not np.allclose(money_values, upper_limit_y):
                # 为TTB_max_delta(mid_path)创建插值函数
                f_tmax = interp1d(
                    money_values, TTB_max_delta_values,
                    kind='linear', bounds_error=False, fill_value="extrapolate"
                )
                mid_path = f_tmax(upper_limit_y)
                
                # 为TTB_max(eqm_limit_x)创建插值函数
                f_teqm = interp1d(
                    money_values, TTB_max_values,
                    kind='linear', bounds_error=False, fill_value="extrapolate"
                )
                eqm_limit_x = f_teqm(upper_limit_y)
            else:
                # 直接使用缓存值
                mid_path = TTB_max_delta_values
                eqm_limit_x = TTB_max_values
                
            print(f"成功加载 zeta={zeta} 的 TTB_max_delta(mid_path) 和 TTB_max(eqm_limit_x) 数据")
            return mid_path, eqm_limit_x
        except Exception as e:
            print(f"从缓存加载数据时出错: {e}")
            print("计算默认值...")
    
    # 文件不存在或加载失败，使用默认值
    print(f"找不到zeta={zeta}的数据文件，使用默认计算值")
    mid_path = upper_limit_x  # 默认使用传入的upper_limit_x作为mid_path(TTB_max_delta)
    eqm_limit_x = calculate_equilibrium_line(left_boundary_x, upper_limit_x, upper_limit_y)
    
    print(f"注意：未找到由plot_path_costs.py生成的数据文件，请先运行该脚本生成tmax_teqm_zeta{zeta}.json")
    return mid_path, eqm_limit_x

# ========================= 散点数据生成/加载功能 =========================

def get_or_generate_scatter_points(mask_reg2, mask_reg, mask_reg3, mask_eqm, F1_GRID, F2_GRID, 
                                 e_value, num_points=3, results_dir='results'):
    """
    获取或生成符合各区域约束的散点数据，并保存到JSON文件
    如果JSON文件中某个区域没有散点数据，会尝试重新生成该区域的散点
    
    注意：由于区域存在包含关系(S⊃BS⊃RBS⊃T_eqm)，将使用区域差集确保各区域散点不重叠
    
    参数:
        mask_reg2: S_0^ζ区域掩码
        mask_reg: BS_0^ζ区域掩码
        mask_reg3: RBS_0^ζ区域掩码
        mask_eqm: T_eqm区域掩码
        F1_GRID: f1网格
        F2_GRID: f2网格
        e_value: 误差上限值
        num_points: 每个区域的点数
        results_dir: 结果存储目录
        
    返回:
        Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]: 四种区域的散点数据
    """
    # 创建存储目录
    os.makedirs(results_dir, exist_ok=True)
    cache_file = os.path.join(results_dir, f"scatter_points_e{int(e_value)}.json")
    
    # 计算区域差集，确保各区域不重叠
    # T_eqm: 使用原始mask_eqm (紫色区域)
    # RBS_0^ζ: 使用mask_reg3 - mask_eqm (绿色区域减去紫色区域)
    # BS_0^ζ: 使用mask_reg - mask_reg3 (蓝色区域减去绿色区域)
    # S_0^ζ: 使用mask_reg2 - mask_reg (红/黄色区域减去蓝色区域)
    mask_rs_only = np.logical_and(mask_reg3, np.logical_not(mask_eqm))
    mask_bs_only = np.logical_and(mask_reg, np.logical_not(mask_reg3))
    mask_s_only = np.logical_and(mask_reg2, np.logical_not(mask_reg))
    
    # 准备各区域的索引，供后续可能需要的重新生成使用
    eqm_indices = np.where(mask_eqm)
    rs_only_indices = np.where(mask_rs_only)
    bs_only_indices = np.where(mask_bs_only)
    s_only_indices = np.where(mask_s_only)
    
    # 定义空的散点数组
    s_constraint_points = np.array([])
    path_constraint_points = np.array([])
    all_constraint_points = np.array([])
    tmax_constraint_points = np.array([])
    
    # 标记是否需要更新JSON文件
    need_update = False
    
    # 检查是否存在缓存文件
    if os.path.exists(cache_file):
        print(f"正在从 {cache_file} 加载散点数据")
        try:
            with open(cache_file, 'r') as f:
                data = json.load(f)
            
            # 从JSON加载数据
            s_constraint_points = np.array(data.get('s_constraint_points', []))
            path_constraint_points = np.array(data.get('path_constraint_points', []))
            all_constraint_points = np.array(data.get('all_constraint_points', []))
            tmax_constraint_points = np.array(data.get('tmax_constraint_points', []))
            
            print(f"成功加载散点数据")
            
            # 检查各区域是否有数据，如果没有则尝试重新生成
            
            # 1. 检查S_0^ζ区域 (仅S区域，不包括其子区域)
            if s_constraint_points.size == 0 and len(s_only_indices[0]) > 0:
                print("S_0^ζ区域没有散点数据，尝试重新生成...")
                points_count = min(num_points, len(s_only_indices[0]))
                if points_count > 0:
                    s_constraint_indices = np.random.choice(len(s_only_indices[0]), points_count, replace=False)
                    s_constraint_points = np.array([
                        [F1_GRID[s_only_indices[0][i], s_only_indices[1][i]], 
                        F2_GRID[s_only_indices[0][i], s_only_indices[1][i]]]
                        for i in s_constraint_indices
                    ])
                    print(f"已为S_0^ζ区域生成 {len(s_constraint_points)} 个散点")
                    need_update = True
                else:
                    print("S_0^ζ区域没有足够的点用于生成散点")
            
            # 2. 检查BS_0^ζ区域 (仅BS区域，不包括其子区域)
            if all_constraint_points.size == 0 and len(bs_only_indices[0]) > 0:
                print("BS_0^ζ区域没有散点数据，尝试重新生成...")
                points_count = min(num_points, len(bs_only_indices[0]))
                if points_count > 0:
                    all_constraint_indices = np.random.choice(len(bs_only_indices[0]), points_count, replace=False)
                    all_constraint_points = np.array([
                        [F1_GRID[bs_only_indices[0][i], bs_only_indices[1][i]], 
                        F2_GRID[bs_only_indices[0][i], bs_only_indices[1][i]]]
                        for i in all_constraint_indices
                    ])
                    print(f"已为BS_0^ζ区域生成 {len(all_constraint_points)} 个散点")
                    need_update = True
                else:
                    print("BS_0^ζ区域没有足够的点用于生成散点")
            
            # 3. 检查RBS_0^ζ区域 (仅RBS区域，不包括其子区域)
            if path_constraint_points.size == 0 and len(rs_only_indices[0]) > 0:
                print("RBS_0^ζ区域没有散点数据，尝试重新生成...")
                points_count = min(num_points, len(rs_only_indices[0]))
                if points_count > 0:
                    path_constraint_indices = np.random.choice(len(rs_only_indices[0]), points_count, replace=False)
                    path_constraint_points = np.array([
                        [F1_GRID[rs_only_indices[0][i], rs_only_indices[1][i]], 
                        F2_GRID[rs_only_indices[0][i], rs_only_indices[1][i]]]
                        for i in path_constraint_indices
                    ])
                    print(f"已为RBS_0^ζ区域生成 {len(path_constraint_points)} 个散点")
                    need_update = True
                else:
                    print("RBS_0^ζ区域没有足够的点用于生成散点")
            
            # 4. 检查T_eqm区域
            if tmax_constraint_points.size == 0 and len(eqm_indices[0]) > 0:
                print("T_eqm区域没有散点数据，尝试重新生成...")
                points_count = min(num_points, len(eqm_indices[0]))
                if points_count > 0:
                    tmax_constraint_indices = np.random.choice(len(eqm_indices[0]), points_count, replace=False)
                    tmax_constraint_points = np.array([
                        [F1_GRID[eqm_indices[0][i], eqm_indices[1][i]], 
                        F2_GRID[eqm_indices[0][i], eqm_indices[1][i]]]
                        for i in tmax_constraint_indices
                    ])
                    print(f"已为T_eqm区域生成 {len(tmax_constraint_points)} 个散点")
                    need_update = True
                else:
                    print("T_eqm区域没有足够的点用于生成散点")
                    
            # 如果有新生成的散点，更新JSON文件
            if need_update:
                data = {
                    's_constraint_points': s_constraint_points.tolist() if len(s_constraint_points) > 0 else [],
                    'path_constraint_points': path_constraint_points.tolist() if len(path_constraint_points) > 0 else [],
                    'all_constraint_points': all_constraint_points.tolist() if len(all_constraint_points) > 0 else [],
                    'tmax_constraint_points': tmax_constraint_points.tolist() if len(tmax_constraint_points) > 0 else [],
                    'description': f"散点数据，用于e={e_value}的区域绘制。部分数据于 {time.strftime('%Y-%m-%d %H:%M:%S')} 更新。"
                }
                
                with open(cache_file, 'w') as f:
                    json.dump(data, f, indent=2)
                
                print(f"已更新散点数据到 {cache_file}")
                
            return s_constraint_points, path_constraint_points, all_constraint_points, tmax_constraint_points
        except Exception as e:
            print(f"从缓存加载散点数据时出错: {e}")
            print("生成新的散点数据...")
    
    # 如果没有缓存文件或加载失败，生成所有新的散点数据
    print(f"为 e={e_value} 生成新的散点数据")
    
    # 各区域独立生成散点，互不影响，按从小到大顺序处理
    # 3. T_eqm 区域 (紫色区域)
    tmax_constraint_points = np.array([])
    if len(eqm_indices[0]) > 0:
        points_count = min(num_points, len(eqm_indices[0]))
        if points_count > 0:
            tmax_constraint_indices = np.random.choice(len(eqm_indices[0]), points_count, replace=False)
            tmax_constraint_points = np.array([
                [F1_GRID[eqm_indices[0][i], eqm_indices[1][i]], 
                F2_GRID[eqm_indices[0][i], eqm_indices[1][i]]]
                for i in tmax_constraint_indices
            ])
            print(f"已为T_eqm区域生成 {len(tmax_constraint_points)} 个散点")
        else:
            print("T_eqm区域没有足够的点用于生成散点")
    else:
        print("T_eqm区域为空，无法生成散点")
        
    # 2. RBS_0^ζ 区域 (绿色区域，不含紫色部分)
    path_constraint_points = np.array([])
    if len(rs_only_indices[0]) > 0:
        points_count = min(num_points, len(rs_only_indices[0]))
        if points_count > 0:
            path_constraint_indices = np.random.choice(len(rs_only_indices[0]), points_count, replace=False)
            path_constraint_points = np.array([
                [F1_GRID[rs_only_indices[0][i], rs_only_indices[1][i]], 
                F2_GRID[rs_only_indices[0][i], rs_only_indices[1][i]]]
                for i in path_constraint_indices
            ])
            print(f"已为RBS_0^ζ区域(不含T_eqm)生成 {len(path_constraint_points)} 个散点")
        else:
            print("RBS_0^ζ区域(不含T_eqm)没有足够的点用于生成散点")
    else:
        print("RBS_0^ζ区域(不含T_eqm)为空，无法生成散点")
    
    # 1. BS_0^ζ 区域 (蓝色区域，不含绿色部分)
    all_constraint_points = np.array([])
    if len(bs_only_indices[0]) > 0:
        points_count = min(num_points, len(bs_only_indices[0]))
        if points_count > 0:
            all_constraint_indices = np.random.choice(len(bs_only_indices[0]), points_count, replace=False)
            all_constraint_points = np.array([
                [F1_GRID[bs_only_indices[0][i], bs_only_indices[1][i]], 
                F2_GRID[bs_only_indices[0][i], bs_only_indices[1][i]]]
                for i in all_constraint_indices
            ])
            print(f"已为BS_0^ζ区域(不含RBS_0^ζ)生成 {len(all_constraint_points)} 个散点")
        else:
            print("BS_0^ζ区域(不含RBS_0^ζ)没有足够的点用于生成散点")
    else:
        print("BS_0^ζ区域(不含RBS_0^ζ)为空，无法生成散点")
    
    # 0. S_0^ζ 区域 (红/黄色区域，不含蓝色部分)
    s_constraint_points = np.array([])
    if len(s_only_indices[0]) > 0:
        points_count = min(num_points, len(s_only_indices[0]))
        if points_count > 0:
            s_constraint_indices = np.random.choice(len(s_only_indices[0]), points_count, replace=False)
            s_constraint_points = np.array([
                [F1_GRID[s_only_indices[0][i], s_only_indices[1][i]], 
                F2_GRID[s_only_indices[0][i], s_only_indices[1][i]]]
                for i in s_constraint_indices
            ])
            print(f"已为S_0^ζ区域(不含BS_0^ζ)生成 {len(s_constraint_points)} 个散点")
        else:
            print("S_0^ζ区域(不含BS_0^ζ)没有足够的点用于生成散点")
    else:
        print("S_0^ζ区域(不含BS_0^ζ)为空，无法生成散点")
    
    # 将数据保存为JSON格式
    import time  # 添加时间模块，记录更新时间
    data = {
        's_constraint_points': s_constraint_points.tolist() if len(s_constraint_points) > 0 else [],
        'path_constraint_points': path_constraint_points.tolist() if len(path_constraint_points) > 0 else [],
        'all_constraint_points': all_constraint_points.tolist() if len(all_constraint_points) > 0 else [],
        'tmax_constraint_points': tmax_constraint_points.tolist() if len(tmax_constraint_points) > 0 else [],
        'description': f"散点数据，用于e={e_value}的区域绘制。生成时间: {time.strftime('%Y-%m-%d %H:%M:%S')}"
    }
    
    with open(cache_file, 'w') as f:
        json.dump(data, f, indent=2)
    
    print(f"已保存散点数据到 {cache_file}")
    return s_constraint_points, path_constraint_points, all_constraint_points, tmax_constraint_points


def _prepare_plot_data(e_val: float, results_dir: str = "results"):
    """根据给定的 epsilon 值计算和更新绘图所需的全局变量。"""
    global E, ZETA_VALUE
    global mask_reg, mask_reg2, mask_reg3, mask_eqm
    global P1_VALS, P2_VALS, P5_VALS
    global mid_path, eqm_limit_x, left_boundary_x
    global s_constraint_points, path_constraint_points, all_constraint_points, tmax_constraint_points

    E = e_val
    ZETA_VALUE = e_val
    
    # 路径成本矩阵
    P1_VALS = path1(F1_GRID, F2_GRID)
    P2_VALS = path2(F1_GRID, F2_GRID)
    P5_VALS = path5(F1_GRID, F2_GRID)

    # 更新可行域掩码
    mask_reg = feasible_mask(F1_GRID, F2_GRID)
    if not np.any(mask_reg):
        raise RuntimeError(f"❌ 可行域为空；请放宽 E 或其他排序约束。")

    # 只考虑可行域内值
    P1_INSIDE = P1_VALS[mask_reg]
    P2_INSIDE = P2_VALS[mask_reg]
    P5_INSIDE = P5_VALS[mask_reg]

    r1_max, r1_min = P1_INSIDE.max(), P1_INSIDE.min()
    r2_max, r2_min = P2_INSIDE.max(), P2_INSIDE.min()
    r5_max, r5_min = P5_INSIDE.max(), P5_INSIDE.min()

    # 计算左边界
    left_boundary_x = np.array([r1_min, r2_min, r5_min])
    default_upper_limit_x = np.array([(r1_max + r1_min) / 2, (r2_max + r2_min) / 2, (r5_max + r5_min) / 2])
    upper_limit_y = np.array([20.0, 15.0, 2.0])

    # 加载或计算路径限制
    mid_path, eqm_limit_x = get_or_compute_path_limits(
        left_boundary_x=left_boundary_x,
        upper_limit_x=default_upper_limit_x,
        upper_limit_y=upper_limit_y,
        e_value=E,
        results_dir=results_dir,
        zeta=ZETA_VALUE
    )

    # 输出计算结果
    print(f"左边界X: {left_boundary_x}")
    print(f"中点路径X (TTB_max_delta): {mid_path}")
    print(f"上限Y (货币成本): {upper_limit_y}")
    print(f"平衡线X (TTB_max): {eqm_limit_x}")
    print(f"使用的zeta值: {ZETA_VALUE}")

    # 更新其他区域掩码
    mask_reg3 = mask_reg & ((P1_VALS <= mid_path[0]) & (P2_VALS <= mid_path[1]) & (P5_VALS <= mid_path[2]))
    mask_reg2 = _box_constraints(F1_GRID, F2_GRID) & ((max_path(F1_GRID, F2_GRID) - min_path(F1_GRID, F2_GRID)) <= E)
    mask_eqm = mask_reg & ((P1_VALS <= eqm_limit_x[0]) & (P2_VALS <= eqm_limit_x[1]) & (P5_VALS <= eqm_limit_x[2]))

    # 获取或生成散点数据
    s_constraint_points, path_constraint_points, all_constraint_points, tmax_constraint_points = get_or_generate_scatter_points(
        mask_reg2=mask_reg2,
        mask_reg=mask_reg,
        mask_reg3=mask_reg3,
        mask_eqm=mask_eqm,
        F1_GRID=F1_GRID,
        F2_GRID=F2_GRID,
        e_value=E,
        num_points=2,
        results_dir=results_dir
    )

def _setup_plot_globals():
    """设置全局绘图样式、颜色和标记。"""
    global color_S0, color_BS0, color_RBS0, color_TTBmax
    global marker_S0, marker_BS0, marker_RBS0, marker_TTBmax

    plt.rcParams.update({
        "font.family": "serif", "font.serif": ["Times New Roman"], "font.size": 11,
        "text.usetex": True, "text.latex.preamble": r"\usepackage{amsmath}\usepackage{amssymb}",
        "axes.linewidth": 0.8, "axes.labelsize": 15, "axes.titlesize": 15,
        "axes.grid": True, "grid.linestyle": ":", "grid.linewidth": 0.5, "grid.alpha": 0.3,
        "xtick.direction": "in", "ytick.direction": "in",
        "xtick.major.size": 3.0, "ytick.major.size": 3.0,
        "xtick.minor.size": 1.5, "ytick.minor.size": 1.5,
        "xtick.major.width": 0.8, "ytick.major.width": 0.8,
        "xtick.minor.width": 0.6, "ytick.minor.width": 0.6,
        "legend.frameon": True, "legend.framealpha": 0.8, "legend.edgecolor": "k",
        "legend.fancybox": False, "legend.fontsize": 15,
        "savefig.dpi": 600, "savefig.format": "pdf", "savefig.bbox": "tight", "savefig.pad_inches": 0.05,
    })

    color_S0 = REGION_STYLES["S"].color
    color_BS0 = REGION_STYLES["BS"].color
    color_RBS0 = REGION_STYLES["RBS"].color
    color_TTBmax = REGION_STYLES["TTB_max"].color

    marker_S0 = REGION_STYLES["S"].marker
    marker_BS0 = REGION_STYLES["BS"].marker
    marker_RBS0 = REGION_STYLES["RBS"].marker
    marker_TTBmax = REGION_STYLES["TTB_max"].marker
    
    print(f"color_S0: {color_S0}, marker_S0: {marker_S0}")
    print(f"color_BS0: {color_BS0}, marker_BS0: {marker_BS0}")
    print(f"color_RBS0: {color_RBS0}, marker_RBS0: {marker_RBS0}")
    print(f"color_TTBmax: {color_TTBmax}, marker_TTBmax: {marker_TTBmax}")

# --------------------------- 创建绘图函数 ---------------------------

def create_plot(
    plot_num, show_reg2=True, show_reg=True, show_reg3=False, show_TTBmax=False, show_points=False, results_dir='results', plot_specific_points=None
):
    """创建不同的区域图
    
    参数:
        plot_num: 图的编号(1-3)
        show_reg2: 是否显示S_0^ζ区域
        show_reg: 是否显示BS_0^ζ区域
        show_reg3: 是否显示RBS_0^ζ区域
        show_TTBmax: 是否显示TTB_max区域
        show_points: 是否显示散点
        plot_specific_points: 指定绘制的散点区域，可选值为['s', 'rs', 'bs', 'TTB_max']
    """
    # 创建图形和轴
    fig, ax = plt.subplots(figsize=(8, 6), constrained_layout=True)
    
    # 定义用于图例的元素列表
    legend_elements = []
    
    # 绘制各个区域
    if show_reg2:
        ax.contourf(
            F1_GRID, F2_GRID, mask_reg2, levels=[0.5, 1.5], colors=[color_S0], alpha=0.7, zorder=1
        )
        ax.contour(F1_GRID, F2_GRID, mask_reg2.astype(int), levels=[0.5], colors="k", linewidths=0.6, zorder=10)
        legend_elements.append(Patch(facecolor=color_S0, edgecolor=color_S0, alpha=0.7, label=r"$S_0^\varepsilon$"))
    
    if show_reg:
        ax.contourf(
            F1_GRID, F2_GRID, mask_reg, levels=[0.5, 1.5], colors=[color_BS0], alpha=0.7, zorder=2
        )
        ax.contour(F1_GRID, F2_GRID, mask_reg.astype(int),  levels=[0.5], colors="k", linewidths=0.6, zorder=10)
        legend_elements.append(Patch(facecolor=color_BS0, edgecolor=color_BS0, alpha=0.7, label=r"$BS_0^\varepsilon$"))
    
    if show_reg3:
        ax.contourf(
            F1_GRID, F2_GRID, mask_reg3, levels=[0.5, 1.5], colors=[color_RBS0], alpha=0.7, zorder=3
        )
        ax.contour(F1_GRID, F2_GRID, mask_reg3.astype(int), levels=[0.5], colors="k", linewidths=0.6, zorder=10)
        legend_elements.append(Patch(facecolor=color_RBS0, edgecolor=color_RBS0, alpha=0.7, label=r"$RBS_0^\varepsilon-TTB_{max}+\delta_k$"))
    
    if show_TTBmax:
        ax.contourf(
            F1_GRID, F2_GRID, mask_eqm, levels=[0.5, 1.5], colors=[color_TTBmax], alpha=0.7, zorder=4
        )
        ax.contour(F1_GRID, F2_GRID, mask_eqm.astype(int),  levels=[0.5], colors="k", linewidths=0.6, zorder=10)
        legend_elements.append(Patch(facecolor=color_TTBmax, edgecolor=color_TTBmax, alpha=0.7, label=r"$RBS_0^\varepsilon-TTB_{max}$"))
    
    # 绘制散点
    point_label_idx = 0  # 用于标记点的索引
    point_labels = ['D', 'C', 'A', 'B']  # 点的标签
    
    if show_points:
        # 新增: 绘制S_0^ζ区域的散点
        if s_constraint_points.size > 0:
            if not plot_specific_points or 's' in plot_specific_points:
                s_pts = s_constraint_points[0:1] if plot_specific_points and len(plot_specific_points) == 1 else s_constraint_points
                # 绘制散点
                ax.scatter(
                    s_pts[:, 0],
                    s_pts[:, 1],
                    color=color_S0,
                    marker=marker_S0,
                    s=120,  # 增大标记尺寸以容纳字母
                    linewidth=0.8,
                    edgecolor='k',
                    label=r"Flows of $S_0^\varepsilon$",
                    zorder=15,  # 确保点在区域之上
                    alpha=0.9,  # 稍微增加透明度使文字更清晰
                )
                # 在标记内部添加字母
                if plot_num == 6:  # 只在fig6中添加标签
                    for i, (x, y) in enumerate(s_pts):
                        if point_label_idx < len(point_labels):
                            ax.annotate(
                                point_labels[point_label_idx],
                                xy=(x, y),
                                xytext=(8, 8),  # 相对于标记的偏移（点）
                                textcoords='offset points',
                                fontsize=13,
                                weight='bold',
                                color='k',
                                ha='left',
                                va='top',
                                zorder=20,
                                bbox=dict(boxstyle="round,pad=0.1", fc="white", alpha=0.8, ec="none"),
                            )
                            point_label_idx += 1
            
        if path_constraint_points.size > 0:
            if not plot_specific_points or 'rs' in plot_specific_points:
                rs_pts = path_constraint_points[0:1] if plot_specific_points and len(plot_specific_points) == 1 else path_constraint_points
                ax.scatter(
                    rs_pts[:, 0],
                    rs_pts[:, 1],
                    color=color_RBS0,
                    marker=marker_RBS0,
                    s=120,  # 增大标记尺寸以容纳字母
                    linewidth=0.8,
                    edgecolor='k',
                    label=r"Flows of $RBS_0^\varepsilon$",
                    zorder=15,  # 确保点在区域之上
                    alpha=0.9,  # 稍微增加透明度使文字更清晰
                )
                # 在标记内部添加字母
                if plot_num == 6:  # 只在fig6中添加标签
                    for i, (x, y) in enumerate(rs_pts):
                        if point_label_idx < len(point_labels):
                            ax.annotate(
                                point_labels[point_label_idx],
                                xy=(x, y),
                                xytext=(8, 8),
                                textcoords='offset points',
                                fontsize=13,
                                weight='bold',
                                color='k',
                                ha='left',
                                va='bottom',
                                zorder=20,
                                bbox=dict(boxstyle="round,pad=0.2", fc="white", alpha=0.8, ec="none"),
                            )
                            point_label_idx += 1
                
        if all_constraint_points.size > 0:
            if not plot_specific_points or 'bs' in plot_specific_points:
                bs_pts = all_constraint_points[0:1] if plot_specific_points and len(plot_specific_points) == 1 else all_constraint_points
                ax.scatter(
                    bs_pts[:, 0],
                    bs_pts[:, 1],
                    color=color_BS0,
                    marker=marker_BS0,
                    s=120,  # 增大标记尺寸以容纳字母
                    linewidth=0.8,
                    edgecolor='k',
                    label=r"Flows of $BS_0^\varepsilon$",
                    zorder=15,
                    alpha=0.9,  # 稍微增加透明度使文字更清晰
                )
                # 在标记内部添加字母
                if plot_num == 6:  # 只在fig6中添加标签
                    for i, (x, y) in enumerate(bs_pts):
                        if point_label_idx < len(point_labels):
                            ax.annotate(
                                point_labels[point_label_idx],
                                xy=(x, y),
                                xytext=(8, 8),
                                textcoords='offset points',
                                fontsize=13,
                                weight='bold',
                                color='k',
                                ha='left',
                                va='bottom',
                                zorder=20,
                                bbox=dict(boxstyle="round,pad=0.1", fc="white", alpha=0.8, ec="none"),
                            )
                            point_label_idx += 1
        
        if tmax_constraint_points.size > 0:
            if not plot_specific_points or 'TTB_max' in plot_specific_points:
                TTB_max_pts = tmax_constraint_points[0:1] if plot_specific_points and len(plot_specific_points) == 1 else tmax_constraint_points
                ax.scatter(
                    TTB_max_pts[:, 0],
                    TTB_max_pts[:, 1],
                    color="orange",
                    marker=marker_TTBmax,
                    s=120,  # 增大标记尺寸以容纳字母
                    linewidth=0.8,
                    edgecolor='k',
                    label=r"Flows of $TTB_{max}$",
                    zorder=15,
                    alpha=0.9,  # 稍微增加透明度使文字更清晰
                )
                # 在标记内部添加字母
                if plot_num == 6:  # 只在fig6中添加标签
                    for i, (x, y) in enumerate(TTB_max_pts):
                        if point_label_idx < len(point_labels):
                            ax.annotate(
                                point_labels[point_label_idx],
                                xy=(x, y),
                                xytext=(8, 8),
                                textcoords='offset points',
                                fontsize=13,
                                weight='bold',
                                color='k',
                                ha='left',
                                va='bottom',
                                zorder=20,
                                bbox=dict(boxstyle="round,pad=0.2", fc="white", alpha=0.8, ec="none"),
                            )
                            point_label_idx += 1
                if plot_num  == 5:
                    ax.annotate(
                        "A",
                        xy=(TTB_max_pts[0, 0], TTB_max_pts[0, 1]),
                        xytext=(8, 8),
                        textcoords='offset points',
                        fontsize=13,
                        weight='bold',
                        color='k',
                        ha='left',
                        va='bottom',
                        zorder=20,
                        bbox=dict(boxstyle="round,pad=0.2", fc="white", alpha=0.8, ec="none"),
                    )
    
    # 设置轴标签和范围
    ax.set_xlabel(r"$f_1$")
    ax.set_ylabel(r"$f_2$")
    ax.set_xlim(F1_MIN, F1_MAX)
    ax.set_ylim(F2_MIN, F2_MAX)
    
    # 添加科学计数法格式的刻度标签
    ax.ticklabel_format(axis='both', style='sci', scilimits=(0,0), useMathText=True)
    
    # 优化图例位置和样式
    if legend_elements:
        ax.legend(
            handles=legend_elements,
            loc="upper right", 
            fontsize=15,
            borderpad=0.4,
            labelspacing=0.3,
            handletextpad=0.5,
        )
    
    # 添加网格线
    ax.grid(True, linestyle=':', alpha=0.3, color='gray', linewidth=0.5)
    
    # 设置图表标题
    plot_titles = {
        1: " ",
        2: " ",
        3: " "
    }
    # ax.set_title(plot_titles.get(plot_num, f"图 {plot_num}"))
    
    # 输出PDF格式
    filename = f"region_path_cost_e{int(E)}_plot{plot_num}.pdf"
    out_path = os.path.join(results_dir, filename)
    plt.savefig(out_path, dpi=600, bbox_inches='tight')
    plt.close()

    print(f"✅ 图 {plot_num} 绘图完成，已保存至 {out_path}")
    
    # # 输出PNG格式
    # png_path = Path(__file__).with_name(f"ee_time_point{int(E)}_plot{plot_num}.png")
    # plt.savefig(png_path, dpi=300, format='png', bbox_inches='tight')
    # print(f"✅ 图 {plot_num} 额外输出PNG格式，已保存至 {png_path.resolve()}")

    
    return fig

def plot_rbs_point_region(F1_GRID, F2_GRID, mask_reg2, mask_reg, mask_reg3, mask_eqm,
                          path_constraint_points, color_S0, color_BS0, color_RBS0, color_TTBmax,
                          F1_MIN, F1_MAX, F2_MIN, F2_MAX, E, results_dir):
    """
    绘制全部区域 + 绿色区域(RBS)的一个点（使用黄色标记为"D"）
    
    参数:
        F1_GRID, F2_GRID: 网格数据
        mask_reg2, mask_reg, mask_reg3, mask_eqm: 各区域掩码
        path_constraint_points: RBS区域的点
        color_S0, color_BS0, color_RBS0, color_TTBmax: 各区域颜色
        F1_MIN, F1_MAX, F2_MIN, F2_MAX: 坐标轴范围
        E: epsilon值
        results_dir: 结果保存目录
    """
    fig7, ax7 = plt.subplots(figsize=(8, 6), constrained_layout=True)
    
    # 绘制各区域
    # S区域
    ax7.contourf(F1_GRID, F2_GRID, mask_reg2, levels=[0.5, 1.5], colors=[color_S0], alpha=0.7, zorder=1)
    ax7.contour(F1_GRID, F2_GRID, mask_reg2.astype(int), levels=[0.5], colors="k", linewidths=0.6, zorder=10)
    
    # BS区域
    ax7.contourf(F1_GRID, F2_GRID, mask_reg, levels=[0.5, 1.5], colors=[color_BS0], alpha=0.7, zorder=2)
    ax7.contour(F1_GRID, F2_GRID, mask_reg.astype(int),  levels=[0.5], colors="k", linewidths=0.6, zorder=10)
    
    # RBS区域
    ax7.contourf(F1_GRID, F2_GRID, mask_reg3, levels=[0.5, 1.5], colors=[color_RBS0], alpha=0.7, zorder=3)
    ax7.contour(F1_GRID, F2_GRID, mask_reg3.astype(int), levels=[0.5], colors="k", linewidths=0.6, zorder=10)
    
    # TTB_max区域
    ax7.contourf(F1_GRID, F2_GRID, mask_eqm, levels=[0.5, 1.5], colors=[color_TTBmax], alpha=0.7, zorder=4)
    ax7.contour(F1_GRID, F2_GRID, mask_eqm.astype(int),  levels=[0.5], colors="k", linewidths=0.6, zorder=10)
    
    # 绘制绿色区域(RBS)的点，使用黄色标记
    if path_constraint_points.size > 0:
        # 只取第一个点
        rs_pt = path_constraint_points[0:1]
        ax7.scatter(
            rs_pt[:, 0],
            rs_pt[:, 1],
            color='yellow',  # 黄色
            marker='D',      # 菱形标记
            s=120,
            linewidth=0.8,
            edgecolor='k',
            zorder=15,
            alpha=0.9,
        )
        # 添加"B"标签
        ax7.annotate(
            "B",
            xy=(rs_pt[0, 0], rs_pt[0, 1]),
            xytext=(8, 8),
            textcoords='offset points',
            fontsize=13,
            weight='bold',
            color='k',
            ha='left',
            va='bottom',
            zorder=20,
            bbox=dict(boxstyle="round,pad=0.2", fc="white", alpha=0.8, ec="none"),
        )
    
    # 添加图例
    legend_elements = [
        Patch(facecolor=color_S0, edgecolor=color_S0, alpha=0.7, label=r"$S_0^\varepsilon$"),
        Patch(facecolor=color_BS0, edgecolor=color_BS0, alpha=0.7, label=r"$BS_0^\varepsilon$"),
        Patch(facecolor=color_RBS0, edgecolor=color_RBS0, alpha=0.7, label=r"$RBS_0^\varepsilon-TTB_{max}+\delta_k$"),
        Patch(facecolor=color_TTBmax, edgecolor=color_TTBmax, alpha=0.7, label=r"$RBS_0^\varepsilon-TTB_{max}$"),
    ]
    
    ax7.legend(
        handles=legend_elements,
        loc="upper right", 
        fontsize=15,
        borderpad=0.4,
        labelspacing=0.3,
        handletextpad=0.5,
    )
    
    # 设置轴标签和范围
    ax7.set_xlabel(r"$f_1$")
    ax7.set_ylabel(r"$f_2$")
    ax7.set_xlim(F1_MIN, F1_MAX)
    ax7.set_ylim(F2_MIN, F2_MAX)
    ax7.ticklabel_format(axis='both', style='sci', scilimits=(0,0), useMathText=True)
    ax7.grid(True, linestyle=':', alpha=0.3, color='gray', linewidth=0.5)
    
    # 保存图7
    filename7 = f"region_path_cost_e{int(E)}_plot7.pdf"
    out_path7 = os.path.join(results_dir, filename7)
    plt.savefig(out_path7, dpi=600, bbox_inches='tight')
    plt.close()
    print(f"✅ 图 7 绘图完成，已保存至 {out_path7}")
    
    return rs_pt[0] if path_constraint_points.size > 0 else None


def plot_ttbmax_point_region(F1_GRID, F2_GRID, mask_reg2, mask_reg, mask_reg3, mask_eqm,
                             tmax_constraint_points, color_S0, color_BS0, color_RBS0, color_TTBmax,
                             F1_MIN, F1_MAX, F2_MIN, F2_MAX, E, results_dir):
    """
    绘制全部区域 + 紫色区域(TTB_max)的一个点（使用橙色标记为"D"）
    
    参数:
        F1_GRID, F2_GRID: 网格数据
        mask_reg2, mask_reg, mask_reg3, mask_eqm: 各区域掩码
        tmax_constraint_points: TTB_max区域的点
        color_S0, color_BS0, color_RBS0, color_TTBmax: 各区域颜色
        F1_MIN, F1_MAX, F2_MIN, F2_MAX: 坐标轴范围
        E: epsilon值
        results_dir: 结果保存目录
    """
    fig8, ax8 = plt.subplots(figsize=(8, 6), constrained_layout=True)
    
    # 绘制各区域
    # S区域
    ax8.contourf(F1_GRID, F2_GRID, mask_reg2, levels=[0.5, 1.5], colors=[color_S0], alpha=0.7, zorder=1)
    ax8.contour(F1_GRID, F2_GRID, mask_reg2.astype(int), levels=[0.5], colors="k", linewidths=0.6, zorder=10)
    
    # BS区域
    ax8.contourf(F1_GRID, F2_GRID, mask_reg, levels=[0.5, 1.5], colors=[color_BS0], alpha=0.7, zorder=2)
    ax8.contour(F1_GRID, F2_GRID, mask_reg.astype(int),  levels=[0.5], colors="k", linewidths=0.6, zorder=10)
    
    # RBS区域
    ax8.contourf(F1_GRID, F2_GRID, mask_reg3, levels=[0.5, 1.5], colors=[color_RBS0], alpha=0.7, zorder=3)
    ax8.contour(F1_GRID, F2_GRID, mask_reg3.astype(int), levels=[0.5], colors="k", linewidths=0.6, zorder=10)
    
    # TTB_max区域
    ax8.contourf(F1_GRID, F2_GRID, mask_eqm, levels=[0.5, 1.5], colors=[color_TTBmax], alpha=0.7, zorder=4)
    ax8.contour(F1_GRID, F2_GRID, mask_eqm.astype(int),  levels=[0.5], colors="k", linewidths=0.6, zorder=10)
    
    # 绘制紫色区域(TTB_max)的点，使用橙色标记
    if tmax_constraint_points.size > 0:
        # 只取第一个点
        tmax_pt = tmax_constraint_points[0:1]
        ax8.scatter(
            tmax_pt[:, 0],
            tmax_pt[:, 1],
            color='orange',  # 橙色
            marker='D',      # 菱形标记
            s=120,
            linewidth=0.8,
            edgecolor='k',
            zorder=15,
            alpha=0.9,
        )
        # 添加"A"标签
        ax8.annotate(
            "A",
            xy=(tmax_pt[0, 0], tmax_pt[0, 1]),
            xytext=(8, 8),
            textcoords='offset points',
            fontsize=13,
            weight='bold',
            color='k',
            ha='left',
            va='bottom',
            zorder=20,
            bbox=dict(boxstyle="round,pad=0.2", fc="white", alpha=0.8, ec="none"),
        )
    
    # 添加图例
    legend_elements = [
        Patch(facecolor=color_S0, edgecolor=color_S0, alpha=0.7, label=r"$S_0^\varepsilon$"),
        Patch(facecolor=color_BS0, edgecolor=color_BS0, alpha=0.7, label=r"$BS_0^\varepsilon$"),
        Patch(facecolor=color_RBS0, edgecolor=color_RBS0, alpha=0.7, label=r"$RBS_0^\varepsilon-TTB_{max}+\delta_k$"),
        Patch(facecolor=color_TTBmax, edgecolor=color_TTBmax, alpha=0.7, label=r"$RBS_0^\varepsilon-TTB_{max}$"),
    ]
    
    ax8.legend(
        handles=legend_elements,
        loc="upper right", 
        fontsize=15,
        borderpad=0.4,
        labelspacing=0.3,
        handletextpad=0.5,
    )
    
    # 设置轴标签和范围
    ax8.set_xlabel(r"$f_1$")
    ax8.set_ylabel(r"$f_2$")
    ax8.set_xlim(F1_MIN, F1_MAX)
    ax8.set_ylim(F2_MIN, F2_MAX)
    ax8.ticklabel_format(axis='both', style='sci', scilimits=(0,0), useMathText=True)
    ax8.grid(True, linestyle=':', alpha=0.3, color='gray', linewidth=0.5)
    
    # 保存图8
    filename8 = f"region_path_cost_e{int(E)}_plot8.pdf"
    out_path8 = os.path.join(results_dir, filename8)
    plt.savefig(out_path8, dpi=600, bbox_inches='tight')
    plt.close()
    print(f"✅ 图 8 绘图完成，已保存至 {out_path8}")
    
    return tmax_pt[0] if tmax_constraint_points.size > 0 else None


# ====================== 批量绘制不同 ε 的图 ======================


def generate_plots_for_e(e_val: int, results_dir: str = "results"):
    """根据给定 ε (zeta) 值重新计算区域并输出全部 6 张图。"""
    _prepare_plot_data(float(e_val), results_dir)

    print(f"\n====================  生成 ε = {e_val} 的区域图  ====================")

    # ---------- 绘制八张图 ----------
    create_plot(
        plot_num=1,
        show_reg2=True,
        show_reg=True,
        show_reg3=False,
        show_TTBmax=False,
        show_points=False,
        results_dir=results_dir,
    )
    create_plot(
        plot_num=2,
        show_reg2=True,
        show_reg=True,
        show_reg3=True,
        show_TTBmax=True,
        show_points=False,
        results_dir=results_dir,
    )
    create_plot(
        plot_num=3,
        show_reg2=True,
        show_reg=True,
        show_reg3=True,
        show_TTBmax=True,
        show_points=True,
        results_dir=results_dir,
    )
    create_plot(
        plot_num=4,
        show_reg2=True,
        show_reg=True,
        show_reg3=False,
        show_TTBmax=True,
        show_points=False,
        results_dir=results_dir,
    )
    create_plot(
        plot_num=5,
        show_reg2=True,
        show_reg=True,
        show_reg3=False,
        show_TTBmax=False,
        show_points=True,
        results_dir=results_dir,
        plot_specific_points=["TTB_max"],
    )
    create_plot(
        plot_num=6,
        show_reg2=True,
        show_reg=True,
        show_reg3=False,
        show_TTBmax=False,
        show_points=True,
        results_dir=results_dir,
        plot_specific_points=["s", "bs"],
    )
    # 调用封装好的函数绘制RBS和TTB_max点图
    plot_rbs_point_region(
        F1_GRID, F2_GRID, 
        mask_reg2, mask_reg, mask_reg3, mask_eqm,
        path_constraint_points, 
        color_S0, color_BS0, color_RBS0, color_TTBmax,
        F1_MIN, F1_MAX, F2_MIN, F2_MAX,
        E, results_dir
    )
    
    plot_ttbmax_point_region(
        F1_GRID, F2_GRID, 
        mask_reg2, mask_reg, mask_reg3, mask_eqm,
        tmax_constraint_points, 
        color_S0, color_BS0, color_RBS0, color_TTBmax,
        F1_MIN, F1_MAX, F2_MIN, F2_MAX,
        E, results_dir
    )

    print(f"✅ ε = {E} 的区域图全部生成完毕\n")

if __name__ == "__main__":
    # ============================== 基本常量 ==============================
    RHO: float = 15
    SIGMA: float = 0.02
    E: float = 24  # 误差上限 (与 Mathematica 变量 e 对应)
    EPS: float = 1e-8  # 严格不等式缓冲 (与 Mathematica 变量 eps 对应)

    # 边界条件 (与 box 变量对应)
    F1_MIN, F1_MAX = 1000.0, 8000.0
    F2_MIN, F2_MAX = 0.0, 7000.0
    TOTAL_MAX = 10000.0

    # --------------------------- 网格采样 ---------------------------

    STEP = 10  # 调节步长可提升/降低精度
    f1_vals = np.arange(F1_MIN, F1_MAX + STEP, STEP)
    f2_vals = np.arange(F2_MIN, F2_MAX + STEP, STEP)
    F1_GRID, F2_GRID = np.meshgrid(f1_vals, f2_vals, indexing="ij")

    mask_reg = feasible_mask(F1_GRID, F2_GRID)
    if not np.any(mask_reg):
        raise RuntimeError("❌ 可行域为空；请放宽 E 或其他排序约束。")
    
    # 为默认 E 值准备数据并设置绘图全局变量
    _prepare_plot_data(E)
    _setup_plot_globals()

    # ============================ 链路/路径函数 ===========================
    # 一次性生成多个 ε 值的所有结果图
    for _eps in [8, 16, 24, 32]:
        generate_plots_for_e(_eps)
    # generate_plots_for_e(32)