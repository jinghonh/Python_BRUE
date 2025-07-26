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
from pathlib import Path
from typing import Callable, Tuple, Dict, Any, List

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch
from scipy.interpolate import interp1d

# ============================== 基本常量 ==============================
RHO: float = 15
SIGMA: float = 0.02
E: float = 24  # 误差上限 (与 Mathematica 变量 e 对应)
EPS: float = 1e-8  # 严格不等式缓冲 (与 Mathematica 变量 eps 对应)

# 边界条件 (与 box 变量对应)
F1_MIN, F1_MAX = 1000.0, 8000.0
F2_MIN, F2_MAX = 0.0, 7000.0
TOTAL_MAX = 10000.0

# ============================ 链路/路径函数 ===========================

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

def get_or_compute_path_limits(left_boundary_x, upper_limit_x, upper_limit_y, e_value, results_dir='results'):
    """
    获取或计算mid_path和eqm_limit_x值，使用JSON格式存储以便手动编辑
    
    参数:
        left_boundary_x: 左边界X值
        upper_limit_x: 上限X值 (mid_path1, mid_path2, mid_path5)
        upper_limit_y: 上限Y值 (货币成本值)
        e_value: 误差上限值
        results_dir: 结果存储目录
        
    返回:
        Tuple[np.ndarray, np.ndarray]: mid_path和eqm_limit_x数组
    """
    # 创建存储目录
    os.makedirs(results_dir, exist_ok=True)
    cache_file = os.path.join(results_dir, f"path_limits_e{int(e_value)}.json")
    
    # 检查是否存在缓存文件
    if os.path.exists(cache_file):
        print(f"正在从 {cache_file} 加载 mid_path 和 eqm_limit_x 数据")
        try:
            with open(cache_file, 'r') as f:
                data = json.load(f)
            
            # 从JSON加载数据
            money_values = np.array(data['money_values'])
            mid_path_values = np.array(data['mid_path'])
            eqm_limit_values = np.array(data['eqm_limit_x'])
            
            # 如果货币值与缓存不同，需要进行插值
            if len(money_values) != len(upper_limit_y) or not np.allclose(money_values, upper_limit_y):
                # 为mid_path创建插值函数
                f_midpath = interp1d(
                    money_values, mid_path_values,
                    kind='linear', bounds_error=False, fill_value="extrapolate"
                )
                mid_path = f_midpath(upper_limit_y)
                
                # 为eqm_limit_x创建插值函数
                f_eqmlimit = interp1d(
                    money_values, eqm_limit_values,
                    kind='linear', bounds_error=False, fill_value="extrapolate"
                )
                eqm_limit_x = f_eqmlimit(upper_limit_y)
            else:
                # 直接使用缓存值
                mid_path = mid_path_values
                eqm_limit_x = eqm_limit_values
                
            print(f"成功加载 e={e_value} 的 mid_path 和 eqm_limit_x 数据")
            return mid_path, eqm_limit_x
        except Exception as e:
            print(f"从缓存加载 mid_path 和 eqm_limit_x 时出错: {e}")
            print("计算新值...")
    
    # 计算mid_path和eqm_limit_x
    print(f"为 e={e_value} 计算新的 mid_path 和 eqm_limit_x 值")
    mid_path = upper_limit_x
    eqm_limit_x = calculate_equilibrium_line(left_boundary_x, upper_limit_x, upper_limit_y)
    
    # 将数据保存为JSON格式（易于编辑）
    data = {
        'money_values': upper_limit_y.tolist(),
        'mid_path': mid_path.tolist(),
        'eqm_limit_x': eqm_limit_x.tolist(),
        'description': f"mid_path和eqm_limit_x数据，用于e={e_value}的区域绘制。可以手动编辑这些值以确保一致性。"
    }
    
    with open(cache_file, 'w') as f:
        json.dump(data, f, indent=2)
    
    print(f"已保存 mid_path 和 eqm_limit_x 到 {cache_file}")
    return mid_path, eqm_limit_x

# ========================= 散点数据生成/加载功能 =========================

def get_or_generate_scatter_points(mask_reg2, mask_reg, mask_reg3, mask_eqm, F1_GRID, F2_GRID, 
                                 e_value, num_points=3, results_dir='results'):
    """
    获取或生成符合各区域约束的散点数据，并保存到JSON文件
    如果JSON文件中某个区域没有散点数据，会尝试重新生成该区域的散点
    
    参数:
        mask_reg2: S_0^ζ区域掩码
        mask_reg: BS_0^ζ区域掩码
        mask_reg3: RS_0^ζ区域掩码
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
    
    # 准备各区域的索引，供后续可能需要的重新生成使用
    reg2_indices = np.where(mask_reg2)
    reg_indices = np.where(mask_reg)
    reg3_indices = np.where(mask_reg3)
    eqm_indices = np.where(mask_eqm)
    
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
            
            # 1. 检查S_0^ζ区域
            if s_constraint_points.size == 0 and len(reg2_indices[0]) > 0:
                print("S_0^ζ区域没有散点数据，尝试重新生成...")
                points_count = min(num_points, len(reg2_indices[0]))
                if points_count > 0:
                    s_constraint_indices = np.random.choice(len(reg2_indices[0]), points_count, replace=False)
                    s_constraint_points = np.array([
                        [F1_GRID[reg2_indices[0][i], reg2_indices[1][i]], 
                        F2_GRID[reg2_indices[0][i], reg2_indices[1][i]]]
                        for i in s_constraint_indices
                    ])
                    print(f"已为S_0^ζ区域生成 {len(s_constraint_points)} 个散点")
                    need_update = True
                else:
                    print("S_0^ζ区域没有足够的点用于生成散点")
            
            # 2. 检查RS_0^ζ区域
            if path_constraint_points.size == 0 and len(reg3_indices[0]) > 0:
                print("RS_0^ζ区域没有散点数据，尝试重新生成...")
                points_count = min(num_points, len(reg3_indices[0]))
                if points_count > 0:
                    path_constraint_indices = np.random.choice(len(reg3_indices[0]), points_count, replace=False)
                    path_constraint_points = np.array([
                        [F1_GRID[reg3_indices[0][i], reg3_indices[1][i]], 
                        F2_GRID[reg3_indices[0][i], reg3_indices[1][i]]]
                        for i in path_constraint_indices
                    ])
                    print(f"已为RS_0^ζ区域生成 {len(path_constraint_points)} 个散点")
                    need_update = True
                else:
                    print("RS_0^ζ区域没有足够的点用于生成散点")
            
            # 3. 检查BS_0^ζ区域
            if all_constraint_points.size == 0 and len(reg_indices[0]) > 0:
                print("BS_0^ζ区域没有散点数据，尝试重新生成...")
                points_count = min(num_points, len(reg_indices[0]))
                if points_count > 0:
                    all_constraint_indices = np.random.choice(len(reg_indices[0]), points_count, replace=False)
                    all_constraint_points = np.array([
                        [F1_GRID[reg_indices[0][i], reg_indices[1][i]], 
                        F2_GRID[reg_indices[0][i], reg_indices[1][i]]]
                        for i in all_constraint_indices
                    ])
                    print(f"已为BS_0^ζ区域生成 {len(all_constraint_points)} 个散点")
                    need_update = True
                else:
                    print("BS_0^ζ区域没有足够的点用于生成散点")
            
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
    
    # 各区域独立生成散点，互不影响
    # 0. S_0^ζ 区域 (黄/红色区域，s_constraint_points)
    s_constraint_points = np.array([])
    if len(reg2_indices[0]) > 0:
        points_count = min(num_points, len(reg2_indices[0]))
        if points_count > 0:
            s_constraint_indices = np.random.choice(len(reg2_indices[0]), points_count, replace=False)
            s_constraint_points = np.array([
                [F1_GRID[reg2_indices[0][i], reg2_indices[1][i]], 
                F2_GRID[reg2_indices[0][i], reg2_indices[1][i]]]
                for i in s_constraint_indices
            ])
            print(f"已为S_0^ζ区域生成 {len(s_constraint_points)} 个散点")
        else:
            print("S_0^ζ区域没有足够的点用于生成散点")
    else:
        print("S_0^ζ区域为空，无法生成散点")
        
    # 1. RS_0^ζ 区域 (绿色区域，path_constraint_points)
    path_constraint_points = np.array([])
    if len(reg3_indices[0]) > 0:
        points_count = min(num_points, len(reg3_indices[0]))
        if points_count > 0:
            path_constraint_indices = np.random.choice(len(reg3_indices[0]), points_count, replace=False)
            path_constraint_points = np.array([
                [F1_GRID[reg3_indices[0][i], reg3_indices[1][i]], 
                F2_GRID[reg3_indices[0][i], reg3_indices[1][i]]]
                for i in path_constraint_indices
            ])
            print(f"已为RS_0^ζ区域生成 {len(path_constraint_points)} 个散点")
        else:
            print("RS_0^ζ区域没有足够的点用于生成散点")
    else:
        print("RS_0^ζ区域为空，无法生成散点")
    
    # 2. BS_0^ζ 区域 (蓝色区域，all_constraint_points)
    all_constraint_points = np.array([])
    if len(reg_indices[0]) > 0:
        points_count = min(num_points, len(reg_indices[0]))
        if points_count > 0:
            all_constraint_indices = np.random.choice(len(reg_indices[0]), points_count, replace=False)
            all_constraint_points = np.array([
                [F1_GRID[reg_indices[0][i], reg_indices[1][i]], 
                F2_GRID[reg_indices[0][i], reg_indices[1][i]]]
                for i in all_constraint_indices
            ])
            print(f"已为BS_0^ζ区域生成 {len(all_constraint_points)} 个散点")
        else:
            print("BS_0^ζ区域没有足够的点用于生成散点")
    else:
        print("BS_0^ζ区域为空，无法生成散点")
    
    # 3. T_eqm 区域 (紫色区域，tmax_constraint_points)
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

# --------------------------- 网格采样 ---------------------------

STEP = 10  # 调节步长可提升/降低精度
f1_vals = np.arange(F1_MIN, F1_MAX + STEP, STEP)
f2_vals = np.arange(F2_MIN, F2_MAX + STEP, STEP)
F1_GRID, F2_GRID = np.meshgrid(f1_vals, f2_vals, indexing="ij")

mask_reg = feasible_mask(F1_GRID, F2_GRID)
if not np.any(mask_reg):
    raise RuntimeError("❌ 可行域为空；请放宽 E 或其他排序约束。")

# ---------------------- 计算 max/min/mid ----------------------

# 在采样网格上近似取最值（若需更精确可结合 SciPy 优化器）
P1_VALS = path1(F1_GRID, F2_GRID)
P2_VALS = path2(F1_GRID, F2_GRID)
P5_VALS = path5(F1_GRID, F2_GRID)

# 只考虑可行域内值
P1_INSIDE = P1_VALS[mask_reg]
P2_INSIDE = P2_VALS[mask_reg]
P5_INSIDE = P5_VALS[mask_reg]

r1_max, r1_min = P1_INSIDE.max(), P1_INSIDE.min()
r2_max, r2_min = P2_INSIDE.max(), P2_INSIDE.min()
r5_max, r5_min = P5_INSIDE.max(), P5_INSIDE.min()

# ---------------------- 计算 mid_path 和平衡线 ----------------------

# 计算左边界 - 使用最小值作为近似
left_boundary_x = np.array([r1_min, r2_min, r5_min])

# 设置默认的upper_limit_x和upper_limit_y
default_upper_limit_x = np.array([(r1_max + r1_min) / 2, (r2_max + r2_min) / 2, (r5_max + r5_min) / 2])
upper_limit_y = np.array([20.0, 15.0, 2.0])

# 使用JSON文件加载或计算mid_path和eqm_limit_x
results_dir = 'results'
mid_path, eqm_limit_x = get_or_compute_path_limits(
    left_boundary_x=left_boundary_x,
    upper_limit_x=default_upper_limit_x,
    upper_limit_y=upper_limit_y,
    e_value=E,
    results_dir=results_dir
)

# 输出计算结果
print(f"左边界X: {left_boundary_x}")
print(f"中点路径X: {mid_path}")
print(f"上限Y (货币成本): {upper_limit_y}")
print(f"平衡线X: {eqm_limit_x}")

# --------------------------- 其他区域 ---------------------------

mask_reg3 = mask_reg & (
    (P1_VALS <= mid_path[0]) & (P2_VALS <= mid_path[1]) & (P5_VALS <= mid_path[2])
)
mask_reg2 = _box_constraints(F1_GRID, F2_GRID) & (
    (max_path(F1_GRID, F2_GRID) - min_path(F1_GRID, F2_GRID)) <= E
)

# 平衡线区域 - Pi_VALS < eqmLimitX
mask_eqm = mask_reg & (
    (P1_VALS <= eqm_limit_x[0]) & 
    (P2_VALS <= eqm_limit_x[1]) & 
    (P5_VALS <= eqm_limit_x[2])
)

# 使用函数获取或生成散点数据
s_constraint_points, path_constraint_points, all_constraint_points, tmax_constraint_points = get_or_generate_scatter_points(
    mask_reg2=mask_reg2,
    mask_reg=mask_reg, 
    mask_reg3=mask_reg3, 
    mask_eqm=mask_eqm,
    F1_GRID=F1_GRID,
    F2_GRID=F2_GRID,
    e_value=E,
    num_points=2,  # 每个区域生成2个点
    results_dir=results_dir
)

# --------------------------- 设置绘图样式 ---------------------------

# 设置学术期刊风格
plt.rcParams.update({
    # 基础字体设置
    "font.family": "serif",
    "font.serif": ["Times New Roman"],
    "font.size": 11,
    
    # LaTeX设置
    "text.usetex": True,
    "text.latex.preamble": r"\usepackage{amsmath}\usepackage{amssymb}",
    
    # 图形风格
    "axes.linewidth": 0.8,
    "axes.labelsize": 11,
    "axes.titlesize": 12,
    "axes.grid": True,
    "grid.linestyle": ":",
    "grid.linewidth": 0.5,
    "grid.alpha": 0.3,
    
    # 刻度设置
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.major.size": 3.0,
    "ytick.major.size": 3.0,
    "xtick.minor.size": 1.5,
    "ytick.minor.size": 1.5,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "xtick.minor.width": 0.6,
    "ytick.minor.width": 0.6,
    
    # 图例设置
    "legend.frameon": True,
    "legend.framealpha": 0.8,
    "legend.edgecolor": "k",
    "legend.fancybox": False,
    "legend.fontsize": 9,
    
    # 保存设置
    "savefig.dpi": 600,
    "savefig.format": "pdf",
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.05,
})

# 定义一套柔和且清晰的学术颜色方案 (基于 ColorBrewer Set2) 别修改注释
color_S0 = "#ea9999"   # 更红的柔和红色 - S_0^ζ（替换原灰色，更适合大面积底色）
color_BS0 = "#1f77b4"  # 蓝色 - BS_0^ζ
color_RS0 = "#4daf4a"  # 绿色 - RS_0^ζ
color_Teqm = "#984ea3" # 紫色 - T_eqm

# --------------------------- 创建绘图函数 ---------------------------

def create_plot(plot_num, show_reg2=True, show_reg=True, show_reg3=False, show_eqm=False, show_points=False):
    """创建不同的区域图
    
    参数:
        plot_num: 图的编号(1-3)
        show_reg2: 是否显示S_0^ζ区域
        show_reg: 是否显示BS_0^ζ区域
        show_reg3: 是否显示RS_0^ζ区域
        show_eqm: 是否显示T_eqm区域
        show_points: 是否显示散点
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
        legend_elements.append(Patch(facecolor=color_S0, edgecolor=color_S0, alpha=0.7, label=r"$S_0^\zeta$"))
    
    if show_reg:
        ax.contourf(
            F1_GRID, F2_GRID, mask_reg, levels=[0.5, 1.5], colors=[color_BS0], alpha=0.7, zorder=2
        )
        ax.contour(F1_GRID, F2_GRID, mask_reg.astype(int),  levels=[0.5], colors="k", linewidths=0.6, zorder=10)
        legend_elements.append(Patch(facecolor=color_BS0, edgecolor=color_BS0, alpha=0.7, label=r"$BS_0^\zeta$"))
    
    if show_reg3:
        ax.contourf(
            F1_GRID, F2_GRID, mask_reg3, levels=[0.5, 1.5], colors=[color_RS0], alpha=0.7, zorder=3
        )
        ax.contour(F1_GRID, F2_GRID, mask_reg3.astype(int), levels=[0.5], colors="k", linewidths=0.6, zorder=10)
        legend_elements.append(Patch(facecolor=color_RS0, edgecolor=color_RS0, alpha=0.7, label=r"$RS_0^\zeta$"))
    
    if show_eqm:
        ax.contourf(
            F1_GRID, F2_GRID, mask_eqm, levels=[0.5, 1.5], colors=[color_Teqm], alpha=0.7, zorder=4
        )
        ax.contour(F1_GRID, F2_GRID, mask_eqm.astype(int),  levels=[0.5], colors="k", linewidths=0.6, zorder=10)
        legend_elements.append(Patch(facecolor=color_Teqm, edgecolor=color_Teqm, alpha=0.7, label=r"$T_{eqm}$"))
    
    # 绘制散点
    if show_points:
        # 新增: 绘制S_0^ζ区域的散点
        if s_constraint_points.size > 0:
            ax.scatter(
                s_constraint_points[:, 0],
                s_constraint_points[:, 1],
                color=color_S0,
                marker="^",  # 使用三角形标记，与其他区域区分开
                s=45,
                linewidth=0.8,
                edgecolor='k',
                label=r"Flows of $S_0^\zeta$",
                zorder=15,  # 确保点在区域之上
            )
            
        if path_constraint_points.size > 0:
            ax.scatter(
                path_constraint_points[:, 0],
                path_constraint_points[:, 1],
                color=color_RS0,
                marker="o",
                s=40,
                linewidth=0.8,
                edgecolor='k',
                label=r"Flows of $RS_0^\zeta$",
                zorder=15,  # 确保点在区域之上
            )
        if all_constraint_points.size > 0:
            ax.scatter(
                all_constraint_points[:, 0],
                all_constraint_points[:, 1],
                color=color_BS0,
                marker="s",
                s=40,
                linewidth=0.8,
                edgecolor='k',
                label=r"Flows of $BS_0^\zeta$",
                zorder=15,
            )
        if tmax_constraint_points.size > 0:
            ax.scatter(
                tmax_constraint_points[:, 0],
                tmax_constraint_points[:, 1],
                color=color_Teqm,
                marker="D",
                s=40,
                linewidth=0.8,
                edgecolor='k',
                label=r"Flows of $T_{eqm}$",
                zorder=15,
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
            fontsize=8,
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
    ax.set_title(plot_titles.get(plot_num, f"图 {plot_num}"))
    
    # 输出PDF格式
    out_path = Path(__file__).with_name(f"ee_time_point{int(E)}_plot{plot_num}.pdf")
    plt.savefig(out_path, dpi=600, bbox_inches='tight')
    print(f"✅ 图 {plot_num} 绘图完成，已保存至 {out_path.resolve()}")
    
    # 输出PNG格式
    png_path = Path(__file__).with_name(f"ee_time_point{int(E)}_plot{plot_num}.png")
    plt.savefig(png_path, dpi=300, format='png', bbox_inches='tight')
    print(f"✅ 图 {plot_num} 额外输出PNG格式，已保存至 {png_path.resolve()}")
    
    return fig

# --------------------------- 创建三个不同的图 ---------------------------

# 图1: 只绘制前两个区域
fig1 = create_plot(plot_num=1, show_reg2=True, show_reg=True, show_reg3=False, show_eqm=False, show_points=False)

# 图2: 绘制全部区域
fig2 = create_plot(plot_num=2, show_reg2=True, show_reg=True, show_reg3=True, show_eqm=True, show_points=False)

# 图3: 绘制全部区域和散点
fig3 = create_plot(plot_num=3, show_reg2=True, show_reg=True, show_reg3=True, show_eqm=True, show_points=True)

# 显示所有图形
plt.show()  