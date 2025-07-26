"""demo.py – Python 版的 BRUE 路径成本区域绘制示例
====================================================
该脚本根据给出的 Mathematica 代码 (见 README) 重写而成，
完成如下功能：
1. 定义链路/路径成本函数；
2. 构造三类可行域 (reg, reg2, reg3)；
3. 计算各路径成本的最值及中点；
4. 采样并绘制区域及给定的方案点；
5. 导出 pdf 文件 `ee_time_point{e}.pdf`。

与 Mathematica 的 `RegionPlot` 不同，Python 这里采用网格采样 +
布尔掩码的方式近似绘制区域。若需更高精度，可适当缩小步长。
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Callable, Tuple

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

# ============================== 基本常量 ==============================
RHO: float = 15
SIGMA: float = 0.02
E: float = 15.0  # 误差上限 (与 Mathematica 变量 e 对应)
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

mid_path1, mid_path2, mid_path5 = (r1_max + r1_min) / 2, (r2_max + r2_min) / 2, (
    r5_max + r5_min
) / 2

# ---------------------- 计算平衡线 ----------------------

# 设置upperLimitX和upperLimitY
upper_limit_x = np.array([mid_path1, mid_path2, mid_path5])
upper_limit_y = np.array([20.0, 15.0, 2.0])

# 计算左边界 - 使用最小值作为近似
left_boundary_x = np.array([r1_min, r2_min, r5_min])

# 计算平衡线
eqm_limit_x = calculate_equilibrium_line(left_boundary_x, upper_limit_x, upper_limit_y)

# 输出计算结果
print(f"左边界X: {left_boundary_x}")
print(f"上限X: {upper_limit_x}")
print(f"上限Y (货币成本): {upper_limit_y}")
print(f"平衡线X: {eqm_limit_x}")

# --------------------------- 其他区域 ---------------------------

mask_reg3 = mask_reg & (
    (P1_VALS <= mid_path1) & (P2_VALS <= mid_path2) & (P5_VALS <= mid_path5)
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

# --------------------------- 绘图 ---------------------------

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

# 定义学术期刊适合的颜色方案
color_S0 = "#0072B2"       # 深蓝色
color_BS0 = "#D55E00"      # 砖红色
color_RS0 = "#009E73"      # 深绿色
color_Teqm = "#CC79A7"     # 紫红色

# 创建图形和轴
fig, ax = plt.subplots(figsize=(4.5, 3.5), constrained_layout=True)

# 使用掩码散点可视化近似区域
ax.scatter(
    F1_GRID[mask_reg2],
    F2_GRID[mask_reg2],
    s=2.5,  # 减小标记大小
    color=color_S0,
    alpha=0.3,
    edgecolors="none",
    rasterized=True,  # 栅格化以减少PDF文件大小
    label=r"$S_0^\zeta$",
)
ax.scatter(
    F1_GRID[mask_reg],
    F2_GRID[mask_reg],
    s=2.5,
    color=color_BS0,
    alpha=0.3,
    edgecolors="none",
    rasterized=True,
    label=r"$BS_0^\zeta$",
)
ax.scatter(
    F1_GRID[mask_reg3],
    F2_GRID[mask_reg3],
    s=2.5,
    color=color_RS0,
    alpha=0.3,
    edgecolors="none",
    rasterized=True,
    label=r"$RS_0^\zeta$",
)
ax.scatter(
    F1_GRID[mask_eqm],
    F2_GRID[mask_eqm],
    s=2.5,
    color=color_Teqm,
    alpha=0.3,
    edgecolors="none",
    rasterized=True,
    label=r"$T_{eqm}$",
)

# 方案点 - 如果需要的话取消注释这段代码
# path_constraint_points = np.array([[5373.80, 4413.76], [4582.07, 3312.18]])
# all_constraint_points = np.array([[4758.09, 3882.45], [4493.79, 3621.15]])
# tmax_constraint_points = np.array([[4250.23, 3855.09], [4181.32, 4064.58]])

# ax.scatter(
#     path_constraint_points[:, 0],
#     path_constraint_points[:, 1],
#     color=color_RS0,
#     marker="o",
#     s=40,
#     linewidth=0.8,
#     edgecolor='k',
#     label=r"Flows of $S_0^\zeta$",
#     zorder=5,  # 确保点在散点之上
# )
# ax.scatter(
#     all_constraint_points[:, 0],
#     all_constraint_points[:, 1],
#     color=color_BS0,
#     marker="s",
#     s=40,
#     linewidth=0.8,
#     edgecolor='k',
#     label=r"Flows of $BS_0^\zeta$",
#     zorder=5,
# )
# ax.scatter(
#     tmax_constraint_points[:, 0],
#     tmax_constraint_points[:, 1],
#     color=color_Teqm,
#     marker="D",
#     s=40,
#     linewidth=0.8,
#     edgecolor='k',
#     label=r"Flows of $RS_0^\zeta$",
#     zorder=5,
# )

# 设置轴标签和范围
ax.set_xlabel(r"$f_1$")
ax.set_ylabel(r"$f_2$")
ax.set_xlim(F1_MIN, F1_MAX)
ax.set_ylim(F2_MIN, F2_MAX)

# 添加科学计数法格式的刻度标签
ax.ticklabel_format(axis='both', style='sci', scilimits=(0,0), useMathText=True)

# 优化图例位置和样式
ax.legend(
    loc="upper right", 
    fontsize=8,
    markerscale=2,  # 增大图例中的标记尺寸
    borderpad=0.4,
    labelspacing=0.3,
    handletextpad=0.5,
)

# 添加网格线
ax.grid(True, linestyle=':', alpha=0.3, color='gray', linewidth=0.5)

# 添加子图标签 (a), (b)等，如果是多子图的一部分
# ax.text(-0.15, 1.05, '(a)', transform=ax.transAxes, fontsize=11, fontweight='bold')

# 输出 pdf，使用严格的学术期刊标准
out_path = Path(__file__).with_name(f"ee_time_point{int(E)}_academic.pdf")
plt.savefig(out_path, dpi=600, bbox_inches='tight')
print(f"✅ 绘图完成，已保存至 {out_path.resolve()}")

# 如果需要PNG格式的输出
png_path = Path(__file__).with_name(f"ee_time_point{int(E)}_academic.png")
plt.savefig(png_path, dpi=300, format='png', bbox_inches='tight')
print(f"✅ 额外输出PNG格式图片，已保存至 {png_path.resolve()}")

# 显示图形
plt.show()
