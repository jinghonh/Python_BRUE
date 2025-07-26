#!/usr/bin/env python
# -*- coding: utf-8 -*-

import json
import matplotlib.pyplot as plt
import os
import sys
from pathlib import Path
import numpy as np
from dataclasses import dataclass
from typing import Tuple


@dataclass
class PlotParams:
    """Configuration parameters for plotting."""
    save_path: str = 'results/'
    font_name: str = 'serif'
    font_size: int = 10
    show_grid: bool = True
    figure_dpi: int = 300
    figure_size: Tuple[int, int] = (8, 6)


def load_tmax_teqm_data(zeta):
    """
    加载指定zeta值的t_max和t_eqm数据
    
    Args:
        zeta (int): zeta值
    
    Returns:
        dict: 包含money_values, t_max, t_eqm的字典
    """
    # 获取项目根目录
    current_dir = Path(__file__).parent
    project_root = current_dir.parent
    
    # 构建JSON文件路径
    json_path = project_root / "results" / f"tmax_teqm_zeta{zeta}.json"
    
    try:
        with open(json_path, 'r', encoding='utf-8') as f:
            data = json.load(f)
        return data
    except FileNotFoundError:
        print(f"错误：找不到zeta={zeta}的数据文件：{json_path}")
        return None


def configure_plot(ax, params, title, xlabel, ylabel):
    """Apply common plot configurations."""
    ax.set_title(title, fontsize=params.font_size + 4, fontweight='bold')
    ax.set_xlabel(xlabel, fontsize=params.font_size + 2, fontweight='bold')
    ax.set_ylabel(ylabel, fontsize=params.font_size + 2, fontweight='bold')
    if params.show_grid:
        ax.grid(True, which='major', alpha=0.2)
        ax.grid(True, which='minor', alpha=0.1)
    ax.tick_params(axis='both', which='major', labelsize=params.font_size)


def plot_tmax_teqm(zeta_values, params=None):
    """
    根据多个zeta值绘制t_max和t_eqm的折线图
    
    Args:
        zeta_values (list): zeta值列表
        params (PlotParams): 绘图参数
    """
    # 设置绘图样式
    plt.rcParams['text.usetex'] = True
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['font.serif'] = ['Computer Modern Roman', 'Times New Roman']
    plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'
    
    # 设置绘图参数
    if params is None:
        params = PlotParams()
    
    # 创建图表
    fig, ax = plt.subplots(figsize=params.figure_size)
    configure_plot(ax, params, r' ', r'Time ($T$)', r'Money Value ($M$)')
    
    # 为不同的zeta值定义不同的标记形状
    markers = ['o', 's', '^', 'v', 'd', '*', 'p', 'h', '+', 'x']
    
    # 定义颜色
    tmax_color = [0.8, 0.2, 0.2]  # 红色系
    teqm_color = [0.5, 0.0, 0.8]  # 紫色系

    # 为自定义图例创建句柄列表
    legend_handles = []
    
    # 绘制图表
    for i, zeta in enumerate(zeta_values):
        data = load_tmax_teqm_data(zeta)
        if data:
            # 选择当前zeta值对应的标记
            marker = markers[i % len(markers)]
            
            # 绘制T_max，使用实线 (-) 和红色系，添加标记
            ax.plot(data['t_max'], data['money_values'], '-', 
                   color=tmax_color, linewidth=2, 
                   marker=marker, markersize=8, markerfacecolor=tmax_color, markeredgecolor='black', markeredgewidth=0.5)
            
            # 绘制T_eqm，使用虚线 (--) 和紫色系，添加标记
            ax.plot(data['t_eqm'], data['money_values'], '--', 
                   color=teqm_color, linewidth=2, 
                   marker=marker, markersize=8, markerfacecolor=teqm_color, markeredgecolor='black', markeredgewidth=0.5)
            
            # 只为第一个zeta值创建zeta图例条目
            if i == 0:
                # 为T_max和T_eqm创建图例条目
                legend_handles.append(plt.Line2D([0], [0], color=tmax_color, lw=2, linestyle='-', 
                                              label=r'$T_{max}$'))
                legend_handles.append(plt.Line2D([0], [0], color=teqm_color, lw=2, linestyle='--', 
                                              label=r'$T_{eqm}$'))
            
            # 为每个zeta值创建图例条目，使用黑色以便区分
            legend_handles.append(plt.Line2D([0], [0], marker=marker, color='black', 
                                          label=r'$\varepsilon=' + str(zeta) + r'$',
                                          markerfacecolor='black', markersize=8))
    
    # 添加自定义图例
    ax.legend(handles=legend_handles, loc='best', fontsize=params.font_size, 
              framealpha=0.9, fancybox=True, shadow=True)
    
    # 保存图像
    os.makedirs(params.save_path, exist_ok=True)
    save_path = Path(params.save_path) / "tmax_teqm_comparison.pdf"
    fig.savefig(save_path, dpi=params.figure_dpi, bbox_inches='tight')
    print(f"图像已保存至: {save_path}")
    
    plt.show()


def main():
    """主函数，获取用户输入的zeta值并绘制图表"""
    # 设置默认绘图参数
    params = PlotParams(save_path='results/', figure_dpi=600)
    
    print("请输入要比较的zeta值（用空格分隔，例如：15 16 17）：")
    try:
        input_str = input().strip()
        zeta_values = [int(x) for x in input_str.split()]
        
        if not zeta_values:
            print("未输入zeta值，将使用默认值：16")
            zeta_values = [16]
        
        plot_tmax_teqm(zeta_values, params)
    
    except ValueError:
        print("输入格式错误，请输入以空格分隔的整数")
        sys.exit(1)


if __name__ == "__main__":
    main() 