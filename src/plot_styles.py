from dataclasses import dataclass

@dataclass(frozen=True)
class RegionStyle:
    """Store color, marker and LaTeX label information for a region."""
    color: str  # Hex color or matplotlib‐compatible color
    marker: str  # Marker symbol used in scatter/line plots
    label: str  # LaTeX label for legends


# 统一的区域样式表，可在此处集中调整
REGION_STYLES = {
    # S_0^ε 区域（红色、正方形）
    "S": RegionStyle(color="#ea9999", marker="s", label=r"$S_0^\varepsilon$"),
    # BS_0^ε 区域（蓝色、三角形）
    "BS": RegionStyle(color="#1f77b4", marker="o", label=r"$BS_0^\varepsilon$"),
    # RBS_0^ε 区域（绿色、圆形）
    "RBS": RegionStyle(color="#4daf4a", marker="^", label=r"$RBS_0^\varepsilon$"),
    # T_eqm 区域（紫色、菱形）
    "T_eqm": RegionStyle(color="#984ea3", marker="D", label=r"$T_{eqm}$"),
} 