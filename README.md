# BRUE (Boundedly Rational User Equilibrium) 交通网络分析系统

## 项目简介

本项目实现了有界理性用户均衡(BRUE)模型，用于分析交通网络中的路径选择和均衡状态。该系统考虑了用户的有界理性决策行为，通过优化模型求解交通流分布，并提供可视化分析功能。

## 主要特性

- **多种求解器支持**：包括基础求解器、集合求解器和迭代求解器
- **多种网络配置**：支持单一起终点对、多起终点对和复杂网络配置
- **双目标优化**：同时考虑时间成本和金钱成本的帕累托分析
- **可视化分析**：提供路径成本分析图表和迭代过程可视化
- **有界理性建模**：考虑用户决策的不完全理性特征
- **迭代求解算法**：逐步识别有效路径的迭代优化方法

## 项目结构

```
Python_BRUE/
├── src/                    # 源代码目录
│   ├── brue_base.py        # 基础抽象类
│   ├── brue_solver.py      # 主要求解器实现
│   ├── brue_set_solver.py  # 集合求解器（帕累托分析）
│   ├── brue_matlab_solver.py # MATLAB接口求解器
│   └── traffic_network_config.py # 网络配置类
├── matlab/                 # MATLAB实现
│   ├── analyzeTrafficNetwork.m  # 网络分析脚本
│   ├── plotPathCosts.m     # 成本绘图脚本
│   ├── main.m              # 主程序
│   └── test_run.m          # 测试运行脚本
├── examples/               # 示例代码
│   └── test.m              # 测试示例
├── output/                 # 输出目录
│   └── results/            # 计算结果
├── cache/                  # 缓存目录
├── pyproject.toml         # 项目配置文件
├── uv.lock                # 依赖锁文件
└── README.md              # 项目说明文档
```

## 依赖环境

### Python 依赖
- Python >= 3.12
- pyomo >= 6.9.2 (优化建模)
- matplotlib >= 3.10.3 (可视化)
- rich >= 14.0.0 (终端美化输出)
- tqdm >= 4.67.1 (进度条)

### 求解器要求
- IPOPT (推荐，用于非线性优化)
- 或其他Pyomo支持的求解器

## 安装方式

### 使用 uv (推荐)
```bash
# 克隆项目
git clone <repository-url>
cd Python_BRUE

# 使用uv安装依赖
uv sync
```

### 使用 pip
```bash
# 克隆项目
git clone <repository-url>
cd Python_BRUE

# 安装依赖
pip install -e .
```

## 核心模型

### BRUE模型数学表述

**目标函数**：
```
minimize ε (有界理性误差)
```

**约束条件**：
1. **流量平衡约束**：所有路径流量之和等于总需求
2. **有界理性约束**：路径选择在感知成本±ε范围内
3. **成本约束**：路径成本由BPR函数计算
4. **非负约束**：所有变量非负

### BPR (Bureau of Public Roads) 函数
```
t_a = t_0 * (1 + α * (v_a/C_a)^β)
```
其中：
- t_a：路段实际旅行时间
- t_0：自由流旅行时间
- v_a：路段流量
- C_a：路段通行能力
- α, β：BPR参数（默认 α=0.15, β=4）

## 使用方法

### 基础使用示例

```python
from src.brue_solver import BRUESolver
from src.traffic_network_config import TrafficNetworkConfig

# 创建网络配置
config = TrafficNetworkConfig.create_basic_network()

# 创建求解器
solver = BRUESolver(config)

# 运行迭代求解
results = solver.run_with_iterations()

# 绘制初始成本分析
solver.plot_initial_costs()
```

### 多起终点对网络分析

```python
# 创建多OD网络配置
multi_config = TrafficNetworkConfig.create_multi_od_network()
multi_solver = BRUESolver(multi_config)

# 运行分析
results = multi_solver.run_with_iterations()
```

### 帕累托分析

```python
from src.brue_set_solver import BRUESetSolver

# 使用集合求解器进行帕累托分析
set_solver = BRUESetSolver(config)
pareto_results = set_solver.solve_pareto_optimal_set()
```

## 网络配置

项目提供四种预定义网络配置：

### 1. 基础网络 (create_basic_network)
- 8条路径，6个起终点对
- 单一OD组，总需求10000
- 适用于基础分析

### 2. 单一起终点网络 (create_single_od_network)
- 简化网络，用于算法验证
- 单一起终点对

### 3. 多起终点网络 (create_multi_od_network)
- 复杂网络，19条路径，14个起终点对
- 多个OD组，不同权重
- 适用于复杂场景分析

### 4. 反例网络 (create_anti_example_network)
- 特殊设计的网络，用于验证算法鲁棒性
- 10条路径，2个OD组
- 测试边界情况

### 自定义网络配置

```python
from src.traffic_network_config import TrafficNetworkConfig

# 创建自定义配置
custom_config = TrafficNetworkConfig(
   num_links=6,
   num_paths=4,
   total_demand=5000,
   od_groups={'Group1': [1, 2], 'Group2': [3, 4]},
   od_demands={'Group1': 3000, 'Group2': 2000},
   free_flow_time={1: 15, 2: 20, 3: 10, 4: 25, 5: 5, 6: 8},
   link_money_cost={1: 10, 2: 15, 3: 0, 4: 5, 5: 0, 6: 2},
   link_capacity={1: 2000, 2: 2500, 3: 1500, 4: 1800, 5: 1200, 6: 1600},
   path_link_matrix={
      # (路径, 路段): 关联关系 (0或1)
      (1, 1): 1, (1, 2): 0,  # 路径1使用路段1
      # ... 其他关联关系
   }
)
```

## 输出结果

### 终端输出
- 迭代过程信息
- 路径成本分析表格
- 有效路径识别结果
- 帕累托前沿分析

### 可视化图表
1. **路径成本散点图**：时间成本vs金钱成本
2. **有效路径识别**：标识帕累托最优路径
3. **迭代过程分析**：显示算法收敛过程
4. **初始状态分析**：不同负载情况下的路径比较

## 核心算法

### 1. 迭代求解算法
```
1. 初始化：从单一路径开始
2. 求解BRUE模型
3. 识别有效路径（在ε范围内的路径）
4. 添加新有效路径到约束集合
5. 重复步骤2-4直到收敛
```

### 2. 帕累托分析
- 识别时间-金钱成本空间中的非支配路径
- 计算帕累托前沿
- 分析路径选择的权衡关系

### 3. 有界理性建模
- 引入感知误差ε
- 允许次优路径的存在
- 反映真实用户决策行为

## MATLAB接口
项目包含MATLAB接口文件，支持：
- `analyzeTrafficNetwork.m`：网络分析主函数
- `plotPathCosts.m`：成本可视化
- `main.m`：主要分析脚本

## 缓存机制
项目实现了结果缓存机制，缓存文件存储在`cache/`目录中

## 测试和算法验证

### 测试文件说明
- **test_paper_algorithm.py**: 实现严格按照论文的多OD对BRUE算法
- **test_corrected.py**: 目前为空，预留用于修正测试
- **demo.py**: 演示帕累托非支配路径的随机搜索算法
- **temp.py**: 临时测试文件，用于4维超平面的可视化

### 算法验证
项目包含多种算法验证方法：
1. 论文算法的严格实现
2. 帕累托最优解的计算
3. 网格搜索验证
4. 边界情况测试

## 扩展功能

### MATLAB接口
项目包含MATLAB接口文件，支持：
- `analyzeTrafficNetwork.m`：网络分析主函数
- `plotPathCosts.m`：成本可视化
- `main.m`：主要分析脚本

### 高级分析
- 敏感性分析
- 不同权重下的均衡状态
- 网络容量影响分析

### 缓存机制
项目实现了结果缓存机制，缓存文件存储在`cache/`目录中：
- `cache_*.mat`：MATLAB格式的缓存文件
- 自动缓存计算结果，提高重复运行效率

## 故障排除

### 常见问题

1. **求解器错误**
   ```bash
   # 安装IPOPT求解器
   conda install -c conda-forge ipopt
   ```

2. **依赖版本冲突**
   ```bash
   # 使用uv解决依赖
   uv sync --refresh
   ```

3. **可视化显示问题**
   ```python
   # 设置matplotlib后端
   import matplotlib
   matplotlib.use('TkAgg')  # 或 'Qt5Agg'
   ```

## 开发指南

### 添加新的网络配置
1. 在`TrafficNetworkConfig`类中添加新的类方法
2. 定义网络参数：路径数、容量、成本等
3. 设置路径-路段关联矩阵

### 扩展求解算法
1. 继承`BRUEBase`抽象类
2. 实现必要的抽象方法
3. 添加特定的约束和目标函数

### 自定义可视化
1. 修改`plot_cost_analysis`方法
2. 添加新的图表类型
3. 调整颜色方案和标签

## 学术引用

如果您在学术研究中使用本项目，请引用相关的BRUE模型文献。

## 许可证

本项目采用MIT许可证，详见LICENSE文件。

## 联系方式

如有问题或建议，请通过GitHub Issues提交反馈。

---

*最后更新：2025年7月14日*

# BRUE项目使用指南

## 基本使用流程

1. 创建交通网络配置
2. 初始化求解器（有多种选择：BRUE, UE, δ-EBR-MUE）
3. 运行迭代求解过程
4. 分析路径成本和流量分配结果

## 网络配置

[traffic_network_config.py](src/traffic_network_config.py)提供了多种预设网络配置：

- `create_basic_network()`: 基础测试网络
- `create_single_od_network()`: 单一OD对网络
- `create_multi_od_network()`: 多OD对网络
- `create_two_od_network()`: 两个OD对网络
- `create_large_single_od_network()`: 大型单OD对网络，具有15条路径

## 求解器选项

### BRUE求解器

```python
from src import BRUESolver
from src import TrafficNetworkConfig

config = TrafficNetworkConfig.create_multi_od_network()
solver = BRUESolver(config)
results = solver.run_with_iterations()
```

### 用户均衡(UE)求解器

```python
from src import UESolver
from src import TrafficNetworkConfig

config = TrafficNetworkConfig.create_multi_od_network()
solver = UESolver(config)
results = solver.run()
```

### δ-EBR-MUE求解器

这是最新添加的求解器，实现了δ-有界ε-非严格支配的多目标用户均衡模型。

```python
from src import DeltaEBRMUESolver
from src import TrafficNetworkConfig

config = TrafficNetworkConfig.create_multi_od_network()
delta = 2.0  # 时间容忍度
epsilon = (0.5, 0.5)  # 支配容忍度向量(时间, 金钱)

# 使用凸优化方法
solver = DeltaEBRMUESolver(config, delta=delta, epsilon=epsilon)
result = solver.run(use_epsilon=True, solver_method='convex')

# 使用变分不等式(VI)方法
vi_result = solver.run(use_epsilon=True, solver_method='vi')
```

## δ-EBR-MUE求解方法比较

δ-EBR-MUE求解器实现了两种子问题求解方法：

1. **凸优化方法**：基于Beckmann公式，通过求解等价的凸优化问题获取均衡流量
2. **变分不等式(VI)方法**：直接基于VI理论，使用投影梯度法求解变分不等式问题

测试结果显示，两种方法在性能和结果上存在差异：

| 对比项 | 凸优化方法 | VI方法 |
| ------ | ---------- | ------ |
| 求解时间 | 更快 | 约为凸优化方法的4倍 |
| 收敛性 | 更好，通常能快速收敛 | 收敛性较差，可能需要更多迭代 |
| 路径选择 | 倾向于选择更少的路径 | 倾向于分散流量到更多路径 |
| 时间成本 | 通常更低 | 通常更高 |
| 金钱成本 | 通常更高 | 通常更低 |

两种方法的主要理论差异：

- **凸优化方法**：适用于对称问题（雅可比矩阵对称），计算效率高
- **VI方法**：适用范围更广，可处理非对称/非单调问题，但计算效率较低

**使用建议**：
- 对于标准交通分配问题，凸优化方法效率更高、收敛性更好
- 对于复杂非对称问题或需要考虑多目标权衡的情况，VI方法可能提供更多样化的解

## 结果分析

求解器提供多种结果分析和可视化方法：

- `display_results()`: 显示计算结果表格
- `analyze_path_costs()`: 分析路径成本和有效路径
- `plot_cost_analysis()`: 可视化路径成本分析
- `plot_initial_costs()`: 比较不同场景下的路径成本
