import numpy as np
import matplotlib.pyplot as plt

# 中文乱码，mac 使用宋体
plt.rcParams['font.sans-serif'] = ['Songti SC']
plt.rcParams['axes.unicode_minus'] = False


# ---------- 参数 ----------
n = 4
num_samples = 800  # 采样点数量

# # ---------- 1) 随机采样位于 Σxi = 1 的 4D 点 ----------
# points = np.random.dirichlet(alpha=np.ones(n), size=num_samples)  # shape (num_samples, 4)

# 平均采样 sum(x)=1 的点，不要随机采样 使用网格采样
# 创建三维网格，第四维由约束条件确定
grid_size = int(np.cbrt(num_samples))  # 立方根，确保总点数接近 num_samples
x = np.linspace(0, 1, grid_size)
y = np.linspace(0, 1, grid_size)
z = np.linspace(0, 1, grid_size)

# 创建网格点
points = []
for xi in x:
    for yi in y:
        for zi in z:
            # 确保 xi + yi + zi <= 1，这样 wi 才能为正
            if xi + yi + zi <= 1:
                wi = 1 - xi - yi - zi
                points.append([xi, yi, zi, wi])

# 转换为 numpy 数组
points = np.array(points)

# 如果点数不足，可以通过插值或其他方法增加点数
if len(points) < num_samples:
    # 简单复制点以达到所需数量
    indices = np.random.choice(len(points), num_samples - len(points))
    extra_points = points[indices]
    points = np.vstack([points, extra_points])
# 如果点数过多，随机选择所需数量的点
elif len(points) > num_samples:
    indices = np.random.choice(len(points), num_samples, replace=False)
    points = points[indices]

# 验证 sum(x)=1
print(np.sum(points, axis=1))


# ---------- 2) 超平面中心 ----------
p0 = np.ones(n) / n  # (0.25,0.25,0.25,0.25)

# ---------- 3) 构造三个与 (1,1,1,1)⊤ 正交且互正交的基向量 ----------
N = np.array([[1, -1,  0,  0],
              [1,  1, -2,  0],
              [1,  1,  1, -3]], dtype=float).T   # 4×3
Q, _ = np.linalg.qr(N)     # 正交化 → Q: 4×3
B = Q                      # 每列单位且与 (1,1,1,1) 正交

# ---------- 4) 3D 可逆投影 ----------
coords3 = (points - p0) @ B        # num_samples × 3
points_rec3 = p0 + coords3 @ B.T   # 逆映射
max_err3 = np.max(np.abs(points - points_rec3))

# ---------- 5) 2D 不可逆投影 ----------
basis2 = B[:, :2]                  # 只取前两列 → 4×2
coords2 = (points - p0) @ basis2   # num_samples × 2
points_rec2 = p0 + coords2 @ basis2.T   # 只能恢复到 4D 子空间
max_err2 = np.max(np.abs(points - points_rec2))

print("Max reconstruction error (3D → 4D):", max_err3)
print("Max reconstruction error (2D → 4D):", max_err2)

# ---------- 6) 绘制 3D 散点 ----------
fig1 = plt.figure()
ax1 = fig1.add_subplot(111, projection='3d')
ax1.scatter(coords3[:, 0], coords3[:, 1], coords3[:, 2], s=12, alpha=0.7)
ax1.set_xlabel('c₁')
ax1.set_ylabel('c₂')
ax1.set_zlabel('c₃')
ax1.set_title('可逆：4D 超平面 → 3D')

# ---------- 7) 绘制 2D 散点 ----------
fig2 = plt.figure()
plt.scatter(coords2[:, 0], coords2[:, 1], s=12, alpha=0.7)
plt.xlabel('c₁')
plt.ylabel('c₂')
plt.title('不可逆：4D 超平面 → 2D')
plt.gca().set_aspect('equal', adjustable='box')

plt.show()