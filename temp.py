import numpy as np
import matplotlib.pyplot as plt

# 定义无差异曲线，例如： m = 10 / (t - 1) + 2
t = np.linspace(1.5, 6, 200)
m = 10 / (t - 1) + 2

# 时间容忍线
t_star = 2  # 最短时间
epsilon = 2
t_tol = t_star + epsilon

# 绘图
plt.figure(figsize=(8, 6))
plt.plot(t, m, label='Indifference Curve $\mathcal{U}_w$', color='blue')
plt.axvline(x=t_tol, color='red', linestyle='--', label=r'Time tolerance $t^* + \epsilon$')

# 填充可接受区域
tt, mm = np.meshgrid(np.linspace(1.5, t_tol, 300), np.linspace(0, 20, 300))
accept_region = mm < 10 / (tt - 1) + 2
plt.contourf(tt, mm, accept_region, levels=[0.5, 1], colors=['lightgreen'], alpha=0.5)

plt.xlabel('Travel Time $t$')
plt.ylabel('Monetary Cost $m$')
plt.title('Acceptable Cost Region $A_w$')
plt.legend()
plt.grid(True)
plt.xlim(1.5, 6)
plt.ylim(0, 20)
plt.show()