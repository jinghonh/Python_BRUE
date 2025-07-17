import random

# ---------- 1. 数据准备 ----------
# 12 条 Link，编号 1~12
LINK_CNT = 12

# 10 条 Path（按题目顺序）
PATH1 = [
    [1, 7],          # OD1-Path1
    [1, 8, 2],       # OD1-Path2
    [2, 9, 12],      # OD1-Path3
    [3, 11],         # OD1-Path4
    [3, 10, 12]]     # OD1-Path5
PATH2=[
    [6, 11],         # OD2-Path1
    [6, 10, 12],     # OD2-Path2
    [5, 9, 12],      # OD2-Path3
    [4, 7],          # OD2-Path4
    [4, 8, 12],      # OD2-Path5
]


def is_pareto_nondominated(pairs):
    """判断一组 (PTTF, PM) 是否两两非支配。默认均为最小化目标。"""
    n = len(pairs)
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            a1, b1 = pairs[i]
            a2, b2 = pairs[j]
            # i 支配 j 的判定：两目标均 <= 且至少一个 <
            if a1 <= a2 and b1 <= b2:
                return False
    return True


# ---------- 2. 随机搜索 ----------
MAX_VAL = 5          # FFT/M 最大值，可按需要调小
RND_TRIALS = 200_000  # 搜索轮数，视机器性能可调

for _ in range(RND_TRIALS):
    # 反相关启发式：FFT 随机取正整数，M 取 “MAX_VAL+1-FFT”
    FFT = [random.randint(1, MAX_VAL) for _ in range(LINK_CNT)]
    M   = [MAX_VAL + 1 - f           for f in FFT]   # 反向映射，保证正整数

    # 计算 10 条 Path 的 PTTF、PM
    pair1 = []
    pair2 = []
    for path in PATH1:
        pttf = sum(FFT[i - 1] for i in path)
        pm   = sum(M[i - 1]   for i in path)
        pair1.append((pttf, pm))

    for path in PATH2:
        pttf = sum(FFT[i - 1] for i in path)
        pm   = sum(M[i - 1]   for i in path)
        pair2.append((pttf, pm))

    if is_pareto_nondominated(pair1) and is_pareto_nondominated(pair2):
        break
else:
    raise RuntimeError("在设定范围内未找到可行解，请增大 MAX_VAL 或提高迭代次数。")


# ---------- 3. 结果输出 ----------
print("★ 找到一组互不支配解：\n")
print("FFT =", FFT)
print("M   =", M, "\n")

print("对应 10 条 Path 的 (PTTF, PM)：")
for idx, (pttf, pm) in enumerate(pair1, 1):
    print(f"Path{idx:<2}: ({pttf:>3}, {pm:>3})")

for idx, (pttf, pm) in enumerate(pair2, 1):
    print(f"Path{idx:<2}: ({pttf:>3}, {pm:>3})")