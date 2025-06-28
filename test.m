%% ---------------------- 参数区 ----------------------
e      = 35;      % 允许的路径差上限
rho    = 15;
sigma  = 0.02;
m1     = 20; 
m2     = 15;
m5     =  2;

% f1、f2 取值范围与栅格分辨率
N      = 800;                               % 采样点数 (>=200 推荐)
f1_vec = linspace(1000, 8000, N);           % f1 ∈ [1000,8000]
f2_vec = linspace(   0, 8000, N);           % f2 ∈ [   0,8000]
[F1, F2] = meshgrid(f1_vec, f2_vec);        % 网格

%% ---------------------- 公式区 ----------------------
X1 = F1;
X2 = F2;
X3 = 10000 - F1 - F2;    % 同时也是 X5 与 X8
X5 = X3;   X8 = X3;

link1 = 18  .* (1 + 0.15 * (X1/3600).^4);
link2 = 22.5.* (1 + 0.15 * (X2/3600).^4);
link3 = 12  .* (1 + 0.15 * (X3/1800).^4);
link5 = 2.4 .* (1 + 0.15 * (X5/1800).^4);
link8 = 12  .* (1 + 0.15 * (X8/1800).^4);

path1 = link1 + rho .* (1 - exp(-sigma .* link1));
path2 = link2 + rho .* (1 - exp(-sigma .* link2));
path5 = link3 + link5 + link8 + ...
        rho .* (1 - exp(-sigma .* (link3 + link5 + link8)));

%% ---------------------- 区域判定 ----------------------
C_eq1 = abs(path1 - path2) <= e;
C_eq2 = abs(path1 - path5) <= e;
C_eq3 = abs(path2 - path5) <= e;

C_sign1 = (path1 - path2) * (m1 - m2) < 0;   % 等价于 path1 - path2 < 0
C_sign2 = (path1 - path5) * (m1 - m5) < 0;   % 等价于 path1 - path5 < 0
C_sign3 = (path2 - path5) * (m2 - m5) < 0;   % 等价于 path2 - path5 < 0

C_nonneg = (10000 - F1 - F2) >= 0;

region_p2 =  C_eq1 & C_eq2 & C_eq3;                % "大" 可行集
region_p1 =  region_p2 & C_sign1 & C_sign2 & C_sign3 & C_nonneg;  % "小" 可行集

%% ---------------------- 绘图 ----------------------
figure('Color','w');
hold on;

% 先画 p2（蓝色），再叠加 p1（红色）
% 采用散点而非 imagesc，可避免 NaN 处理并保留坐标比例
scatter(F1(region_p2), F2(region_p2), 6, [0.25 0.55 1.00], 'filled'); % blue
scatter(F1(region_p1), F2(region_p1), 6, [1.00 0.10 0.10], 'filled'); % red

xlabel('f_1','FontSize',13);
ylabel('f_2','FontSize',13);
xlim([1000 8000]); ylim([0 8000]);
axis square; box on;
set(gca,'FontSize',12,'LineWidth',1);
legend({'p2 区域','p1 区域'},'Location','best');
title('RegionPlot 对应 Matlab 结果');