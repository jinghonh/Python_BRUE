% 交通网络分析主脚本
% 用于运行不同参数组合的交通网络分析

% 清理工作空间
clear;
clc;
close all;

% 设置参数
zeta = 31;          % 可选值：15 或 31
subset_index = 2;   % zeta=15时可选0,1；zeta=31时可选0,1,2

% 验证输入参数
if ~ismember(zeta, [15, 31])
    error('zeta必须为15或31');
end

if zeta == 15 && ~ismember(subset_index, [0, 1])
    error('当zeta=15时，subset_index必须为0或1');
end

if zeta == 31 && ~ismember(subset_index, [0, 1, 2])
    error('当zeta=31时，subset_index必须为0、1或2');
end

% 根据zeta和subset_index选择搜索范围
if zeta == 15
    switch subset_index
        case 0  % zeta 15 0
            rangeMin = [3201, 3151, 1201];
            rangeMax = [5002, 4800, 2400];
        case 1  % zeta 15 1
            rangeMin = [4200, 3000, 0, 0, 0];
            rangeMax = [4900, 4200, 1610, 1610, 2000];
    end
else  % zeta == 31
    switch subset_index
        case 0  % zeta 31 0
            rangeMin = [1588, 2800, 1198];
            rangeMax = [4903, 5823, 2840];
        case 1  % zeta 31 1
            rangeMin = [500, 1200, 0, 0, 0];
            rangeMax = [4899, 6001, 2507, 2503, 2717];
        case 2  % zeta 31 2
            rangeMin = [4000, 3000, 0, 0, 902, 0];
            rangeMax = [5007, 4201, 3909, 3909, 1353, 1601];
    end
end

% 显示当前参数设置
fprintf('当前参数设置：\n');
fprintf('zeta = %d\n', zeta);
fprintf('subset_index = %d\n', subset_index);
fprintf('rangeMin = [%s]\n', num2str(rangeMin));
fprintf('rangeMax = [%s]\n', num2str(rangeMax));

% 运行分析
analyzeTrafficNetwork(zeta, rangeMin, rangeMax, subset_index); 