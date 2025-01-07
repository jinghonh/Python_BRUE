clear
clc
close all
%%  0
orderedPair = [1,1;2,2;3,3;3,5;3,8];
%%  1
% orderedPair=[1,1;2,2;3,3;3,7;4,4;4,8;5,3;5,5;5,8];
%%  2
% orderedPair=[1,1;2,2;3,3;5,3;4,4;6,4;5,5;6,6;3,7;6,7;4,8;5,8];
%%
relationMatrix = pairsToMatrix(orderedPair);
total = 0;
%%
n = size(relationMatrix,1);
tic
%% 修正
bound = 0;
% zeta = 15 0
% rangeMin = [3201, 3151, 1201];  % 每个维度的最小值
% rangeMax = [5002, 4800, 2400]; % 每个维度的最大值
%zeta 15 1
% rangeMin = [4200, 3000, 0, 0, 0];  % 每个维度的最小值
% rangeMax = [4900, 4200, 1610, 1610, 2000]; % 每个维度的最大值
% zeta 31 0
% rangeMin = [1588, 2800, 1198];  % 每个维度的最小值
% rangeMax = [4903, 5823, 2840]; % 每个维度的最大值
% zeta 31 1 _1
rangeMin = [500, 1200, 0, 0, 0];  % 每个维度的最小值
rangeMax = [4899, 6001, 2507, 2503, 2717]; % 每个维度的最大值
% zeta 31 2
% rangeMin = [4000, 3000, 0, 0, 902, 0];  % 每个维度的最小值
% rangeMax = [5007, 4201, 3909, 3909, 1353, 1601]; % 每个维度的最大值
for ii = 10:65
    dimNUm = ones(1,n)*ii;
    % dimNUm(4:5) = [1,1];
    samples = cell(1, n); % 使用 cell 数组来存储每个维度的样本
    rangeMin=rangeMin-bound;
    rangeMin(rangeMin<0)=0;
    rangeMax=rangeMax+bound;

    for i = 1:n
        samples{i} = linspace(rangeMin(i), rangeMax(i), dimNUm(i));
    end

    samplesMat = samples{1}';
    
    for i = 2:n-1
        mat = samples{i};
        lenMat = size(samplesMat, 1);
        lenNewCol = length(mat);
        samplesMat = [repmat(samplesMat, lenNewCol, 1), kron(mat', ones(lenMat, 1))];
        samplesMat(sum(samplesMat, 2) > 10000, :) = [];
    end
    % 计算最后一个维度值并附加到矩阵中
    lastDim = 10000 - sum(samplesMat, 2);
    samplesMat = [samplesMat, lastDim];

    % 计算目标函数和惩罚量
    [ff, sum_err] = FevalS1(samplesMat,relationMatrix);
    toc
    % sum_err = sum(err,2);
    valid = find(sum_err==0);
    disp(min(samplesMat(valid,:)))
    disp(max(samplesMat(valid,:)))
    temp_ff = ff(valid,:);
    if total==0
        total = temp_ff;
    else
        total = [total;temp_ff];
    end
    disp(ii)
    if ~isempty(valid)
        rangeMin = min(samplesMat(valid,:));
        rangeMax = max(samplesMat(valid,:));
        bound = (rangeMax-rangeMin)/ii;
        bound(bound<20)=20;
    else
        bound = 0;
    end
end
s = scatter(total (:,1),total(:,2),'.');
s.SizeData = 5;


%%
function [ff, err] = FevalS1(f,M)
    zeta = 31;
    n = size(M,1);
    % 定义常数
    money = [20, 15, 1, 0, 0, 0, 0, 1];
    freeFlowTime = [18,22.5,12,24,2.4,6,24,12];
    maxCapacity = [3600,3600,1800,1800,1800,1800,1800,1800];
    money = money*M';
    disp("block0")
    toc
    err = zeros(size(f, 1), 1);
        % 计算总值限制惩罚
    % err(:) = err(:)+PD(sum(f, 2) - 10001);
    % err(:) = err(:)+PD(9999 - sum(f, 2));
    f(err>0,:) = [];
    err(err>0) = [];

    disp("block1")
    toc
    x = f*M;
    disp("block2")
    toc
    % 计算 RT 值
    RT = realTime(x,M,freeFlowTime,maxCapacity);
    disp("block3")
    toc

    % 获取两两列的组合索引
    combs = nchoosek(1:n, 2);
    % 根据组合索引进行列的两两相减
    RTDiffs = abs(RT(:, combs(:,1)) - RT(:, combs(:,2)));
    disp("block4")
    toc
    % 计算路径差异的惩罚
    
    % valid_index1 = err==0;
    for i = 1:size(RTDiffs,2)
        % valid_index = find(err(valid_index1)==0);
        err(:) = err(:)+PD(RTDiffs(:, i) - zeta);
    end
    disp("block5")
    toc
    % 计算基于 money 差异的惩罚
    for i = 1:n
        for j = i+1:n
            err(:) = err(:)+PD((RT(:,i)- RT(:,j)) .* (money(i) - money(j)));
        end
    end
    disp("block6")
    toc
    % 计算目标函数
    ff(:, 1) = sum(f.*RT,2);
    ff(:, 2) = f * money';
end



function t = ta(x, t0, C)
    % 交通时间计算公式
    t = t0 .* (1 + 0.15 .* (x ./ C).^4);
end

function p = PD(x)
    % 惩罚函数，如果约束不满足则引入惩罚
    p = max(0, x);
end

function M = pairsToMatrix(pairs)
    % 将二元组转换为矩阵
    % 输入:
    %   pairs: 一个 n x 2 的矩阵，其中每一行是一个二元组
    % 输出:
    %   M: 一个矩阵，包含二元组的值

    % 初始化矩阵
    M = zeros(max(pairs));
    
    for i  = 1:size(pairs,1)
        M(pairs(i,1),pairs(i,2)) = 1;
    end

end

function RT = realTime(x,M,freeFlowTime,maxCapacity)
    index = find(sum(M)~=0);
    time = ta(x(:,index),freeFlowTime(index),maxCapacity(index));
    pathTime = time* M(:,index)';
    RT = pathTime + 15 * (1 - exp(-0.02 *pathTime));
end