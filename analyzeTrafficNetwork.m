function analyzeTrafficNetwork(zeta, rangeMin, rangeMax, subset_index)
    % 交通网络分析主函数
    % 实现交通流量分配和路径时间计算
    %
    % 输入参数:
    %   zeta         - 路径时间差异约束值
    %   rangeMin     - 每个维度的最小值数组 [d1_min, d2_min, ...]
    %   rangeMax     - 每个维度的最大值数组 [d1_max, d2_max, ...]
    %   subset_index - 子集索引，用于选择不同的路径配对
    
    % 检查是否存在缓存文件
    cacheFileName = sprintf('cache/cache_zeta%d_subset%d.mat', zeta, subset_index);
    pathConstraintCacheFileName = sprintf('cache/cache_path_only_zeta%d_subset%d.mat', zeta, subset_index);
    
    if exist(cacheFileName, 'file') && exist(pathConstraintCacheFileName, 'file')
        fprintf('加载缓存数据...\n');
        load(cacheFileName, 'totalValidCost', 'totalValidFlow', 'relationMatrix');
        load(pathConstraintCacheFileName, 'totalPathValidCost', 'totalPathValidFlow');
        fprintf('加载缓存数据完成...\n');
        toc;
        
        % 直接绘图
        q = 10000;
        if ~isempty(totalValidFlow)
            selectedIndices = randperm(size(totalValidFlow, 1), min(q, size(totalValidFlow, 1)));
            % plotResults(totalValidCost, selectedIndices);
            plotPathCosts(totalValidFlow, relationMatrix, selectedIndices);
            % 绘制流量向量可视化
            plotFlowVectors(totalValidFlow, totalPathValidFlow, relationMatrix, selectedIndices);
            % 绘制三维流量可视化
            % plot3DFlowVectors(totalValidFlow, totalPathValidFlow, relationMatrix, selectedIndices, subset_index);
        end
        return;
    end
    
    % 参数验证
    validateInputs(zeta, rangeMin, rangeMax);
    
    %% 根据subset_index选择orderedPair
    orderedPair = selectOrderedPair(zeta, subset_index);
    relationMatrix = pairsToMatrix(orderedPair);
    
    %% 初始化参数
    n = size(relationMatrix, 1);
    totalValidCost = [];
    totalValidFlow = [];
    totalPathValidCost = []; % 只满足路径约束的成本
    totalPathValidFlow = []; % 只满足路径约束的流量
    bound = 0;
    
    %% 确保rangeMin和rangeMax的长度与n一致
    if length(rangeMin) < n
        rangeMin(end+1:n) = rangeMin(end);
    elseif length(rangeMin) > n
        rangeMin = rangeMin(1:n);
    end
    
    if length(rangeMax) < n
        rangeMax(end+1:n) = rangeMax(end);
    elseif length(rangeMax) > n
        rangeMax = rangeMax(1:n);
    end
    
    %% 主循环
    tic
    % 添加进度显示
    h = waitbar(0, '开始计算...');
    minIt = 30;
    maxIt = 40;
    
    % 设置数据压缩参数
    maxPointsPerIteration = 1e5; % 每次迭代最大保留点数

    % 传统网格采样策略
    for ii = minIt:maxIt
        waitbar((ii-minIt)/(maxIt-minIt), h, sprintf('正在计算第 %d/%d 次迭代...', ii-minIt+1, maxIt-minIt+1));
        [samplesMat, totalValidCost, totalValidFlow, totalPathValidCost, totalPathValidFlow] = processIteration(ii, n, rangeMin, rangeMax, bound, relationMatrix, totalValidCost, totalValidFlow, totalPathValidCost, totalPathValidFlow, zeta);
        
        % 增量数据压缩：当有效点数量超过阈值时进行压缩
        if size(totalValidCost, 1) > maxPointsPerIteration
            fprintf('迭代 %d: 压缩全部约束数据点，从 %d 个点压缩...\n', ii, size(totalValidCost, 1));
            [totalValidCost, totalValidFlow] = reduceDataPoints(totalValidCost, totalValidFlow, maxPointsPerIteration);
        end
        
        % 增量数据压缩：当路径约束点数量超过阈值时进行压缩
        if size(totalPathValidCost, 1) > maxPointsPerIteration
            fprintf('迭代 %d: 压缩路径约束数据点，从 %d 个点压缩...\n', ii, size(totalPathValidCost, 1));
            [totalPathValidCost, totalPathValidFlow] = reduceDataPoints(totalPathValidCost, totalPathValidFlow, maxPointsPerIteration);
        end
        
        % 更新搜索范围
        [~, path_err, ~] = evaluateObjective(samplesMat, relationMatrix, zeta);
        valid = find(path_err == 0);
        
        if ~isempty(valid)
            [rangeMin, rangeMax, bound] = updateSearchRange(samplesMat(valid,:), ii);
        else
            bound = 0;
        end
        
        fprintf('完成迭代 %d\n', ii);
    end
    close(h);
    toc
    
    % 保存满足所有约束的结果
    [totalValidCost,Ia,~]=unique(totalValidCost,"rows");
    totalValidFlow = totalValidFlow(Ia,:);
    
    % 使用网格采样减少数据点数量
    targetSize = 5e5; % 目标数据点数量：十万级别
    if size(totalValidCost, 1) > targetSize
        fprintf('正在进行全约束数据压缩，从 %d 个点压缩到目标 %d 个点左右...\n', size(totalValidCost, 1), targetSize);
        [totalValidCost, totalValidFlow] = reduceDataPoints(totalValidCost, totalValidFlow, targetSize);
    end
    
    % 保存满足路径约束的结果
    [totalPathValidCost,Ia,~]=unique(totalPathValidCost,"rows");
    totalPathValidFlow = totalPathValidFlow(Ia,:);
    
    % 使用网格采样减少数据点数量
    if size(totalPathValidCost, 1) > targetSize
        fprintf('正在进行路径约束数据压缩，从 %d 个点压缩到目标 %d 个点左右...\n', size(totalPathValidCost, 1), targetSize);
        [totalPathValidCost, totalPathValidFlow] = reduceDataPoints(totalPathValidCost, totalPathValidFlow, targetSize);
    end
    
    % 保存结果到缓存文件
    save(cacheFileName, 'totalValidCost', 'totalValidFlow', 'relationMatrix');
    fprintf('全部约束结果已保存到缓存文件%s\n', cacheFileName);
    
    % 保存只满足路径约束的数据到缓存文件
    save(pathConstraintCacheFileName, 'totalPathValidCost', 'totalPathValidFlow', 'relationMatrix');
    fprintf('只满足路径约束的结果已保存到缓存文件%s\n', pathConstraintCacheFileName);
    
    % 设置随机选择的流量向量数量
    q = 30;  % 可以根据需要修改
    
    % 随机选择q个流量向量的索引
    if ~isempty(totalValidFlow)
        selectedIndices = randperm(size(totalValidFlow, 1), min(q, size(totalValidFlow, 1)));
        
        % 先绘制总体结果
        % plotResults(totalValidCost, selectedIndices);
        
        % 调用路径成本绘图函数
        plotPathCosts(totalValidFlow, relationMatrix, selectedIndices);
        
        % 绘制流量向量可视化
        plotFlowVectors(totalValidFlow, totalPathValidFlow, relationMatrix, selectedIndices);
        
        % 绘制三维流量可视化
        plot3DFlowVectors(totalValidFlow, totalPathValidFlow, relationMatrix, selectedIndices, subset_index);
    else
        warning('没有找到可行解');
    end
end

function orderedPair = selectOrderedPair(zeta, subset_index)
    % 根据zeta和subset_index选择orderedPair
    if zeta == 15
        switch subset_index
            case 0  % zeta 15 0
                orderedPair = [1,1;2,2;3,3;3,5;3,8];
            case 1  % zeta 15 1
                orderedPair=[1,1;2,2;3,3;3,7;4,4;4,8;5,3;5,5;5,8];
        end
    else  % zeta == 31
        switch subset_index
            case 0  % zeta 31 0
                orderedPair = [1,1;2,2;3,3;3,5;3,8];
            case 1  % zeta 31 1
                orderedPair=[1,1;2,2;3,3;3,7;4,4;4,8;5,3;5,5;5,8];
            case 2  % zeta 31 2
                orderedPair = [1,1;2,2;3,3;5,3;4,4;6,4;5,5;6,6;3,7;6,7;4,8;5,8];
        end
    end
end

function validateInputs(zeta, rangeMin, rangeMax)
    % 验证输入参数
    if ~isnumeric(zeta) || ~isscalar(zeta) || zeta <= 0
        error('zeta必须是一个正数');
    end
    
    if ~isnumeric(rangeMin) || ~isvector(rangeMin)
        error('rangeMin必须是一个数值向量');
    end
    
    if ~isnumeric(rangeMax) || ~isvector(rangeMax)
        error('rangeMax必须是一个数值向量');
    end
    
    if length(rangeMin) ~= length(rangeMax)
        error('rangeMin和rangeMax的长度必须相同');
    end
    
    if any(rangeMin >= rangeMax)
        error('rangeMin的值必须小于对应的rangeMax值');
    end
end

function [samplesMat, totalValidCost, totalValidFlow, totalPathValidCost, totalPathValidFlow] = processIteration(ii, n, rangeMin, rangeMax, bound, relationMatrix, totalValidCost, totalValidFlow, totalPathValidCost, totalPathValidFlow, zeta)
    % 处理单次迭代
    dimNum = ones(1,n)*ii;
    samples = generateSamples(n, rangeMin-bound, rangeMax+bound, dimNum);
    samplesMat = combineSamples(samples, n);
    
    % 计算目标函数和约束违反
    [ff, path_err, money_err] = evaluateObjective(samplesMat, relationMatrix, zeta);
    
    % 满足所有约束的流量向量
    valid = path_err == 0 & money_err == 0;
    
    % 只满足路径约束的流量向量
    path_valid = path_err == 0;
    
    % 更新满足所有约束的结果
    if any(valid)
        if isempty(totalValidCost)
            totalValidCost = ff(valid,:);
        else
            totalValidCost = [totalValidCost; ff(valid,:)];
        end
        if isempty(totalValidFlow)
            totalValidFlow = samplesMat(valid,:);
        else
            totalValidFlow = [totalValidFlow;samplesMat(valid,:)];
        end
    end
    
    % 更新只满足路径约束的结果
    if any(path_valid)
        if isempty(totalPathValidCost)
            totalPathValidCost = ff(path_valid,:);
        else
            totalPathValidCost = [totalPathValidCost; ff(path_valid,:)];
        end
        if isempty(totalPathValidFlow)
            totalPathValidFlow = samplesMat(path_valid,:);
        else
            totalPathValidFlow = [totalPathValidFlow;samplesMat(path_valid,:)];
        end
    end
end

function samples = generateSamples(n, rangeMin, rangeMax, dimNum)
    % 生成采样点
    samples = cell(1, n);
    
    % 确保rangeMin不小于0且长度与n一致
    rangeMin = max(rangeMin, 0);
    if length(rangeMin) < n
        rangeMin(end+1:n) = rangeMin(end);
    elseif length(rangeMin) > n
        rangeMin = rangeMin(1:n);
    end
    
    % 确保rangeMax长度与n一致
    if length(rangeMax) < n
        rangeMax(end+1:n) = rangeMax(end);
    elseif length(rangeMax) > n
        rangeMax = rangeMax(1:n);
    end
    
    for i = 1:n
        samples{i} = linspace(rangeMin(i), rangeMax(i), dimNum(i));
    end
end

function samplesMat = combineSamples(samples, n)
    % 组合采样点
    
    samplesMat = samples{1}';
    
    for i = 2:n-1
        mat = samples{i};
        lenMat = size(samplesMat, 1);
        lenNewCol = length(mat);
        
        % 分批处理以优化大规模计算
        if lenMat * lenNewCol > 1e6
            batchSize = min(1e6, lenMat);
            numBatches = ceil(lenMat / batchSize);
            tempResult = zeros(lenMat * lenNewCol, i);
            
            for b = 1:numBatches
                startIdx = (b-1)*batchSize + 1;
                endIdx = min(b*batchSize, lenMat);
                currBatch = samplesMat(startIdx:endIdx, :);
                
                tempBatch = [repmat(currBatch, lenNewCol, 1), kron(mat', ones(size(currBatch, 1), 1))];
                startIdxResult = (startIdx-1) * lenNewCol + 1;
                endIdxResult = startIdxResult + size(tempBatch, 1) - 1;
                tempResult(startIdxResult:endIdxResult, :) = tempBatch;
            end
            
            samplesMat = tempResult;
        else
            samplesMat = [repmat(samplesMat, lenNewCol, 1), kron(mat', ones(lenMat, 1))];
        end
        
        % 移除超过总量的样本
        samplesMat(sum(samplesMat, 2) > 10000, :) = [];
    end
    
    % 计算最后一个维度
    lastDim = 10000 - sum(samplesMat, 2);
    samplesMat = [samplesMat, lastDim];
end

function [ff, path_err, money_err] = evaluateObjective(f, M, zeta)
    % 评估目标函数和约束违反
    n = size(M,1);
    
    % 使用persistent变量缓存常量以加快计算
    persistent money_cached freeFlowTime_cached maxCapacity_cached
    
    if isempty(money_cached)
        money_cached = [20, 15, 1, 0, 0, 0, 0, 1];
        freeFlowTime_cached = [18,22.5,12,24,2.4,6,24,12];
        maxCapacity_cached = [3600,3600,1800,1800,1800,1800,1800,1800];
    end
    
    % 使用缓存的常量
    money = money_cached;
    freeFlowTime = freeFlowTime_cached;
    maxCapacity = maxCapacity_cached;
    
    % 计算货币成本向量
    money = money*M';
    
    % 初始化误差
    path_err = zeros(size(f, 1), 1);
    money_err = zeros(size(f, 1), 1);
    
    % 计算实际行驶时间
    x = f*M;
    RT = calculateRealTime(x, M, freeFlowTime, maxCapacity);
    
    % 检查路径时间差异约束
    path_err = checkPathConstraints(RT, zeta, n, path_err);
    
    % 检查货币约束
    money_err = checkMoneyConstraints(RT, money, n, money_err);
    
    % 计算目标函数
    ff = calculateObjectives(f, RT, money);
end

function RT = calculateRealTime(x, M, freeFlowTime, maxCapacity)
    % 计算实际行驶时间
    index = find(sum(M)~=0);
    time = calculateTravelTime(x(:,index), freeFlowTime(index), maxCapacity(index));
    pathTime = time * M(:,index)';
    RT = pathTime + 15 * (1 - exp(-0.02 * pathTime));
end

function t = calculateTravelTime(x, t0, C)
    % 计算行驶时间
    t = t0 .* (1 + 0.15 .* (x ./ C).^4);
end

function err = checkPathConstraints(RT, zeta, n, err)
    % 检查路径约束（矢量化实现）
    % 对所有路径对一次性计算 |RT_i - RT_j| - zeta，并累加违反量
    combs = nchoosek(1:n, 2);                  % 所有(i,j)组合，i<j
    % 计算差值并扣减 zeta
    RTDiffs = abs(RT(:, combs(:,1)) - RT(:, combs(:,2))) - zeta;  % (samples × pairs)
    % 仅累加违反部分 (>0)
    err = err + sum(max(RTDiffs, 0), 2);
end

function err = checkMoneyConstraints(RT, money, n, err)
    % 检查货币约束（矢量化实现）
    % 对所有路径对一次性计算 (RT_i - RT_j)*(money_i - money_j)
    combs = nchoosek(1:n, 2);                          % 所有(i,j)组合
    moneyDiffs = money(combs(:,1)) - money(combs(:,2));% 1 × pairs
    % RT 差 (samples × pairs)
    RTDiffs = RT(:, combs(:,1)) - RT(:, combs(:,2));
    % 违反量：只保留正值再累加
    err = err + sum(max(RTDiffs .* moneyDiffs, 0), 2);
end

function ff = calculateObjectives(f, RT, money)
    % 计算目标函数值
    ff = zeros(size(f,1), 2);
    ff(:,1) = sum(f.*RT, 2);
    ff(:,2) = f * money';
end

function [rangeMin, rangeMax, bound] = updateSearchRange(validSamples, ii)
    % 更新搜索范围
    rangeMin = min(validSamples);
    rangeMax = max(validSamples);
    bound = (rangeMax-rangeMin)/ii;
    bound(bound<20) = 20;
end

function plotResults(totalValidCost, selectedIndices)
    % 绘制结果并使用boundary函数添加科研风格的非凸边界
    fig = figure('Name', 'Traffic Network Analysis', 'NumberTitle', 'off', 'Position', [100, 100, 800, 600]);
    set(fig, 'Color', 'white');  % 白色背景
    set(gca, 'FontName', 'Arial', 'FontSize', 10, 'Box', 'on', 'LineWidth', 1.2);
    
    % 创建优雅的色彩方案
    baseColor = [0.2, 0.4, 0.8];  % 基础蓝色
    highlightColor = [0.8, 0.3, 0.3];  % 高亮红色
    boundaryColor = [0.4, 0.5, 0.8];  % 中等蓝色边界
    
    hold on;
    
        % 绘制所有点（浅蓝色，半透明）
    scatter(totalValidCost(:,1), totalValidCost(:,2), 10, baseColor, 'filled', ...
        'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', 'none');
    % 应用boundary函数提取非凸边界（仅使用sf=0.3）
    try
        cleanData = totalValidCost;  % 使用全部数据
        
        % 固定sf值为0.3
        sf = 1;
        
        % 提取边界
        idx = boundary(cleanData(:,1), cleanData(:,2), sf);
        
        % 绘制边界
        if ~isempty(idx) && length(idx) > 3
            % 计算边界面积
            boundaryArea = polyarea(cleanData(idx,1), cleanData(idx,2));
            
            % 创建边界闭合多边形
            boundaryX = cleanData(idx,1);
            boundaryY = cleanData(idx,2);
            
            % 添加渐变填充
            patch('XData', boundaryX, 'YData', boundaryY, ...
                'FaceColor', [0.9, 0.95, 1], ... % 非常浅的蓝色
                'EdgeColor', boundaryColor, ... % 中等蓝色边缘
                'LineWidth', 1.5, ...
                'FaceAlpha', 0.9);  % 半透明填充
                
            % % 在图形右上角添加边界信息
            % infoText = sprintf('边界参数: sf=%.1f\n边界面积: %.1f', sf, boundaryArea);
            % annotation('textbox', [0.65, 0.8, 0.3, 0.1], ...
            %     'String', infoText, 'FitBoxToText', 'on', ...
            %     'BackgroundColor', [1, 1, 1, 0.7], ... % 半透明白色背景
            %     'EdgeColor', boundaryColor, ...
            %     'LineWidth', 1.0, ...
            %     'FontName', 'Arial', 'FontSize', 9);
        end
    catch e
        fprintf('boundary方法错误 (sf=%.2f): %s\n', sf, e.message);
    end

    % 绘制选定点（深红色，突出显示）
    scatter(totalValidCost(selectedIndices,1), totalValidCost(selectedIndices,2), 20, highlightColor, 'filled', ...
        'MarkerEdgeColor', 'black', 'MarkerEdgeAlpha', 0.5);

    
    % 设置坐标轴标签和标题
    xlabel('Total Travel Time', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Total Cost', 'FontSize', 12, 'FontWeight', 'bold');
    % title('交通网络分析结果与非凸边界', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial');
    
    % 添加科学风格网格
    grid on;
    grid minor;
    set(gca, 'GridAlpha', 0.1, 'MinorGridAlpha', 0.05, 'Layer', 'top');
    
    % 坐标轴美化
    ax = gca;
    ax.XAxis.TickLength = [0.01, 0.01];
    ax.YAxis.TickLength = [0.01, 0.01];
    
    % Add legend 
    legend({'Feasible Solutions', 'Concave Boundary', 'Selected Solutions'}, ...
        'Location', 'best', ...
        'FontName', 'Arial', 'FontSize', 9, ...
        'EdgeColor', [0.7, 0.7, 0.7], ...
        'Box', 'on');
    
    % 自动调整坐标轴，留出边距
    axis tight;
    axisLimits = axis;
    axisRange = [axisLimits(2)-axisLimits(1), axisLimits(4)-axisLimits(3)];
    axis([axisLimits(1)-0.02*axisRange(1), axisLimits(2)+0.02*axisRange(1), ...
          axisLimits(3)-0.02*axisRange(2), axisLimits(4)+0.02*axisRange(2)]);
    
    hold off;
    
    % 保存图像为高分辨率PNG文件
    figFile = sprintf('results/concave_boundary_%s.png', datestr(now, 'yyyymmdd_HHMMSS'));
    print(figFile, '-dpng', '-r300');
end

function plotFlowVectors(totalValidFlow, totalPathValidFlow, relationMatrix, selectedIndices)
    % Draw 2D projection visualization of flow vectors
    % 
    % Input parameters:
    %   totalValidFlow     - Flow vectors meeting all constraints
    %   totalPathValidFlow - Flow vectors meeting only path constraints
    %   relationMatrix     - Relation matrix
    %   selectedIndices    - Selected flow vector indices
    
    % Create figure with scientific style
    fig = figure('Name', 'Flow Vectors Visualization', 'NumberTitle', 'off', 'Position', [100, 100, 700, 600]);
    set(fig, 'Color', 'white'); % White background
    set(gca, 'FontName', 'Arial', 'FontSize', 10, 'Box', 'on', 'LineWidth', 1.2);
    
    % Create elegant color scheme
    fullConstraintColor = [0.2, 0.6, 0.8]; % Blue: flow meeting all constraints
    pathConstraintColor = [0.8, 0.4, 0.2]; % Orange: flow meeting only path constraints
    selectedColor = [0.2, 0.8, 0.4];       % Green: selected flows
    feasibleColor = [0.2, 0.7, 0.9];       % Bright blue: flows meeting cost upper limit
    infeasibleColor = [0.9, 0.3, 0.3];     % Red: flows not meeting cost upper limit
    
    % 计算时间和金钱成本，用于绘制成本上限
    allPathTimeCosts = [];
    allPathMoneyCosts = [];
    flowCostsMap = containers.Map('KeyType', 'double', 'ValueType', 'any');
    
    % Define constants for cost calculation
    money = [20, 15, 1, 0, 0, 0, 0, 1];
    freeFlowTime = [18,22.5,12,24,2.4,6,24,12];
    maxCapacity = [3600,3600,1800,1800,1800,1800,1800,1800];
    money = money * relationMatrix';  % 1 x n
    
    % 计算每个流量向量的成本
    if ~isempty(totalValidFlow)
        for i = 1:size(totalValidFlow, 1)
            currentFlow = totalValidFlow(i, :);
            
            % 计算时间成本
            x = currentFlow * relationMatrix;
            RT = calculateRealTime(x, relationMatrix, freeFlowTime, maxCapacity);
            
            % 获取每条路径的时间和金钱成本
            pathTimeCosts = RT(1, :);  % 1 x n
            pathMoneyCosts = money;    % 1 x n
            
            % 移除零成本路径
            validPaths = pathTimeCosts > 0;
            validPathMoneyCosts = pathMoneyCosts(validPaths);
            validPathTimeCosts = pathTimeCosts(validPaths);
            
            % 存储所有路径数据
            allPathTimeCosts = [allPathTimeCosts; validPathTimeCosts'];
            allPathMoneyCosts = [allPathMoneyCosts; validPathMoneyCosts'];
            
            % 组合和排序成本（按金钱成本排序）
            costs = [validPathTimeCosts', validPathMoneyCosts'];
            [~, sortIdx] = sort(validPathMoneyCosts);
            costs = costs(sortIdx, :);
            
            % 存储该流量向量的成本
            flowCostsMap(i) = costs;
        end
    end
    
    % 计算成本上限
    % 按金钱成本分组，为每个金钱成本找到时间的最小值和最大值
    % 找到所有唯一的金钱成本值
    if ~isempty(allPathMoneyCosts)
        uniqueMoneyValues = unique(allPathMoneyCosts);
        leftBoundaryX = [];
        leftBoundaryY = [];
        rightBoundaryX = [];
        rightBoundaryY = [];
        upperLimitX = [];
        upperLimitY = [];
        
        % 为每个金钱成本找到对应的最小和最大时间成本
        for i = 1:length(uniqueMoneyValues)
            currMoney = uniqueMoneyValues(i);
            % 找到具有相同金钱成本的所有点
            sameMoneyIdx = abs(allPathMoneyCosts - currMoney) < 0.001;
            
            if sum(sameMoneyIdx) > 0
                timesForMoney = allPathTimeCosts(sameMoneyIdx);
                minTimeForMoney = min(timesForMoney);
                maxTimeForMoney = max(timesForMoney);
                
                % 计算该金钱成本的时间上限（使用中点）
                midTimeForMoney = (minTimeForMoney + maxTimeForMoney) / 2;
                upperLimitX = [upperLimitX; midTimeForMoney];
                upperLimitY = [upperLimitY; currMoney];
                
                % 添加到边界数组
                leftBoundaryX = [leftBoundaryX; minTimeForMoney];
                leftBoundaryY = [leftBoundaryY; currMoney];
                rightBoundaryX = [rightBoundaryX; maxTimeForMoney];
                rightBoundaryY = [rightBoundaryY; currMoney];
            end
        end
        
        % % 按金钱成本排序边界点
        % [~, sortIdx] = sort(leftBoundaryY);
        % leftBoundaryX = leftBoundaryX(sortIdx);
        
        % [~, sortIdx] = sort(rightBoundaryY);
        % rightBoundaryX = rightBoundaryX(sortIdx);
        
        % 排序上限点
        [upperLimitY, sortIdx] = sort(upperLimitY);
        upperLimitX = upperLimitX(sortIdx);
    end
    
    % Prepare data for visualization
    allFlow = [];
    colors = [];
    markerSizes = [];
    
    % Process flows meeting all constraints
    if ~isempty(totalValidFlow)
        allFlow = [allFlow; totalValidFlow];
        colors = [colors; repmat(fullConstraintColor, size(totalValidFlow, 1), 1)];
        markerSizes = [markerSizes; repmat(10, size(totalValidFlow, 1), 1)];
    end
    
    % Process flows meeting only path constraints
    uniquePathFlows = [];
    if ~isempty(totalPathValidFlow)
        % Exclude flows already in totalValidFlow
        if ~isempty(totalValidFlow)
            uniquePathFlows = setdiff(totalPathValidFlow, totalValidFlow, 'rows');
        else
            uniquePathFlows = totalPathValidFlow;
        end
        
        if ~isempty(uniquePathFlows)
            allFlow = [allFlow; uniquePathFlows];
            % colors = [colors; repmat(pathConstraintColor, size(uniquePathFlows, 1), 1)];
            % markerSizes = [markerSizes; repmat(10, size(uniquePathFlows, 1), 1)];
        end
    end
    
    if isempty(allFlow)
        warning('No flow vectors available for visualization');
        return;
    end
    
    % Verify all flow vectors sum to 10000
    sumFlow = sum(allFlow, 2);
    if any(abs(sumFlow - 10000) > 1e-6)
        warning('Some flow vectors do not sum to 10000');
    end
    
    hold on;

    sf=0.3;
    
    % Project to hyperplane
    projectedData = projectToHyperplane(allFlow);
    
    % Define point sizes based on data count for better visibility
    pointSize = min(50, max(20, 500/sqrt(size(projectedData,1))));
    
    % 检查每个流量向量是否满足成本上限
    feasibleFlowIndices = [];
    infeasibleFlowIndices = [];
    
    if ~isempty(upperLimitX) && ~isempty(totalValidFlow)
        for i = 1:size(totalValidFlow, 1)
            if flowCostsMap.isKey(i)
                costs = flowCostsMap(i);
                isPointFeasible = zeros(size(costs, 1), 1);
                
                for j = 1:size(costs, 1)
                    currPoint = costs(j, :);  % [时间成本, 金钱成本]
                    
                    % 找到最接近的金钱成本点
                    [~, idx] = min(abs(upperLimitY - currPoint(2)));
                    
                    % 检查时间成本是否低于或等于上限
                    isPointFeasible(j) = currPoint(1) <= upperLimitX(idx);
                end
                
                % 如果所有点都可行，则整个流量向量可行
                if all(isPointFeasible)
                    feasibleFlowIndices = [feasibleFlowIndices, i];
                else
                    infeasibleFlowIndices = [infeasibleFlowIndices, i];
                end
            end
        end
    end
    
    % Create empty arrays for legend handles and names
    legendHandles = [];
    legendNames = {};

    % Draw boundaries and regions for different categories
    
    % 1. Project data for different categories
    if ~isempty(totalValidFlow)
        fullIdx = 1:size(totalValidFlow, 1);
        projectedFull = projectedData(fullIdx, :);
        
        % 将可行和不可行的流量分开
        feasibleIdx = ismember(fullIdx, feasibleFlowIndices);
        infeasibleIdx = ismember(fullIdx, infeasibleFlowIndices);
        
        projectedFeasible = projectedFull(feasibleIdx, :);
        projectedInfeasible = projectedFull(infeasibleIdx, :);
    else
        % projectedFull = [];
        projectedFeasible = [];
        projectedInfeasible = [];
    end
    
    if ~isempty(uniquePathFlows)
        pathIdx = size(totalValidFlow,1)+1:size(allFlow,1);
        projectedPath = projectedData(pathIdx, :);
    else
        projectedPath = [];
    end

    % Draw points meeting only path constraints
    if ~isempty(uniquePathFlows) && ~isempty(pathIdx)
        h_path_points = scatter(projectedData(pathIdx, 1), projectedData(pathIdx, 2), ...
            pointSize, pathConstraintColor, 'o', 'filled', ...
            'MarkerFaceAlpha', 0.1, 'MarkerEdgeColor', 'none');
        
        % Add to legend
        legendHandles = [legendHandles, h_path_points];
        legendNames{end+1} = 'Path Constraint Only';
    end
    
    % 绘制不满足成本上限的点
    if ~isempty(projectedInfeasible)
        h_infeasible = scatter(projectedInfeasible(:, 1), projectedInfeasible(:, 2), ...
            pointSize, infeasibleColor, 'o', 'filled', ...
            'MarkerFaceAlpha', 0.7, 'MarkerEdgeColor', 'none');
        
        % Add to legend
        legendHandles = [legendHandles, h_infeasible];
        legendNames{end+1} = 'Not Meeting Cost Upper Limit';
    end
    
    % 绘制满足成本上限的点
    if ~isempty(projectedFeasible)
        h_feasible = scatter(projectedFeasible(:, 1), projectedFeasible(:, 2), ...
            pointSize, feasibleColor, 'o', 'filled', ...
            'MarkerFaceAlpha', 0.7, 'MarkerEdgeColor', 'none');
        
        % Add to legend
        legendHandles = [legendHandles, h_feasible];
        legendNames{end+1} = 'Meeting Cost Upper Limit';
    end

   
    % 2. Draw boundary for "Path Constraint Only" category
    if size(projectedPath, 1) > 3
        try
            % Calculate boundary
            idx_path = boundary(projectedPath(:,1), projectedPath(:,2), sf);
            if ~isempty(idx_path) && length(idx_path) > 3
                % Create boundary polygon
                pathBoundaryX = projectedPath(idx_path,1);
                pathBoundaryY = projectedPath(idx_path,2);
                
                % Fill area with light color
                h_path = patch('XData', pathBoundaryX, 'YData', pathBoundaryY, ...
                    'FaceColor', pathConstraintColor, ... 
                    'FaceAlpha', 0.9, ...
                    'EdgeColor', pathConstraintColor, ...
                    'LineWidth', 1.5);
                
                % Add to legend
                legendHandles = [legendHandles, h_path];
                legendNames{end+1} = 'Path Constraint Region';
            end
        catch e
            fprintf('Boundary calculation error (path): %s\n', e.message);
        end
    end
    
    % 3. 绘制满足成本上限的流量点区域
    if size(projectedFeasible, 1) > 3
        try
            % Calculate boundary
            idx_feasible = boundary(projectedFeasible(:,1), projectedFeasible(:,2), sf);
            if ~isempty(idx_feasible) && length(idx_feasible) > 3
                % Create boundary polygon
                feasibleBoundaryX = projectedFeasible(idx_feasible,1);
                feasibleBoundaryY = projectedFeasible(idx_feasible,2);
                
                % Fill area with light color
                h_feasible_region = patch('XData', feasibleBoundaryX, 'YData', feasibleBoundaryY, ...
                    'FaceColor', feasibleColor, ... 
                    'FaceAlpha', 0.4, ...
                    'EdgeColor', feasibleColor, ...
                    'LineWidth', 1.5);
                
                % Add to legend
                legendHandles = [legendHandles, h_feasible_region];
                legendNames{end+1} = 'Cost Upper Limit Region';
            end
        catch e
            fprintf('Boundary calculation error (feasible): %s\n', e.message);
        end
    end
    
    % 4. 绘制不满足成本上限的流量点区域
    if size(projectedInfeasible, 1) > 3
        try
            % Calculate boundary
            idx_infeasible = boundary(projectedInfeasible(:,1), projectedInfeasible(:,2), sf);
            if ~isempty(idx_infeasible) && length(idx_infeasible) > 3
                % Create boundary polygon
                infeasibleBoundaryX = projectedInfeasible(idx_infeasible,1);
                infeasibleBoundaryY = projectedInfeasible(idx_infeasible,2);
                
                % Fill area with light color
                h_infeasible_region = patch('XData', infeasibleBoundaryX, 'YData', infeasibleBoundaryY, ...
                    'FaceColor', infeasibleColor, ... 
                    'FaceAlpha', 0.4, ...
                    'EdgeColor', infeasibleColor, ...
                    'LineWidth', 1.5);
                
                % Add to legend
                legendHandles = [legendHandles, h_infeasible_region];
                legendNames{end+1} = 'Beyond Cost Upper Limit Region';
            end
        catch e
            fprintf('Boundary calculation error (infeasible): %s\n', e.message);
        end
    end
    
    % Highlight selected points (drawn last to be on top)
    if ~isempty(totalValidFlow) && ~isempty(selectedIndices)
        if max(selectedIndices) <= size(totalValidFlow, 1)
            selectedPoints = totalValidFlow(selectedIndices, :);
            % Project selected points
            projectedSelected = projectToHyperplane(selectedPoints);
            h_selected = scatter(projectedSelected(:, 1), projectedSelected(:, 2), ...
                pointSize*1.5, selectedColor, 'o', 'filled', ...
                'MarkerEdgeColor', 'black', 'LineWidth', 1);
            
            % Add to legend
            legendHandles = [legendHandles, h_selected];
            legendNames{end+1} = 'Selected Flows';
        end
    end
    
    % Add scientific style grid
    grid on;
    grid minor;
    set(gca, 'GridAlpha', 0.1, 'MinorGridAlpha', 0.05, 'Layer', 'top');
    
    % Set title and labels
    title('Flow Vectors Projected on Hyperplane with Cost Upper Limit', 'FontSize', 14, 'FontWeight', 'bold');
    xlabel('Projection Dimension 1', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Projection Dimension 2', 'FontSize', 12, 'FontWeight', 'bold');
    
    % Add legend
    if ~isempty(legendHandles)
        legend(legendHandles, legendNames, ...
            'Location', 'best', ...
            'FontName', 'Arial', 'FontSize', 9, ...
            'EdgeColor', [0.7, 0.7, 0.7], ...
            'Box', 'on');
    end
    
    % Adjust axis limits to provide margin
    axis tight;
    axisLimits = axis;
    axisRange = [axisLimits(2)-axisLimits(1), axisLimits(4)-axisLimits(3)];
    axis([axisLimits(1)-0.03*axisRange(1), axisLimits(2)+0.03*axisRange(1), ...
          axisLimits(3)-0.03*axisRange(2), axisLimits(4)+0.03*axisRange(2)]);
    
    % 输出信息
    if ~isempty(feasibleFlowIndices) || ~isempty(infeasibleFlowIndices)
        fprintf('满足成本上限的流量索引: %s\n', mat2str(feasibleFlowIndices));
        fprintf('超出成本上限的流量索引: %s\n', mat2str(infeasibleFlowIndices));
    end
    
    % Save figure as high-resolution PNG
    figFile = sprintf('results/flow_vectors_with_cost_upper_limit_%s.png', datestr(now, 'yyyymmdd_HHMMSS'));
    print(figFile, '-dpng', '-r300');
end

function projectedData = projectToHyperplane(data)
    % 将数据投影到超平面上
    % 
    % 输入参数:
    %   data - 要投影的数据矩阵，每行一个样本
    %
    % 输出参数:
    %   projectedData - 投影后的二维数据
    
    % 获取数据维度
    [numPoints, n] = size(data);
    
    % 维度检查
    if n < 3
        warning('输入数据维度应大于2以便在超平面上进行二维投影');
        if n == 2
            % 如果是2维，则投影结果是一条线，将第二维设为0
            projectedData = [data(:,1)-data(:,2), zeros(numPoints,1)];
            return;
        elseif n == 1
            % 一维数据，无法投影，直接返回原数据和零列
            projectedData = [data, zeros(numPoints,1)];
            return;
        end
    end
    
    % 将数据转为double类型以保证数值精度
    data = double(data);
    
    % 首先，确保所有点在总和为10000的超平面上
    sumFlow = sum(data, 2);
    if any(abs(sumFlow - 10000) > 1e-6)
        % 法向量：归一化的(1,1,...,1)
        normal = ones(1, n) / sqrt(n);
        
        % 一次性计算所有点到超平面的距离并投影
        distances = (sumFlow - 10000) / sqrt(n);
        data = data - distances * normal;
    end
    
    % 使用质心作为基点（坐标原点）
    basePoint = ones(1,n) * (10000/n);
    
    % 使用null函数构造超平面内的正交基
    % 超平面法向量为(1,1,...,1)，方程为sum(x)=10000
    if n >= 3
        % 对于n>=3，使用null函数获取零空间基
        B = null(ones(1, n));  % 返回 n x (n-1) 矩阵，列是正交基
        basis = B(:, 1:2);     % 取前两列作为基，形成 n x 2 矩阵
    else
        % n=2的情况（前面已经处理了n=1的情况）
        % 手动构造正交基
        basis = [-1, 1; 1, 1] / sqrt(2);  % n x 2 矩阵，正交化的基向量
    end
    
    % 计算相对于基点的向量
    diff = data - basePoint;
    
    % 使用矩阵乘法一次性计算所有投影坐标
    projectedData = diff * basis;  % numPoints x 2 矩阵
end

function M = pairsToMatrix(pairs)
    % 将有序对转换为矩阵
    M = zeros(max(pairs));
    for i = 1:size(pairs,1)
        M(pairs(i,1), pairs(i,2)) = 1;
    end
end

function [reducedCost, reducedFlow] = reduceDataPoints(totalValidCost, totalValidFlow, targetSize)
    % 使用自适应网格采样和保留Pareto前沿减少数据点数量
    % 同时保持数据分布特性和重要点
    
    % 标准化数据到[0,1]区间
    minVals = min(totalValidCost, [], 1);
    maxVals = max(totalValidCost, [], 1);
    normCost = (totalValidCost - minVals) ./ (maxVals - minVals);
    
    % 计算每个维度上需要的网格数量，使总网格数接近目标数量
    gridSize = ceil(sqrt(targetSize));
    
    % 计算每个点所在的网格索引
    gridIndices = floor(normCost * gridSize) + 1;
    
    % 确保网格索引在有效范围内
    gridIndices(gridIndices > gridSize) = gridSize;
    gridIndices(gridIndices < 1) = 1;
    
    % 将二维网格索引转换为一维索引
    cellIndices = sub2ind([gridSize, gridSize], gridIndices(:,1), gridIndices(:,2));
    
    % 对于每个网格，选择一个代表点
    [uniqueCells, ~, ic] = unique(cellIndices);
    
    % 如果网格数量少于目标数量，直接使用所有网格
    if length(uniqueCells) <= targetSize
        fprintf('网格采样后得到 %d 个点\n', length(uniqueCells));
        
        % 为每个网格选择一个代表点
        selectedIndices = zeros(length(uniqueCells), 1);
        for i = 1:length(uniqueCells)
            % 找到属于当前网格的所有点
            gridPoints = find(cellIndices == uniqueCells(i));
            
            % 选择网格中心点作为代表
            if length(gridPoints) > 1
                % 计算网格内所有点的平均位置
                meanPos = mean(normCost(gridPoints,:), 1);
                % 找到最接近平均位置的点
                [~, minIdx] = min(sum((normCost(gridPoints,:) - meanPos).^2, 2));
                selectedIndices(i) = gridPoints(minIdx);
            else
                selectedIndices(i) = gridPoints(1);
            end
        end
    else
        % 如果网格数量超过目标数量，选择更重要的网格
        
        % 计算每个网格的点密度
        cellCounts = histcounts(ic, 1:length(uniqueCells)+1);
        
        % 标记Pareto前沿点所在的网格（保证保留这些网格）
        paretoIndices = findParetoFrontier(totalValidCost);
        paretoCells = unique(cellIndices(paretoIndices));
        
        % 优先选择Pareto前沿点所在的网格
        selectedCells = paretoCells;
        
        % 计算剩余需要的网格数量
        remainingCells = targetSize - length(selectedCells);
        
        if remainingCells > 0
            % 对于非Pareto前沿的网格，按点密度排序
            nonParetoCells = setdiff(uniqueCells, paretoCells);
            nonParetoCellCounts = cellCounts(ismember(uniqueCells, nonParetoCells));
            
            % 按密度排序
            [~, sortIdx] = sort(nonParetoCellCounts, 'descend');
            sortedNonParetoCells = nonParetoCells(sortIdx);
            
            % 选择密度最高的网格
            selectedCells = [selectedCells; sortedNonParetoCells(1:min(remainingCells, length(sortedNonParetoCells)))];
        end
        
        % 为每个选中的网格选择一个代表点
        selectedIndices = zeros(length(selectedCells), 1);
        for i = 1:length(selectedCells)
            % 找到属于当前网格的所有点
            gridPoints = find(cellIndices == selectedCells(i));
            
            % 优先选择Pareto前沿点
            if any(ismember(gridPoints, paretoIndices))
                % 如果网格中有Pareto前沿点，选择第一个
                paretoPoints = intersect(gridPoints, paretoIndices);
                selectedIndices(i) = paretoPoints(1);
            else
                % 否则选择网格中心点作为代表
                if length(gridPoints) > 1
                    % 计算网格内所有点的平均位置
                    meanPos = mean(normCost(gridPoints,:), 1);
                    % 找到最接近平均位置的点
                    [~, minIdx] = min(sum((normCost(gridPoints,:) - meanPos).^2, 2));
                    selectedIndices(i) = gridPoints(minIdx);
                else
                    selectedIndices(i) = gridPoints(1);
                end
            end
        end
    end
    
    % 提取结果
    reducedCost = totalValidCost(selectedIndices, :);
    reducedFlow = totalValidFlow(selectedIndices, :);
    
    fprintf('数据点成功从 %d 减少到 %d\n', size(totalValidCost, 1), length(selectedIndices));
end

function paretoIndices = findParetoFrontier(points)
    % 找到Pareto前沿点的索引
    n = size(points, 1);
    isDominated = false(n, 1);
    
    % 对于每个点，检查它是否被其他点支配
    for i = 1:n
        for j = 1:n
            if i ~= j && all(points(j,:) <= points(i,:)) && any(points(j,:) < points(i,:))
                isDominated(i) = true;
                break;
            end
        end
    end
    
    % 返回非支配点的索引
    paretoIndices = find(~isDominated);
end

function plot3DFlowVectors(totalValidFlow, totalPathValidFlow, relationMatrix, selectedIndices, subset_index)
    % 绘制三维流量向量可视化
    % 
    % 输入参数:
    %   totalValidFlow     - 满足所有约束的流量向量
    %   totalPathValidFlow - 仅满足路径约束的流量向量
    %   relationMatrix     - 关系矩阵
    %   selectedIndices    - 选择的流量向量索引
    %   subset_index       - 子集索引，用于确定使用哪些路径
    
    % 检查是否有有效流量数据
    if isempty(totalValidFlow)
        warning('没有满足所有约束的流量数据可供可视化');
        return;
    end
    
    % 创建科学风格的图形
    fig = figure('Name', '3D Flow Vectors Visualization', 'NumberTitle', 'off', 'Position', [150, 150, 800, 700]);
    set(fig, 'Color', 'white'); % 白色背景
    set(gca, 'FontName', 'Arial', 'FontSize', 10, 'Box', 'on', 'LineWidth', 1.2);
    
    % 创建优雅的颜色方案
    fullConstraintColor = [0.2, 0.6, 0.8]; % 蓝色: 满足所有约束的流量
    pathConstraintColor = [0.8, 0.4, 0.2]; % 橙色: 仅满足路径约束的流量
    selectedColor = [0.2, 0.8, 0.4];       % 绿色: 选中的流量
    feasibleColor = [0.2, 0.7, 0.9];       % 亮蓝色: 满足成本上限的流量
    infeasibleColor = [0.9, 0.3, 0.3];     % 红色: 不满足成本上限的流量
    
    % 计算时间和金钱成本，用于绘制成本上限
    allPathTimeCosts = [];
    allPathMoneyCosts = [];
    flowCostsMap = containers.Map('KeyType', 'double', 'ValueType', 'any');
    
    % Define constants for cost calculation
    money = [20, 15, 1, 0, 0, 0, 0, 1];
    freeFlowTime = [18,22.5,12,24,2.4,6,24,12];
    maxCapacity = [3600,3600,1800,1800,1800,1800,1800,1800];
    money = money * relationMatrix';  % 1 x n
    
    % 计算每个流量向量的成本
    if ~isempty(totalValidFlow)
        for i = 1:size(totalValidFlow, 1)
            currentFlow = totalValidFlow(i, :);
            
            % 计算时间成本
            x = currentFlow * relationMatrix;
            RT = calculateRealTime(x, relationMatrix, freeFlowTime, maxCapacity);
            
            % 获取每条路径的时间和金钱成本
            pathTimeCosts = RT(1, :);  % 1 x n
            pathMoneyCosts = money;    % 1 x n
            
            % 移除零成本路径
            validPaths = pathTimeCosts > 0;
            validPathMoneyCosts = pathMoneyCosts(validPaths);
            validPathTimeCosts = pathTimeCosts(validPaths);
            
            % 存储所有路径数据
            allPathTimeCosts = [allPathTimeCosts; validPathTimeCosts'];
            allPathMoneyCosts = [allPathMoneyCosts; validPathMoneyCosts'];
            
            % 组合和排序成本（按金钱成本排序）
            costs = [validPathTimeCosts', validPathMoneyCosts'];
            [~, sortIdx] = sort(validPathMoneyCosts);
            costs = costs(sortIdx, :);
            
            % 存储该流量向量的成本
            flowCostsMap(i) = costs;
        end
    end
    
    % 计算成本上限
    % 按金钱成本分组，为每个金钱成本找到时间的最小值和最大值
    % 找到所有唯一的金钱成本值
    upperLimitX = [];
    upperLimitY = [];
    
    if ~isempty(allPathMoneyCosts)
        uniqueMoneyValues = unique(allPathMoneyCosts);
        leftBoundaryX = [];
        leftBoundaryY = [];
        rightBoundaryX = [];
        rightBoundaryY = [];
        
        % 为每个金钱成本找到对应的最小和最大时间成本
        for i = 1:length(uniqueMoneyValues)
            currMoney = uniqueMoneyValues(i);
            % 找到具有相同金钱成本的所有点
            sameMoneyIdx = abs(allPathMoneyCosts - currMoney) < 0.001;
            
            if sum(sameMoneyIdx) > 0
                timesForMoney = allPathTimeCosts(sameMoneyIdx);
                minTimeForMoney = min(timesForMoney);
                maxTimeForMoney = max(timesForMoney);
                
                % 计算该金钱成本的时间上限（使用中点）
                midTimeForMoney = (minTimeForMoney + maxTimeForMoney) / 2;
                upperLimitX = [upperLimitX; midTimeForMoney];
                upperLimitY = [upperLimitY; currMoney];
                
                % 添加到边界数组
                leftBoundaryX = [leftBoundaryX; minTimeForMoney];
                leftBoundaryY = [leftBoundaryY; currMoney];
                rightBoundaryX = [rightBoundaryX; maxTimeForMoney];
                rightBoundaryY = [rightBoundaryY; currMoney];
            end
        end
        
        % 按金钱成本排序边界点
        [upperLimitY, sortIdx] = sort(upperLimitY);
        upperLimitX = upperLimitX(sortIdx);
    end
    
    % 检查每个流量向量是否满足成本上限
    feasibleFlowIndices = [];
    infeasibleFlowIndices = [];
    
    if ~isempty(upperLimitX) && ~isempty(totalValidFlow)
        for i = 1:size(totalValidFlow, 1)
            if flowCostsMap.isKey(i)
                costs = flowCostsMap(i);
                isPointFeasible = zeros(size(costs, 1), 1);
                
                for j = 1:size(costs, 1)
                    currPoint = costs(j, :);  % [时间成本, 金钱成本]
                    
                    % 找到最接近的金钱成本点
                    [~, idx] = min(abs(upperLimitY - currPoint(2)));
                    
                    % 检查时间成本是否低于或等于上限
                    isPointFeasible(j) = currPoint(1) <= upperLimitX(idx);
                end
                
                % 如果所有点都可行，则整个流量向量可行
                if all(isPointFeasible)
                    feasibleFlowIndices = [feasibleFlowIndices, i];
                else
                    infeasibleFlowIndices = [infeasibleFlowIndices, i];
                end
            end
        end
    end
    
    % 定义所有可能的路径组合
    pathCombinations = {};
    pathCombinationNames = {};
    
    if subset_index == 0
        % 对于subset_index=0，使用path1, path2, path5
        pathCombinations{1} = [1, 2, 3];
        pathCombinationNames{1} = 'Path 1,2,5';
    elseif subset_index == 1
        % 对于subset_index=1，使用不同的路径组合
        pathCombinations{1} = [1, 2, 3, 5];
        pathCombinations{2} = [1, 2, 4, 5];
        pathCombinations{3} = [1, 3, 5];
        pathCombinations{4} = [2, 3, 5];
        
        pathCombinationNames{1} = 'Path 1,2,3,5';
        pathCombinationNames{2} = 'Path 1,2,4,5';
        pathCombinationNames{3} = 'Path 1,3,5';
        pathCombinationNames{4} = 'Path 2,3,5';
    end
    
    % 如果没有组合，使用默认值
    if isempty(pathCombinations)
        pathCombinations{1} = [1, 2, 3];
        pathCombinationNames{1} = 'Path 1,2,3';
    end
    
    % 创建图形控制面板
    panel = uipanel('Title', '可视化控制', 'Position', [0.01, 0.01, 0.35, 0.3], ...
        'FontSize', 10, 'FontWeight', 'bold');
    
    % 面板内部边距
    padding = 10;
    controlHeight = 22;
    labelWidth = 80;
    controlWidth = 120;
    
    % 计算垂直位置（从上往下排列）
    panelHeight = panel.Position(4) * fig.Position(4);
    row1 = panelHeight - padding - controlHeight;
    row2 = row1 - controlHeight - padding/2;
    row3 = row2 - controlHeight - padding/2;
    row4 = row3 - controlHeight - padding/2;
    row5 = row4 - controlHeight - padding/2;
    
    % ===== 第一行：路径组合选择 =====
    uicontrol('Parent', panel, 'Style', 'text', 'String', '路径组合:', ...
        'Position', [padding, row1, labelWidth, controlHeight], ...
        'HorizontalAlignment', 'left');
    pathComboDropdown = uicontrol('Parent', panel, 'Style', 'popupmenu', ...
        'String', pathCombinationNames, 'Position', [padding+labelWidth, row1, controlWidth, controlHeight], ...
        'Callback', @updatePathCombination);
    
    % ===== 第二行：显示选项 =====
    % 成本上限显示控制
    uicontrol('Parent', panel, 'Style', 'checkbox', 'String', '显示成本上限区分', ...
        'Position', [padding, row2, controlWidth, controlHeight], ...
        'Value', 1, 'Callback', @toggleCostLimit);
    
    % 显示内部线条控制
    uicontrol('Parent', panel, 'Style', 'checkbox', 'String', '显示内部线条', ...
        'Position', [padding+labelWidth+50, row2, controlWidth, controlHeight], ...
        'Value', 0, 'Callback', @toggleInnerLines);
    
    % ===== 第三行：边界类型控制 =====
    % 边界类型标签
    uicontrol('Parent', panel, 'Style', 'text', 'String', '边界类型:', ...
        'Position', [padding, row3, labelWidth, controlHeight], ...
        'HorizontalAlignment', 'left');
    
    % 添加凸包选项
    uicontrol('Parent', panel, 'Style', 'radiobutton', 'String', '凸包', ...
        'Position', [padding+labelWidth, row3, 70, controlHeight], 'Tag', 'convex', ...
        'Value', 1, 'Callback', @(src,~)selectBoundaryType(src,'convex'));
    
    % 添加Alpha形状选项
    uicontrol('Parent', panel, 'Style', 'radiobutton', 'String', 'Alpha形状', ...
        'Position', [padding+labelWidth+80, row3, 100, controlHeight], 'Tag', 'alpha', ...
        'Value', 0, 'Callback', @(src,~)selectBoundaryType(src,'alpha'));
    
    % ===== 第四行：Alpha值滑动条 =====
    % Alpha值滑动条
    uicontrol('Parent', panel, 'Style', 'text', 'String', 'Alpha值:', ...
        'Position', [padding, row4, labelWidth, controlHeight], ...
        'HorizontalAlignment', 'left');
    alphaSlider = uicontrol('Parent', panel, 'Style', 'slider', ...
        'Min', 0.1, 'Max', 3, 'Value', 1.0, ...
        'Position', [padding+labelWidth, row4, controlWidth, controlHeight], ...
        'Callback', @updateAlphaValue);
    
    % ===== 第五行：操作按钮 =====
    % 添加旋转控制按钮
    uicontrol('Parent', panel, 'Style', 'pushbutton', 'String', '旋转视图', ...
        'Position', [padding, row5, labelWidth, controlHeight], 'Callback', @toggleRotation);
    
    % 添加保存图像按钮
    uicontrol('Parent', panel, 'Style', 'pushbutton', 'String', '保存图像', ...
        'Position', [padding+labelWidth+10, row5, labelWidth, controlHeight], 'Callback', @saveImage);
    
    % 初始化旋转标志
    isRotating = false;
    rotationTimer = [];
    
    % 初始化边界类型和Alpha值
    boundaryType = 'convex';
    alphaValue = 1.0;
    showCostLimit = true;
    showInnerLines = false;
    
    % 边界类型选择回调函数
    function selectBoundaryType(src, type)
        % 如果当前按钮被选中，取消选择其他按钮
        if src.Value == 1
            boundaryType = type;
            % 查找并取消选择其他单选按钮
            btns = findobj(panel, 'Style', 'radiobutton');
            for i = 1:length(btns)
                if ~strcmp(btns(i).Tag, type)
                    btns(i).Value = 0;
                end
            end
            updatePathCombination();
        else
            % 防止取消选择当前按钮（必须始终有一个选中）
            src.Value = 1;
        end
    end
    
    % 绘制初始组合的可视化
    updatePathCombination();
    
    % 回调函数：切换成本上限显示
    function toggleCostLimit(src, ~)
        showCostLimit = get(src, 'Value') == 1;
        updatePathCombination();
    end
    
    % 回调函数：切换内部线条显示
    function toggleInnerLines(src, ~)
        showInnerLines = get(src, 'Value') == 1;
        updatePathCombination();
    end
    
    % 回调函数：更新Alpha值
    function updateAlphaValue(~, ~)
        alphaValue = get(alphaSlider, 'Value');
        updatePathCombination();
    end
    
    % 回调函数：更新路径组合
    function updatePathCombination(~, ~)
        % 清除当前图形
        cla;
        
        % 获取选择的路径组合
        comboIdx = get(pathComboDropdown, 'Value');
        pathIndices = pathCombinations{comboIdx};
        
        % 获取路径名称
        pathNames = cell(1, length(pathIndices));
        for i = 1:length(pathIndices)
            pathNames{i} = sprintf('Path %d', pathIndices(i));
        end
        
        % 处理满足所有约束的流量
        if ~isempty(totalValidFlow)
            % 提取选择的路径数据
            validFlowSelected = totalValidFlow(:, pathIndices);
            
            % 处理仅满足路径约束的流量
            uniquePathFlows = [];
            pathFlowSelected = [];
            if ~isempty(totalPathValidFlow)
                % 排除已经在totalValidFlow中的流量
                uniquePathFlows = setdiff(totalPathValidFlow, totalValidFlow, 'rows');
                
                if ~isempty(uniquePathFlows)
                    % 提取选择的路径数据
                    pathFlowSelected = uniquePathFlows(:, pathIndices);
                end
            end
            
            % 准备3D可视化
            hold on;
            
            % 创建图例句柄和名称数组
            legendHandles = [];
            legendNames = {};
            
            % 根据成本上限区分流量
            feasibleFlow = [];
            infeasibleFlow = [];
            
            if showCostLimit && ~isempty(feasibleFlowIndices)
                % 分离满足和不满足成本上限的流量
                feasibleMask = ismember(1:size(totalValidFlow, 1), feasibleFlowIndices);
                infeasibleMask = ismember(1:size(totalValidFlow, 1), infeasibleFlowIndices);
                
                feasibleFlow = totalValidFlow(feasibleMask, pathIndices);
                infeasibleFlow = totalValidFlow(infeasibleMask, pathIndices);
            else
                feasibleFlow = validFlowSelected; % 如果不显示成本上限，则所有点都视为可行点
            end
            
            % 决定是否显示内部线条
            if showInnerLines
                edgeColor = 'k';  % 黑色边缘
                edgeAlpha = 0.3;  % 半透明
            else
                edgeColor = 'none';  % 无边缘
                edgeAlpha = 0;
            end
            
            % 绘制满足路径约束的边界
            if ~isempty(uniquePathFlows) && ~isempty(pathFlowSelected) && size(pathFlowSelected, 1) >= 4
                % 尝试绘制路径约束点的边界
                try
                    if strcmp(boundaryType, 'convex')
                        % 使用凸包绘制边界
                        K_path = convhull(pathFlowSelected(:, 1), pathFlowSelected(:, 2), pathFlowSelected(:, 3));
                        h_path_hull = trisurf(K_path, pathFlowSelected(:, 1), pathFlowSelected(:, 2), pathFlowSelected(:, 3), ...
                            'FaceColor', pathConstraintColor, ...
                            'FaceAlpha', 0.5, ...
                            'EdgeColor', edgeColor, ...
                            'EdgeAlpha', edgeAlpha);
                        legendHandles = [legendHandles, h_path_hull];
                        legendNames{end+1} = '路径约束边界(凸包)';
                    else
                        % 使用Alpha形状绘制边界（如果可用）
                        if exist('alphaShape', 'file')
                            % 使用Alpha形状
                            shp = alphaShape(pathFlowSelected(:, 1), pathFlowSelected(:, 2), pathFlowSelected(:, 3), ...
                                alphaValue * criticalAlpha(pathFlowSelected));
                            h_path_alpha = plot(shp, 'FaceColor', pathConstraintColor, 'FaceAlpha', 0.5, ...
                                'EdgeColor', edgeColor, 'EdgeAlpha', edgeAlpha);
                            legendHandles = [legendHandles, h_path_alpha];
                            legendNames{end+1} = '路径约束边界(Alpha形状)';
                        else
                            % 如果Alpha形状不可用，退回到凸包
                            K_path = convhull(pathFlowSelected(:, 1), pathFlowSelected(:, 2), pathFlowSelected(:, 3));
                            h_path_hull = trisurf(K_path, pathFlowSelected(:, 1), pathFlowSelected(:, 2), pathFlowSelected(:, 3), ...
                                'FaceColor', pathConstraintColor, ...
                                'FaceAlpha', 0.5, ...
                                'EdgeColor', edgeColor, ...
                                'EdgeAlpha', edgeAlpha);
                            legendHandles = [legendHandles, h_path_hull];
                            legendNames{end+1} = '路径约束边界(凸包)';
                            warning('Alpha形状函数不可用，已退回到凸包');
                        end
                    end
                catch e
                    fprintf('路径约束边界计算错误: %s\n', e.message);
                end
            end
            
            % 绘制不满足成本上限的边界
            if showCostLimit && ~isempty(infeasibleFlow) && size(infeasibleFlow, 1) >= 4
                % 尝试绘制不满足成本上限点的边界
                try
                    if strcmp(boundaryType, 'convex') && size(infeasibleFlow, 1) > 3
                        % 使用凸包绘制边界
                        K_infeasible = convhull(infeasibleFlow(:, 1), infeasibleFlow(:, 2), infeasibleFlow(:, 3));
                        h_infeasible_hull = trisurf(K_infeasible, infeasibleFlow(:, 1), infeasibleFlow(:, 2), infeasibleFlow(:, 3), ...
                            'FaceColor', infeasibleColor, ...
                            'FaceAlpha', 0.5, ...
                            'EdgeColor', edgeColor, ...
                            'EdgeAlpha', edgeAlpha);
                        legendHandles = [legendHandles, h_infeasible_hull];
                        legendNames{end+1} = '不满足成本上限区域';
                    elseif size(infeasibleFlow, 1) > 3
                        % 使用Alpha形状绘制边界（如果可用）
                        if exist('alphaShape', 'file')
                            % 计算推荐的Alpha值并应用用户调整因子
                            alpha = alphaValue * criticalAlpha(infeasibleFlow);
                            % 创建Alpha形状
                            shp = alphaShape(infeasibleFlow(:, 1), infeasibleFlow(:, 2), infeasibleFlow(:, 3), alpha);
                            h_infeasible_hull = plot(shp, 'FaceColor', infeasibleColor, 'FaceAlpha', 0.5, ...
                                'EdgeColor', edgeColor, 'EdgeAlpha', edgeAlpha);
                            legendHandles = [legendHandles, h_infeasible_hull];
                            legendNames{end+1} = '不满足成本上限区域(Alpha形状)';
                        end
                    end
                catch e
                    fprintf('不满足成本上限边界计算错误: %s\n', e.message);
                end
            end
            
            % 绘制满足成本上限的边界
            if ~isempty(feasibleFlow) && size(feasibleFlow, 1) >= 4
                % 绘制满足成本上限点的边界
                try
                    if strcmp(boundaryType, 'convex') && size(feasibleFlow, 1) > 3
                        % 使用凸包绘制边界
                        K = convhull(feasibleFlow(:, 1), feasibleFlow(:, 2), feasibleFlow(:, 3));
                        h_hull = trisurf(K, feasibleFlow(:, 1), feasibleFlow(:, 2), feasibleFlow(:, 3), ...
                            'FaceColor', feasibleColor, ...
                            'FaceAlpha', 0.5, ...
                            'EdgeColor', edgeColor, ...
                            'EdgeAlpha', edgeAlpha);
                        legendHandles = [legendHandles, h_hull];
                        if showCostLimit
                            legendNames{end+1} = '满足成本上限区域';
                        else
                            legendNames{end+1} = '可行区域边界(凸包)';
                        end
                    elseif size(feasibleFlow, 1) > 3
                        % 使用Alpha形状绘制边界（如果可用）
                        if exist('alphaShape', 'file')
                            % 计算推荐的Alpha值并应用用户调整因子
                            alpha = alphaValue * criticalAlpha(feasibleFlow);
                            % 创建Alpha形状
                            shp = alphaShape(feasibleFlow(:, 1), feasibleFlow(:, 2), feasibleFlow(:, 3), alpha);
                            h_hull = plot(shp, 'FaceColor', feasibleColor, 'FaceAlpha', 0.5, ...
                                'EdgeColor', edgeColor, 'EdgeAlpha', edgeAlpha);
                            legendHandles = [legendHandles, h_hull];
                            if showCostLimit
                                legendNames{end+1} = '满足成本上限区域(Alpha形状)';
                            else
                                legendNames{end+1} = '可行区域边界(Alpha形状)';
                            end
                            
                            % 显示当前Alpha值
                            titleStr = sprintf('三维流量向量边界可视化 (Alpha=%.2f)', alpha);
                            title(titleStr, 'FontSize', 14, 'FontWeight', 'bold');
                        end
                    end
                catch e
                    fprintf('边界计算错误: %s\n', e.message);
                end
            end
            
            % 确保添加坐标轴标签，并设置字体加粗和增大字号
            if length(pathIndices) >= 3
                xlabel(pathNames{1}, 'FontSize', 14, 'FontWeight', 'bold');
                ylabel(pathNames{2}, 'FontSize', 14, 'FontWeight', 'bold');
                zlabel(pathNames{3}, 'FontSize', 14, 'FontWeight', 'bold');
            end
            
            % 设置标题（如果没有在前面设置）
            if ~strcmp(boundaryType, 'alpha') || ~exist('alphaShape', 'file')
                if showCostLimit
                    titlePrefix = '三维流量向量边界可视化(含成本上限)';
                else
                    titlePrefix = '三维流量向量边界可视化';
                end
                title(titlePrefix, 'FontSize', 14, 'FontWeight', 'bold');
            end
            
            % 如果有4个维度，则添加说明
            if length(pathIndices) > 3
                % 计算第4维度的流量值（总和为10000）
                fourthDimValues = 10000 - sum(validFlowSelected(:,1:3), 2);
                
                % 检查所有值是否相同
                if all(abs(fourthDimValues - fourthDimValues(1)) < 1e-6)
                    if showCostLimit
                        titlePrefix = '三维流量向量边界可视化(含成本上限)';
                    else
                        titlePrefix = '三维流量向量边界可视化';
                    end
                    title([sprintf('%s (%s=%.0f)', titlePrefix, pathNames{4}, fourthDimValues(1))], ...
                        'FontSize', 14, 'FontWeight', 'bold');
                else
                    if showCostLimit
                        titlePrefix = '三维流量向量边界可视化(含成本上限)';
                    else
                        titlePrefix = '三维流量向量边界可视化';
                    end
                    title([sprintf('%s (%s值为变量)', titlePrefix, pathNames{4})], ...
                        'FontSize', 14, 'FontWeight', 'bold');
                end
            end
            
            % 添加图例，确保所有图形句柄都有效
            validHandles = ishandle(legendHandles);
            if any(validHandles) && ~isempty(legendNames)
                try
                    legend(legendHandles(validHandles), legendNames(validHandles), ...
                        'Location', 'best', ...
                        'FontName', 'Arial', 'FontSize', 11, ...
                        'EdgeColor', [0.7, 0.7, 0.7], ...
                        'Box', 'on');
                catch e
                    warning(e.identifier, '图例创建错误: %s', e.message);
                end
            end
            
            % 添加网格
            grid on;
            
            % 设置视图角度
            view(45, 30);
            
            % 调整坐标轴限制以提供边距
            axis tight;
            axisLimits = axis;
            axisRange = [axisLimits(2)-axisLimits(1), axisLimits(4)-axisLimits(3), axisLimits(6)-axisLimits(5)];
            axis([axisLimits(1)-0.03*axisRange(1), axisLimits(2)+0.03*axisRange(1), ...
                  axisLimits(3)-0.03*axisRange(2), axisLimits(4)+0.03*axisRange(2), ...
                  axisLimits(5)-0.03*axisRange(3), axisLimits(6)+0.03*axisRange(3)]);
            
            % 添加一个文本注释，说明约束条件
            if showCostLimit && (~isempty(feasibleFlowIndices) || ~isempty(infeasibleFlowIndices))
                sumText = sprintf('总流量约束: 所有路径流量之和 = 10000 | 满足成本上限的流量: %d个 | 不满足成本上限的流量: %d个', ...
                    length(feasibleFlowIndices), length(infeasibleFlowIndices));
            else
                sumText = sprintf('总流量约束: 所有路径流量之和 = 10000');
            end
            
            hold off;
        else
            warning('没有满足所有约束的流量数据可供可视化');
        end
    end
    
    % 辅助函数：计算Alpha形状的临界Alpha值
    function critAlpha = criticalAlpha(points)
        % 计算适合给定点云的推荐Alpha值
        if size(points, 1) < 4
            critAlpha = 1;  % 默认值
            return;
        end
        
        % 计算点云的典型尺度
        xRange = max(points(:,1)) - min(points(:,1));
        yRange = max(points(:,2)) - min(points(:,2));
        zRange = max(points(:,3)) - min(points(:,3));
        
        % 使用点云尺度的10%作为推荐Alpha值
        critAlpha = 0.1 * mean([xRange, yRange, zRange]);
    end
    
    % 回调函数：切换旋转
    function toggleRotation(~, ~)
        if isRotating
            % 停止旋转
            stop(rotationTimer);
            delete(rotationTimer);
            isRotating = false;
        else
            % 开始旋转
            isRotating = true;
            rotationTimer = timer('ExecutionMode', 'fixedRate', 'Period', 0.1, ...
                'TimerFcn', @rotateView);
            start(rotationTimer);
        end
    end
    
    % 旋转函数
    function rotateView(~, ~)
        % 获取当前视角
        [az, el] = view;
        % 更新方位角
        az = az + 2;
        if az > 360
            az = az - 360;
        end
        % 应用新视角
        view(az, el);
        drawnow;
    end
    
    % 保存图像函数
    function saveImage(~, ~)
        % 生成文件名
        comboIdx = get(pathComboDropdown, 'Value');
        pathStr = strrep(pathCombinationNames{comboIdx}, ',', '_');
        pathStr = strrep(pathStr, ' ', '');
        if showCostLimit
            costLimitStr = 'with_cost_limit';
        else
            costLimitStr = '';
        end
        figFile = sprintf('results/3d_flow_boundaries_%s_%s_%s.png', pathStr, costLimitStr, datestr(now, 'yyyymmdd_HHMMSS'));
        
        % 保存图像
        print(figFile, '-dpng', '-r300');
        msgbox(['图像已保存为: ' figFile], '保存成功');
    end
    
    % 清理函数
    set(fig, 'DeleteFcn', @cleanupFunc);
    function cleanupFunc(~, ~)
        % 如果旋转定时器存在，停止并删除它
        if ~isempty(rotationTimer) && isvalid(rotationTimer)
            stop(rotationTimer);
            delete(rotationTimer);
        end
    end
end