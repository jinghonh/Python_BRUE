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
    cacheFileName = sprintf('cache_zeta%d_subset%d.mat', zeta, subset_index);
    pathConstraintCacheFileName = sprintf('cache_path_only_zeta%d_subset%d.mat', zeta, subset_index);
    
    if exist(cacheFileName, 'file') && exist(pathConstraintCacheFileName, 'file')
        fprintf('加载缓存数据...\n');
        load(cacheFileName, 'totalValidCost', 'totalValidFlow', 'relationMatrix');
        load(pathConstraintCacheFileName, 'totalPathValidCost', 'totalPathValidFlow');
        fprintf('加载缓存数据完成...\n');
        toc;
        
        % 直接绘图
        q = 200;
        if ~isempty(totalValidFlow)
            selectedIndices = randperm(size(totalValidFlow, 1), min(q, size(totalValidFlow, 1)));
            % plotResults(totalValidCost, selectedIndices);
            % plotPathCosts(totalValidFlow, relationMatrix, selectedIndices);
            % 绘制流量向量可视化
            plotFlowVectors(totalValidFlow, totalPathValidFlow, relationMatrix, selectedIndices);
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
    minIt = 10;
    maxIt = 40;
    
    % 设置数据压缩参数
    maxPointsPerIteration = 5e4; % 每次迭代最大保留点数
    
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
    [totalValidCost,Ia,Ic]=unique(totalValidCost,"rows");
    totalValidFlow = totalValidFlow(Ia,:);
    
    % 使用网格采样减少数据点数量
    targetSize = 5e5; % 目标数据点数量：十万级别
    if size(totalValidCost, 1) > targetSize
        fprintf('正在进行全约束数据压缩，从 %d 个点压缩到目标 %d 个点左右...\n', size(totalValidCost, 1), targetSize);
        [totalValidCost, totalValidFlow] = reduceDataPoints(totalValidCost, totalValidFlow, targetSize);
    end
    
    % 保存满足路径约束的结果
    [totalPathValidCost,Ia,Ic]=unique(totalPathValidCost,"rows");
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
        plotResults(totalValidCost, selectedIndices);
        
        % 调用路径成本绘图函数
        plotPathCosts(totalValidFlow, relationMatrix, selectedIndices);
        
        % 绘制流量向量可视化
        plotFlowVectors(totalValidFlow, totalPathValidFlow, relationMatrix, selectedIndices);
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
    
    % 估计结果大小并预分配内存
    estSize = 1;
    for i = 1:min(n-1, numel(samples))
        estSize = estSize * length(samples{i});
    end
    
    % 限制大小以避免内存溢出
    if estSize > 1e7
        estSize = 1e7;
    end
    
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
    % 检查路径约束
    combs = nchoosek(1:n, 2);
    RTDiffs = abs(RT(:, combs(:,1)) - RT(:, combs(:,2)));
    
    for i = 1:size(RTDiffs,2)
        err = err + max(0, RTDiffs(:,i) - zeta);
    end
end

function err = checkMoneyConstraints(RT, money, n, err)
    % 检查货币约束
    for i = 1:n
        for j = i+1:n
            err = err + max(0, (RT(:,i)- RT(:,j)) .* (money(i) - money(j)));
        end
    end
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
            colors = [colors; repmat(pathConstraintColor, size(uniquePathFlows, 1), 1)];
            markerSizes = [markerSizes; repmat(10, size(uniquePathFlows, 1), 1)];
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
    
    % Create empty arrays for legend handles and names
    legendHandles = [];
    legendNames = {};

    % Draw boundaries and regions for different categories
    
    % 1. Project data for different categories
    if ~isempty(totalValidFlow)
        fullIdx = 1:size(totalValidFlow, 1);
        projectedFull = projectedData(fullIdx, :);
    else
        projectedFull = [];
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
            'MarkerFaceAlpha', 0.7, 'MarkerEdgeColor', 'none');
        
        % Add to legend
        legendHandles = [legendHandles, h_path_points];
        legendNames{end+1} = 'Path Constraint Only';
    end
    
    % Draw points meeting all constraints
    if ~isempty(totalValidFlow)
        h_full_points = scatter(projectedData(fullIdx, 1), projectedData(fullIdx, 2), ...
            pointSize, fullConstraintColor, 'o', 'filled', ...
            'MarkerFaceAlpha', 0.7, 'MarkerEdgeColor', 'none');
        
        % Add to legend
        legendHandles = [legendHandles, h_full_points];
        legendNames{end+1} = 'All Constraints';
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
    
    % 3. Draw boundary for "All Constraints" category
    if size(projectedFull, 1) > 3
        try
            % Calculate boundary
            idx_full = boundary(projectedFull(:,1), projectedFull(:,2), sf);
            if ~isempty(idx_full) && length(idx_full) > 3
                % Create boundary polygon
                fullBoundaryX = projectedFull(idx_full,1);
                fullBoundaryY = projectedFull(idx_full,2);
                
                % Fill area with light color
                h_full = patch('XData', fullBoundaryX, 'YData', fullBoundaryY, ...
                    'FaceColor', fullConstraintColor, ... 
                    'FaceAlpha', 0.8, ...
                    'EdgeColor', fullConstraintColor, ...
                    'LineWidth', 1.5);
                
                % Add to legend
                legendHandles = [legendHandles, h_full];
                legendNames{end+1} = 'All Constraints Region';
            end
        catch e
            fprintf('Boundary calculation error (full): %s\n', e.message);
        end
    end
    

    
    % % 4. Draw overall boundary if needed
    % if size(projectedData, 1) > 3 && ~isempty(legendHandles)
    %     try
    %         % Calculate overall boundary
    %         idx_all = boundary(projectedData(:,1), projectedData(:,2), sf);
    %         if ~isempty(idx_all) && length(idx_all) > 3
    %             % Create boundary polygon
    %             boundaryX = projectedData(idx_all,1);
    %             boundaryY = projectedData(idx_all,2);
    % 
    %             % Draw boundary line only (no fill)
    %             h_all = plot(boundaryX, boundaryY, 'Color', [0.4, 0.5, 0.8], 'LineWidth', 1.0, 'LineStyle', '-');
    % 
    %             % Add to legend
    %             legendHandles = [legendHandles, h_all];
    %             legendNames{end+1} = 'Overall Boundary';
    %         end
    %     catch e
    %         fprintf('Boundary calculation error (overall): %s\n', e.message);
    %     end
    % end
    
    % 5. Draw points on top of regions
    
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
    title('Flow Vectors Projected on Hyperplane', 'FontSize', 14, 'FontWeight', 'bold');
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
    
    % Save figure as high-resolution PNG
    figFile = sprintf('results/flow_vectors_%s.png', datestr(now, 'yyyymmdd_HHMMSS'));
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