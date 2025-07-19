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
        load(cacheFileName, 'totalValidFlow', 'relationMatrix');
        load(pathConstraintCacheFileName, 'totalPathValidFlow');
        fprintf('加载缓存数据完成...\n');
        toc;
        
        % 直接绘图
        q = 10000;
        if ~isempty(totalValidFlow)
            selectedIndices = randperm(size(totalValidFlow, 1), min(q, size(totalValidFlow, 1)));
            plotPathCosts(totalValidFlow, relationMatrix, selectedIndices);
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
    totalValidFlow = [];
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
    maxPointsPerIteration = 1e6; % 每次迭代最大保留点数

    % 传统网格采样策略
    for ii = minIt:maxIt
        waitbar((ii-minIt)/(maxIt-minIt), h, sprintf('正在计算第 %d/%d 次迭代...', ii-minIt+1, maxIt-minIt+1));
        [samplesMat, totalValidFlow, totalPathValidFlow] = processIteration(ii, n, rangeMin, rangeMax, bound, relationMatrix, totalValidFlow, totalPathValidFlow, zeta);
        
        % 增量数据压缩：当有效点数量超过阈值时进行压缩
        if size(totalValidFlow, 1) > maxPointsPerIteration
            fprintf('迭代 %d: 压缩全部约束数据点，从 %d 个点压缩...\n', ii, size(totalValidFlow, 1));
            totalValidFlow = reduceDataPoints(totalValidFlow, maxPointsPerIteration);
        end
        
        % 增量数据压缩：当路径约束点数量超过阈值时进行压缩
        if size(totalPathValidFlow, 1) > maxPointsPerIteration
            fprintf('迭代 %d: 压缩路径约束数据点，从 %d 个点压缩...\n', ii, size(totalPathValidFlow, 1));
            totalPathValidFlow = reduceDataPoints(totalPathValidFlow, maxPointsPerIteration);
        end
        
        % 更新搜索范围
        [path_err, ~] = evaluateObjective(samplesMat, relationMatrix, zeta);
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
    [totalValidFlow, Ia, ~] = unique(totalValidFlow, "rows");
    
    % 使用网格采样减少数据点数量
    targetSize = 5e5; % 目标数据点数量：十万级别
    if size(totalValidFlow, 1) > targetSize
        fprintf('正在进行全约束数据压缩，从 %d 个点压缩到目标 %d 个点左右...\n', size(totalValidFlow, 1), targetSize);
        totalValidFlow = reduceDataPoints(totalValidFlow, targetSize);
    end
    
    % 保存满足路径约束的结果
    [totalPathValidFlow, Ia, ~] = unique(totalPathValidFlow, "rows");
    
    % 使用网格采样减少数据点数量
    if size(totalPathValidFlow, 1) > targetSize
        fprintf('正在进行路径约束数据压缩，从 %d 个点压缩到目标 %d 个点左右...\n', size(totalPathValidFlow, 1), targetSize);
        totalPathValidFlow = reduceDataPoints(totalPathValidFlow, targetSize);
    end
    
    % 保存结果到缓存文件
    save(cacheFileName, 'totalValidFlow', 'relationMatrix');
    fprintf('全部约束结果已保存到缓存文件%s\n', cacheFileName);
    
    % 保存只满足路径约束的数据到缓存文件
    save(pathConstraintCacheFileName, 'totalPathValidFlow', 'relationMatrix');
    fprintf('只满足路径约束的结果已保存到缓存文件%s\n', pathConstraintCacheFileName);
    
    % 设置随机选择的流量向量数量
    q = 30;  % 可以根据需要修改
    
    % 随机选择q个流量向量的索引
    if ~isempty(totalValidFlow)
        selectedIndices = randperm(size(totalValidFlow, 1), min(q, size(totalValidFlow, 1)));
        % 调用路径成本绘图函数
        plotPathCosts(totalValidFlow, relationMatrix, selectedIndices);
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

function [samplesMat, totalValidFlow, totalPathValidFlow] = processIteration(ii, n, rangeMin, rangeMax, bound, relationMatrix, totalValidFlow, totalPathValidFlow, zeta)
    % 处理单次迭代
    dimNum = ones(1,n)*ii;
    samples = generateSamples(n, rangeMin-bound, rangeMax+bound, dimNum);
    samplesMat = combineSamples(samples, n);
    
    % 计算目标函数和约束违反
    [path_err, money_err] = evaluateObjective(samplesMat, relationMatrix, zeta);
    
    % 满足所有约束的流量向量
    valid = path_err == 0 & money_err == 0;
    
    % 只满足路径约束的流量向量
    path_valid = path_err == 0;
    
    % 更新满足所有约束的结果
    if any(valid)
        if isempty(totalValidFlow)
            totalValidFlow = samplesMat(valid,:);
        else
            totalValidFlow = [totalValidFlow;samplesMat(valid,:)];
        end
    end
    
    % 更新只满足路径约束的结果
    if any(path_valid)
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

function [path_err, money_err] = evaluateObjective(f, M, zeta)
    % 评估约束违反
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
    
    % 不再计算和返回ff（目标函数值）
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

function [rangeMin, rangeMax, bound] = updateSearchRange(validSamples, ii)
    % 更新搜索范围
    rangeMin = min(validSamples);
    rangeMax = max(validSamples);
    bound = (rangeMax-rangeMin)/ii;
    bound(bound<20) = 20;
end

function M = pairsToMatrix(pairs)
    % 将有序对转换为矩阵
    M = zeros(max(pairs));
    for i = 1:size(pairs,1)
        M(pairs(i,1), pairs(i,2)) = 1;
    end
end

function reducedFlow = reduceDataPoints(totalValidFlow, targetSize)
    % 使用自适应网格采样减少数据点数量
    % 同时保持数据分布特性和重要点
    
    % 标准化数据到[0,1]区间
    minVals = min(totalValidFlow, [], 1);
    maxVals = max(totalValidFlow, [], 1);
    normFlow = (totalValidFlow - minVals) ./ (maxVals - minVals);
    
    % 计算每个维度上需要的网格数量，使总网格数接近目标数量
    gridSize = ceil(sqrt(targetSize));
    
    % 计算每个点所在的网格索引
    gridIndices = floor(normFlow * gridSize) + 1;
    
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
                meanPos = mean(normFlow(gridPoints,:), 1);
                % 找到最接近平均位置的点
                [~, minIdx] = min(sum((normFlow(gridPoints,:) - meanPos).^2, 2));
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
        paretoIndices = findParetoFrontier(totalValidFlow);
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
                    meanPos = mean(normFlow(gridPoints,:), 1);
                    % 找到最接近平均位置的点
                    [~, minIdx] = min(sum((normFlow(gridPoints,:) - meanPos).^2, 2));
                    selectedIndices(i) = gridPoints(minIdx);
                else
                    selectedIndices(i) = gridPoints(1);
                end
            end
        end
    end
    
    % 提取结果
    reducedFlow = totalValidFlow(selectedIndices, :);
    
    fprintf('数据点成功从 %d 减少到 %d\n', size(totalValidFlow, 1), length(selectedIndices));
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