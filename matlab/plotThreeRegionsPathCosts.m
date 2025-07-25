function plotThreeRegionsPathCosts(totalValidFlow, totalPathValidFlow, relationMatrix, varargin)
    % 绘制三个区域的路径成本对比图：
    % 1. 满足全部约束的流量方案 (BS_0^{zeta})
    % 2. 只满足路径约束的流量方案 (S_0^{zeta})
    % 3. 满足T_max上限约束的流量方案
    %
    % 输入参数:
    %   totalValidFlow     - 满足所有约束的流量数据矩阵
    %   totalPathValidFlow - 只满足路径约束的流量数据矩阵
    %   relationMatrix     - 关系矩阵
    %   varargin           - 可选参数:
    %                         'SavePath': 图像保存路径 (默认: 'results/')
    %                         'FigurePosition': 图像位置和大小 (默认: [100, 100, 800, 600])
    %                         'FontName': 字体名称 (默认: 'Arial')
    %                         'FontSize': 字体大小 (默认: 10)
    %                         'ShowGrid': 是否显示网格 (默认: true)
    %                         'NumFlows': 每种类型选择的流量数量 (默认: 2)
    
    % 解析可选参数
    p = inputParser;
    defaultSavePath = 'results/';
    defaultFigPosition = [100, 100, 800, 600];
    defaultFontName = 'Arial';
    defaultFontSize = 10;
    defaultShowGrid = true;
    defaultNumFlows = 2; % 默认每种类型选择2个流量
    
    addParameter(p, 'SavePath', defaultSavePath);
    addParameter(p, 'FigurePosition', defaultFigPosition);
    addParameter(p, 'FontName', defaultFontName);
    addParameter(p, 'FontSize', defaultFontSize);
    addParameter(p, 'ShowGrid', defaultShowGrid);
    addParameter(p, 'NumFlows', defaultNumFlows);
    
    parse(p, varargin{:});
    params = p.Results;
    
    % 在函数开头添加，从varargin中提取zeta和subset_index
    zeta = [];
    subset_index = [];
    
    % 查找varargin中是否有zeta和subset_index参数
    for i = 1:2:length(varargin)
        if strcmpi(varargin{i}, 'zeta')
            zeta = varargin{i+1};
        elseif strcmpi(varargin{i}, 'subset_index')
            subset_index = varargin{i+1};
        end
    end
    
    % 如果参数中没有，尝试从应用数据获取
    if isempty(zeta)
        zeta = getappdata(0, 'current_zeta');
    end
    
    if isempty(subset_index)
        subset_index = getappdata(0, 'current_subset_index');
    end
    
    % 确保结果目录存在
    if ~exist(params.SavePath, 'dir')
        mkdir(params.SavePath);
    end
    
    % 确定选择的流量数量
    numFlows = params.NumFlows;
    
    % 标准化流量数据和关系矩阵
    % 创建标准的完整关系矩阵（6条路径×8个链路）
    fullOrderedPair = [1,1;2,2;3,3;5,3;4,4;6,4;5,5;6,6;3,7;6,7;4,8;5,8];
    fullRelationMatrix = createFullRelationMatrix(fullOrderedPair);
    
    % 获取当前数据中的路径索引映射
    [pathIndices, fullToCurrentMap] = getPathIndicesMap(relationMatrix);
    
    % 将流量数据标准化为n×6的矩阵
    standardTotalValidFlow = standardizeFlowData(totalValidFlow, pathIndices, 6);
    standardTotalPathValidFlow = standardizeFlowData(totalPathValidFlow, pathIndices, 6);

    % 固定种子
    rng(4);
    
    % 随机选择流量
    if ~isempty(standardTotalValidFlow) && ~isempty(standardTotalPathValidFlow)
        % 计算所有满足全部约束的流量的路径成本（用于确定边界和T_max上限）
        allTimeCosts = [];
        allMoneyCosts = [];
        
        % 对所有standardTotalValidFlow进行计算
        for i = 1:size(standardTotalValidFlow, 1)
            flow = standardTotalValidFlow(i, :);
            pathCosts = calculateAllPathCosts(flow, fullRelationMatrix);
            
            % 只收集非零流量路径的成本
            nonZeroFlowIdx = find(flow > 0);
            if ~isempty(nonZeroFlowIdx)
                % 提取非零流量路径的成本
                nonZeroCosts = pathCosts(nonZeroFlowIdx, :);
                % 排除金钱成本为0的路径
                nonZeroMoneyIdx = nonZeroCosts(:,2) > 0;
                if any(nonZeroMoneyIdx)
                    allTimeCosts = [allTimeCosts; nonZeroCosts(nonZeroMoneyIdx, 1)];
                    allMoneyCosts = [allMoneyCosts; nonZeroCosts(nonZeroMoneyIdx, 2)];
                end
            end
        end
        
        % 计算可行域边界（仅基于满足全部约束的流量）
        boundary = calculateFeasibleRegionBoundary(allTimeCosts, allMoneyCosts);
        
        % 计算T_max上限（使用右边界最大时间）
        upperLimitX = (boundary.rightX+boundary.leftX)/2;   % 直接使用最大时间作为上限
        upperLimitY = boundary.rightY;   % 与计算方法保持一致
        % 确保 upperLimitX 与 upperLimitY 长度一致
        if length(upperLimitX) ~= length(upperLimitY)
            minLen = min(length(upperLimitX), length(upperLimitY));
            upperLimitX = upperLimitX(1:minLen);
            upperLimitY = upperLimitY(1:minLen);
        end
        
        % 对满足全部约束的所有流量进行 T_max 检查，划分为可行(区域3)和不可行(区域2)
        tmaxFeasibleAC = [];
        tmaxInfeasibleAC = [];
        tmaxFeasibleCostsAC = {};
        tmaxInfeasibleCostsAC = {};
        for i = 1:size(standardTotalValidFlow, 1)
            flow = standardTotalValidFlow(i, :);
            pathCosts = calculateAllPathCosts(flow, fullRelationMatrix);
            [isFeasible,~] = checkPathFeasibility(pathCosts, upperLimitX, upperLimitY);
            if isFeasible
                tmaxFeasibleAC = [tmaxFeasibleAC; flow];
                tmaxFeasibleCostsAC{end+1} = pathCosts; %#ok<AGROW>
            else
                tmaxInfeasibleAC = [tmaxInfeasibleAC; flow];
                tmaxInfeasibleCostsAC{end+1} = pathCosts; %#ok<AGROW>
            end
        end

        % 从区域3（满足全部约束且满足T_max）随机选择
        feasibleIndicesSel = randperm(size(tmaxFeasibleAC,1), min(numFlows, size(tmaxFeasibleAC,1)));
        selectedRegion3 = tmaxFeasibleAC(feasibleIndicesSel,:);
        allRegion3Costs = tmaxFeasibleCostsAC(feasibleIndicesSel);

        % 从区域2（满足全部约束但违反T_max）随机选择
        infeasibleIndicesSel = randperm(size(tmaxInfeasibleAC,1), min(numFlows, size(tmaxInfeasibleAC,1)));
        selectedRegion2 = tmaxInfeasibleAC(infeasibleIndicesSel,:);
        allRegion2Costs = tmaxInfeasibleCostsAC(infeasibleIndicesSel);

        % -------- 区域1：只满足路径约束 ----------
        pathOnlyValidFlow = filterOutOverlappingFlows(standardTotalPathValidFlow, standardTotalValidFlow);
        if ~isempty(pathOnlyValidFlow)
            region1Idx = randperm(size(pathOnlyValidFlow,1), min(numFlows, size(pathOnlyValidFlow,1)));
            selectedRegion1 = pathOnlyValidFlow(region1Idx,:);
            allRegion1Costs = cell(size(selectedRegion1,1),1);
            for i = 1:size(selectedRegion1,1)
                allRegion1Costs{i} = calculateAllPathCosts(selectedRegion1(i,:), fullRelationMatrix);
            end
        else
            selectedRegion1 = [];
            allRegion1Costs = {};
        end

        % 重新组合所有区域的数据 (确保顺序: 区域1,2,3)
        selectedFlows = [selectedRegion1; selectedRegion2; selectedRegion3];
        % 将各区域成本 cell 数组转换为列向量，避免维度不一致
        allRegion1Costs = allRegion1Costs(:);
        allRegion2Costs = allRegion2Costs(:);
        allRegion3Costs = allRegion3Costs(:);
        allPathCosts   = [allRegion1Costs; allRegion2Costs; allRegion3Costs];

        % --- 新增: 将选中的流量方案存入全局 appdata，供其它绘图函数复用 ---
        setappdata(0, 'region1_flows', selectedRegion1);
        setappdata(0, 'region2_flows', selectedRegion2);
        setappdata(0, 'region3_flows', selectedRegion3);
        % --------------------------------------------------------------

        % 重新分配颜色与图例（一次性精确匹配）
        colorRegion1 = [0.0, 0.6, 0.3];   % 绿色系
        colorRegion2 = [0.85, 0.4, 0.0];  % 橙色系
        colorRegion3 = [0.0, 0.4, 0.8];   % 蓝色系

        % 为每个区域生成颜色矩阵，确保与对应流量数完全一致
        colorVariations = [
            repmat(colorRegion1, size(selectedRegion1,1), 1);
            repmat(colorRegion2, size(selectedRegion2,1), 1);
            repmat(colorRegion3, size(selectedRegion3,1), 1)
        ];

        legendLabels = {};
        for i=1:size(selectedRegion1,1)
            legendLabels{end+1} = sprintf('$S_0^{\\varepsilon}$ %d',i); %#ok<AGROW>
        end
        for i=1:size(selectedRegion2,1)
            legendLabels{end+1} = sprintf('$BS_0^{\\varepsilon}$ %d',i); %#ok<AGROW>
        end
        for i=1:size(selectedRegion3,1)
            legendLabels{end+1} = sprintf('$RS_0^{\\varepsilon}$ %d',i); %#ok<AGROW>
        end

        % 打印信息
        printSelectedFlows_New(selectedRegion1, selectedRegion2, selectedRegion3, allRegion1Costs, allRegion2Costs, allRegion3Costs);

        % ---- 调整函数调用，传递区域流量数量以便统一标记样式 ----
        plotThreeRegionsComparison(allPathCosts, colorVariations, legendLabels, params, selectedFlows, upperLimitX, upperLimitY, boundary, size(selectedRegion1,1), size(selectedRegion2,1));
        % --------------------------------------------------------------

        % 先存储到应用数据中
        setappdata(0, 'current_zeta', zeta);
        setappdata(0, 'current_subset_index', subset_index);

        % 从应用数据中重新获取值，确保变量已定义
        zeta = getappdata(0, 'current_zeta');
        subset_index = getappdata(0, 'current_subset_index');

        % 现在调用saveTmaxData时，zeta和subset_index已定义
        saveTmaxData(upperLimitX, upperLimitY, zeta, subset_index);
    else
        warning('输入的流量数据为空，无法绘图');
    end
end

function plotThreeRegionsComparison(allPathCosts, colorVariations, legendLabels, params, selectedFlows, upperLimitX, upperLimitY, boundary, region1Count, region2Count)
    % 绘制三个区域的路径成本对比图 (区域1: S0^zeta, 区域2: BS0^zeta\Tmax, 区域3: BS0^zeta∩Tmax)

    % 创建图形
    fig = figure('Name', '三区域路径成本对比图', 'NumberTitle', 'off', ...
          'Position', params.FigurePosition);
    configureFigure(fig, params);

    % 区域固定标记样式
    markerRegion1 = 'o';   % 圆形
    markerRegion2 = 's';   % 方形
    markerRegion3 = 'd';   % 菱形

    % 创建不同的线形样式（仍保持多样）
    lineStyles = {'-', '-', '-', '-', '-', '-', '-', '-'};

    % 绘制所有点和连接线
    flowLegendHandles = [];
    flowLegendLabels  = {};

    % 绘制边界和 Tmax
    % h_boundary = plot([boundary.leftX; boundary.rightX], [boundary.leftY; boundary.rightY], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
    hold on;
    h_tmax = plot(upperLimitX, upperLimitY, '-.', 'Color', [0.8 0.2 0.2], 'LineWidth', 2);

    totalFlows = numel(allPathCosts);
    for i = 1:totalFlows
        costs = allPathCosts{i};
        if ~isempty(costs)
            % 选择区域对应的标记
            if i <= region1Count
                marker = markerRegion1;
            elseif i <= region1Count + region2Count
                marker = markerRegion2;
            else
                marker = markerRegion3;
            end
            % 选择线形
            lineIdx = mod(i-1, length(lineStyles)) + 1;
            % 画散点
            h = scatter(costs(:,1), costs(:,2), 60, colorVariations(i,:), marker, 'filled', 'MarkerEdgeColor', 'none');
            hold on;

            % 获取当前流量方案
            flow = selectedFlows(i, :);
            positiveIdx = find(flow > 0);
            if ~isempty(positiveIdx)
                positiveCosts = costs(positiveIdx, :);
                [~, sortIdx] = sort(positiveCosts(:,2));
                sortedPositiveCosts = positiveCosts(sortIdx, :);
                plot(sortedPositiveCosts(:,1), sortedPositiveCosts(:,2), lineStyles{lineIdx}, 'Color', colorVariations(i,:), 'LineWidth', 1.5);
            end

            flowLegendHandles = [flowLegendHandles; h];
            if i <= length(legendLabels)
                flowLegendLabels = [flowLegendLabels; legendLabels{i}];
            end

            % 添加路径编号
            for j = 1:size(costs,1)
                text(costs(j,1), costs(j,2)-1, sprintf('P%d', j), 'FontSize', 8, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Color', [0.3 0.3 0.3]);
            end
        end
    end

    % 坐标轴与网格
    xlabel('Time Cost', 'FontName', params.FontName, 'FontSize', params.FontSize);
    ylabel('Money Cost', 'FontName', params.FontName, 'FontSize', params.FontSize);
    if params.ShowGrid
        grid on; grid minor;
    end

    % 图例 (添加 Tmax)
    flowLegendHandles = [flowLegendHandles; h_tmax];
    flowLegendLabels  = [flowLegendLabels; '$T_{max}$'];
    legend(flowLegendHandles, flowLegendLabels, 'Location', 'northeast', 'FontName', params.FontName, 'FontSize', params.FontSize-1, 'EdgeColor', [0.7 0.7 0.7], 'Box', 'on', 'Interpreter', 'latex');

    % 保存并关闭
    saveFigure(fig);
end

function filteredFlow = filterOutOverlappingFlows(flowA, flowB)
    % 从flowA中过滤掉与flowB重叠的流量向量
    % 将矩阵四舍五入到小数点后3位，以避免浮点误差
    roundedFlowA = round(flowA * 1000) / 1000;
    roundedFlowB = round(flowB * 1000) / 1000;
    
    % 找出不重叠的索引
    nonOverlappingIndices = [];
    for i = 1:size(roundedFlowA, 1)
        isOverlap = false;
        for j = 1:size(roundedFlowB, 1)
            if all(abs(roundedFlowA(i,:) - roundedFlowB(j,:)) < 1e-6)
                isOverlap = true;
                break;
            end
        end
        if ~isOverlap
            nonOverlappingIndices = [nonOverlappingIndices; i];
        end
    end
    
    % 返回不重叠的流量向量
    filteredFlow = flowA(nonOverlappingIndices, :);
end

function boundary = calculateFeasibleRegionBoundary(allPathTimeCosts, allPathMoneyCosts)
    % 计算可行域边界，基于路径成本
    
    % 如果没有成本数据，返回默认边界
    if isempty(allPathTimeCosts) || isempty(allPathMoneyCosts)
        boundary.leftX = [0; 10];
        boundary.leftY = [10; 0];
        boundary.rightX = [20; 10];
        boundary.rightY = [0; 10];
        return;
    end
    
    % 找到唯一的金钱成本值
    uniqueMoneyValues = unique(allPathMoneyCosts);
    leftBoundaryX = [];
    leftBoundaryY = [];
    rightBoundaryX = [];
    rightBoundaryY = [];
    
    % 对于每个金钱成本值，找到最小和最大时间成本
    for i = 1:length(uniqueMoneyValues)
        currMoney = uniqueMoneyValues(i);
        sameMoneyIdx = abs(allPathMoneyCosts - currMoney) < 0.001;
        
        if sum(sameMoneyIdx) > 0
            timesForMoney = allPathTimeCosts(sameMoneyIdx);
            minTimeForMoney = min(timesForMoney);
            maxTimeForMoney = max(timesForMoney);
            
            % 添加到边界数组
            leftBoundaryX = [leftBoundaryX; minTimeForMoney];
            leftBoundaryY = [leftBoundaryY; currMoney];
            rightBoundaryX = [rightBoundaryX; maxTimeForMoney];
            rightBoundaryY = [rightBoundaryY; currMoney];
        end
    end
    
    % 按金钱成本排序边界点
    [leftBoundaryY, sortIdx] = sort(leftBoundaryY);
    leftBoundaryX = leftBoundaryX(sortIdx);
    
    [rightBoundaryY, sortIdx] = sort(rightBoundaryY);
    rightBoundaryX = rightBoundaryX(sortIdx);
    
    % 创建边界结构
    boundary = struct(...
        'leftX', leftBoundaryX, ...
        'leftY', leftBoundaryY, ...
        'rightX', rightBoundaryX, ...
        'rightY', rightBoundaryY ...
    );
end

function upperLimitX = calculateUpperLimit(leftX, rightX)
    % 计算上限曲线
    % 使用边界的中点
    midX = (leftX + rightX) / 2;
    upperLimitX = midX;
end

function [isFeasible, infeasibleIndices] = checkPathFeasibility(pathCosts, upperLimitX, upperLimitY)
    % 检查路径成本是否满足T_max上限约束
    % 
    % 输入:
    %   pathCosts - 路径成本 [时间成本, 金钱成本]
    %   upperLimitX - T_max上限曲线的x坐标
    %   upperLimitY - T_max上限曲线的y坐标
    % 
    % 输出:
    %   isFeasible - 是否所有路径都满足约束
    %   infeasibleIndices - 不满足约束的路径索引
    
    infeasibleIndices = [];
    
    % 对每个路径成本点进行检查
    for i = 1:size(pathCosts, 1)
        timeCost = pathCosts(i, 1);
        moneyCost = pathCosts(i, 2);
        
        % 找到对应金钱成本的T_max上限值
        idx = find(upperLimitY==moneyCost);
        if ~isempty(idx)
            % 如果时间成本超过上限，则路径不可行
            if timeCost > upperLimitX(idx)
                infeasibleIndices = [infeasibleIndices; i];
            end
        end
    end
    
    % 如果存在不可行路径，则整个方案不可行
    isFeasible = isempty(infeasibleIndices);
end

function configureFigure(fig, params)
    % 配置图形
    set(fig, 'Color', 'white');
    set(gca, 'FontName', params.FontName, 'FontSize', params.FontSize);
    box on;
end

function saveFigure(fig)
    % 简化的保存图形函数
    
    % 获取zeta和subset_index值
    zeta = getappdata(0, 'current_zeta');
    subset_index = getappdata(0, 'current_subset_index');
    
    % 简化处理：直接使用固定路径
    outputDir = 'results/pdf_outputs/';
    
    % 确保输出目录存在
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end
    
    % 生成文件名
    if isempty(zeta) || isempty(subset_index)
        % 使用时间戳作为默认名称
        fileName = [outputDir 'three_regions_comparison_' datestr(now, 'yyyymmdd_HHMMSS') '.pdf'];
    else
        % 确保zeta和subset_index是数值
        if ~isnumeric(zeta)
            zeta = str2double(char(zeta));
        end
        if ~isnumeric(subset_index)
            subset_index = str2double(char(subset_index));
        end
        
        % 使用zeta和subset_index生成名称
        fileName = sprintf('%sthree_regions_comparison_zeta%d_subset%d.pdf', outputDir, zeta, subset_index);
    end
    
    % 获取图形的实际大小
    figPos = get(fig, 'Position');
    figWidth = figPos(3);
    figHeight = figPos(4);
    
    % 设置纸张大小与图形大小一致
    set(fig, 'PaperPositionMode', 'manual');
    set(fig, 'PaperUnits', 'points');
    set(fig, 'PaperSize', [figWidth figHeight]);
    set(fig, 'PaperPosition', [0 0 figWidth figHeight]);

    % 保存图形为PDF
    print(fig, fileName, '-dpdf', '-r300');
    fprintf('图形已保存为: %s\n', fileName);
    
    % 关闭绘图保持
    hold off;
end

function pathCosts = calculateAllPathCosts(flow, relationMatrix)
    % 计算所有路径成本（包括流量为0的路径）
    
    % 获取相关常量
    money = [20, 15, 1, 0, 0, 0, 0, 1];
    freeFlowTime = [18, 22.5, 12, 24, 2.4, 6, 24, 12];
    maxCapacity = [3600, 3600, 1800, 1800, 1800, 1800, 1800, 1800];
    
    % 计算链路流量
    x = flow * relationMatrix;
    
    % 计算实际行驶时间
    index = find(sum(relationMatrix) ~= 0);
    time = freeFlowTime(index) .* (1 + 0.15 .* (x(index) ./ maxCapacity(index)).^4);
    
    % 计算每条路径的时间成本
    pathTime = time * relationMatrix(:, index)';
    RT = pathTime + 15 * (1 - exp(-0.02 * pathTime));
    
    % 计算每条路径的货币成本
    pathMoney = money * relationMatrix';
    
    % 组合时间和货币成本
    pathCosts = [RT', pathMoney'];
end

function standardFlow = standardizeFlowData(flowData, pathIndices, totalPaths)
    % 将流量数据标准化为n×totalPaths格式
    
    if isempty(flowData)
        standardFlow = [];
        return;
    end
    
    n = size(flowData, 1);
    standardFlow = zeros(n, totalPaths);
    
    % 将现有流量数据填入对应位置
    for i = 1:length(pathIndices)
        standardFlow(:, pathIndices(i)) = flowData(:, i);
    end
end

function [pathIndices, fullToCurrentMap] = getPathIndicesMap(relationMatrix)
    % 确定当前关系矩阵中的路径对应全局路径的哪些索引
    % subset_index=0时: 1,2,5路径 -> 索引[1,2,5]
    % subset_index=1时: 1,2,3,4,5路径 -> 索引[1,2,3,4,5]
    
    % 检查矩阵行数来确定当前使用的子集
    n = size(relationMatrix, 1);
    
    switch n
        case 3  % subset_index = 0: 路径1,2,5
            pathIndices = [1, 2, 5];
            fullToCurrentMap = zeros(6, 1);
            fullToCurrentMap([1, 2, 5]) = [1, 2, 3];
        case 5  % subset_index = 1: 路径1,2,3,4,5
            pathIndices = [1, 2, 3, 4, 5];
            fullToCurrentMap = zeros(6, 1);
            fullToCurrentMap([1, 2, 3, 4, 5]) = [1, 2, 3, 4, 5];
        case 6  % subset_index = 2: 所有路径
            pathIndices = 1:6;
            fullToCurrentMap = (1:6)';
        otherwise
            error('未知的关系矩阵维度: %d', n);
    end
end

function M = createFullRelationMatrix(fullOrderedPair)
    % 创建完整的关系矩阵 (6条路径 × 8个链路)
    M = zeros(6, 8);  % 6条路径, 8个链路
    
    % 根据用户提供的信息，创建正确的路径-链路关系
    % 路径1: 链路1
    M(1, 1) = 1;
    
    % 路径2: 链路2
    M(2, 2) = 1;
    
    % 路径3: 链路3,7
    M(3, 3) = 1;
    M(3, 7) = 1;
    
    % 路径4: 链路4,8
    M(4, 4) = 1;
    M(4, 8) = 1;
    
    % 路径5: 链路3,5,8
    M(5, 3) = 1;
    M(5, 5) = 1;
    M(5, 8) = 1;
    
    % 路径6: 链路4,6,7
    M(6, 4) = 1;
    M(6, 6) = 1;
    M(6, 7) = 1;
end 

function printSelectedFlows_New(region1, region2, region3, costs1, costs2, costs3)
    fprintf('\n========= 选择的流量方案 =========\n');
    fprintf('区域1: 只满足路径约束 (S_0^{zeta})\n');
    for i=1:size(region1,1)
        fprintf('方案 %d: [',i); fprintf(' %.2f',region1(i,:)); fprintf(' ]\n');
        pathCosts = costs1{i};
        for j = 1:size(pathCosts, 1)
            fprintf('  路径 %d 成本: [时间: %.2f, 金钱: %.2f]\n', j, pathCosts(j, 1), pathCosts(j, 2));
        end
    end
    fprintf('\n区域2: 满足全部约束但违反T_{max} (BS_0^{zeta} \\ T_{max})\n');
    for i=1:size(region2,1)
        fprintf('方案 %d: [',i); fprintf(' %.2f',region2(i,:)); fprintf(' ]\n');
        pathCosts = costs2{i};
        for j = 1:size(pathCosts, 1)
            fprintf('  路径 %d 成本: [时间: %.2f, 金钱: %.2f]\n', j, pathCosts(j, 1), pathCosts(j, 2));
        end
    end
    fprintf('\n区域3: 满足全部约束且满足T_{max} (BS_0^{zeta} ∩ T_{max})\n');
    for i=1:size(region3,1)
        fprintf('方案 %d: [',i); fprintf(' %.2f',region3(i,:)); fprintf(' ]\n');
        pathCosts = costs3{i};
        for j = 1:size(pathCosts, 1)
            fprintf('  路径 %d 成本: [时间: %.2f, 金钱: %.2f]\n', j, pathCosts(j, 1), pathCosts(j, 2));
        end
    end
    fprintf('==================================\n\n');
end 

function saveTmaxData(limitX, limitY, zeta, subset_index)
    % 如果 zeta 和 subset_index 为空，使用默认值
    if isempty(zeta) || isempty(subset_index)
        zeta = 0;
        subset_index = 0;
    end
    
    % 确保是数值
    if ~isnumeric(zeta)
        zeta = str2double(char(zeta));
    end
    if ~isnumeric(subset_index)
        subset_index = str2double(char(subset_index));
    end
    
    % 创建 T_max 数据表 - 使用英文字段名
    tmax_data = table();
    tmax_data.zeta = repmat(zeta, length(limitX), 1);
    tmax_data.type = repmat({'T_max'}, length(limitX), 1);
    tmax_data.time_cost = limitX;
    tmax_data.money_cost = limitY;
    
    % 保存文件路径
    csvFile = sprintf('results/tmax_data_zeta%d_subset%d.csv', zeta, subset_index);
    
    % 设置列名
    tmax_data.Properties.VariableNames = {'zeta值', '数据类型', '时间成本', '金钱成本'};
    
    % 确保目录存在
    if ~exist('results', 'dir')
        mkdir('results');
    end
    
    % 写入CSV
    writetable(tmax_data, csvFile);
    fprintf('T_max 数据已保存至: %s\n', csvFile);
    
    % 同时将结果附加到 result.csv，确保格式兼容
    resultFile = 'results/result.csv';
    if exist(resultFile, 'file')
        % 读取现有数据
        existing_data = readtable(resultFile);
        
        % 准备T_max数据，与结果文件格式保持一致
        point_count = length(limitX);
        tmax_for_result = table();
        tmax_for_result.zeta = repmat(zeta, point_count, 1);
        tmax_for_result.region = repmat(4, point_count, 1); 
        tmax_for_result.desc = repmat({'T_max上限曲线'}, point_count, 1);
        tmax_for_result.scheme = repmat(1, point_count, 1);
        tmax_for_result.path = (1:point_count)';
        tmax_for_result.flow = zeros(point_count, 1);
        tmax_for_result.time = limitX;
        tmax_for_result.money = limitY;
        
        % 设置列名与result.csv匹配
        tmax_for_result.Properties.VariableNames = {'zeta值', '区域', '区域描述', '方案编号', '路径编号', '流量分配', '时间成本', '金钱成本'};
        
        % 合并并重新保存
        updated_data = [existing_data; tmax_for_result];
        writetable(updated_data, resultFile);
        fprintf('T_max 数据已附加至结果文件: %s\n', resultFile);
    else
        fprintf('警告: 结果文件 %s 不存在，无法附加 T_max 数据\n', resultFile);
    end
end 