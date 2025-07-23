function plotComparisonPathCosts(totalValidFlow, totalPathValidFlow, relationMatrix, varargin)
    % 绘制同时满足双目标约束和仅满足路径约束的流量的路径成本对比图
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

    % 如果数据为空，直接返回
    if isempty(standardTotalValidFlow) || isempty(standardTotalPathValidFlow)
        warning('输入的流量数据为空，无法绘图');
        return;
    end

    % 固定种子
    % 如果三区域函数已存储选中的方案，则直接复用，以保证一致性
    storedValidFlows = getappdata(0, 'region3_flows');   % 满足全部约束且满足Tmax (区域3)
    storedPathFlows  = getappdata(0, 'region1_flows');   % 只满足路径约束 (区域1)

    % ----- 选择流量方案 -----
    selectedTotalValid = [];
    selectedPathValid  = [];

    % 辅助函数: 在数据集中匹配存储的流量
    findMatchingRows = @(data, targets) arrayfun(@(rowIdx) helperMatchRow(data, targets(rowIdx,:)), 1:size(targets,1))';

    % 内部辅助函数：匹配单行
    function idx = helperMatchRow(dataMat, targetRow)
        % 由于浮点和采样可能导致轻微差异，这里放宽容差至 0.1
        idxFound = find(all(abs(dataMat - targetRow) < 1e-1, 2), 1, 'first');
        if isempty(idxFound)
            idx = NaN;
        else
            idx = idxFound;
        end
    end

    % 初始化索引变量
    matchIdx = [];
    matchIdxPath = [];

    % 1) 尝试使用已保存的满足约束方案
    if ~isempty(storedValidFlows)
        matchIdx = findMatchingRows(standardTotalValidFlow, storedValidFlows);
        matchIdx = matchIdx(~isnan(matchIdx));
        if ~isempty(matchIdx)
            selectedTotalValid = standardTotalValidFlow(matchIdx(1:min(numFlows, length(matchIdx))), :);
        end
    end

    % 如果仍不足，补充前若干行
    if isempty(selectedTotalValid)
        selectedTotalValid = standardTotalValidFlow(1:min(numFlows, size(standardTotalValidFlow, 1)), :);
    elseif size(selectedTotalValid,1) < numFlows
        remaining = setdiff(1:size(standardTotalValidFlow,1), matchIdx, 'stable');
        needed = numFlows - size(selectedTotalValid,1);
        selectedTotalValid = [selectedTotalValid; standardTotalValidFlow(remaining(1:needed), :)];
    end

    % 2) 尝试使用已保存的路径约束方案
    if ~isempty(storedPathFlows)
        matchIdxPath = findMatchingRows(standardTotalPathValidFlow, storedPathFlows);
        matchIdxPath = matchIdxPath(~isnan(matchIdxPath));
        if ~isempty(matchIdxPath)
            selectedPathValid = standardTotalPathValidFlow(matchIdxPath(1:min(numFlows, length(matchIdxPath))), :);
        end
    end

    % 如果仍不足，补充前若干行
    if isempty(selectedPathValid)
        selectedPathValid = standardTotalPathValidFlow(1:min(numFlows, size(standardTotalPathValidFlow, 1)), :);
    elseif size(selectedPathValid,1) < numFlows
        remainingP = setdiff(1:size(standardTotalPathValidFlow,1), matchIdxPath, 'stable');
        neededP = numFlows - size(selectedPathValid,1);
        selectedPathValid = [selectedPathValid; standardTotalPathValidFlow(remainingP(1:neededP), :)];
    end

    % -----------------------------------

    % 组合选择的流量
    selectedFlows = [selectedTotalValid; selectedPathValid];
    numSelected = size(selectedFlows, 1);

    % 计算路径成本
    allPathCosts = cell(numSelected, 1);
    allTotalValidCosts = {};
    allPathValidCosts = {};
    for i = 1:size(selectedTotalValid, 1)
        flow = selectedTotalValid(i, :);
        costs = calculateAllPathCosts(flow, fullRelationMatrix);
        allPathCosts{i} = costs;
        allTotalValidCosts{end+1} = costs;
    end
    for i = 1:size(selectedPathValid, 1)
        flow = selectedPathValid(i, :);
        costs = calculateAllPathCosts(flow, fullRelationMatrix);
        allPathCosts{size(selectedTotalValid, 1) + i} = costs;
        allPathValidCosts{end+1} = costs;
    end

    % 打印所选流量方案
    printSelectedFlows(selectedTotalValid, selectedPathValid, allTotalValidCosts, allPathValidCosts);

    % -- 颜色保持不变，与之前一致 --
    colorGroup1 = repmat([0.0, 0.4, 0.8], numFlows, 1); % 满足全部约束 - 蓝色
    colorGroup2 = repmat([0.0, 0.6, 0.3], numFlows, 1); % 只满足路径约束 - 绿色
    colorVariations = [colorGroup1(1:size(selectedTotalValid,1), :); colorGroup2(1:size(selectedPathValid,1), :)];

    % 图例标签
    legendLabels = cell(numSelected,1);
    for i=1:size(selectedTotalValid,1)
        legendLabels{i} = sprintf('$BS_0^{\\zeta}$ %d', i);
    end
    for i=1:size(selectedPathValid,1)
        legendLabels{i+size(selectedTotalValid,1)} = sprintf('$S_0^{\\zeta}$ %d', i);
    end

    % ---- 调整绘图调用，传递满足全部约束方案数量以设置标记 ----
    plotPathComparison(allPathCosts, colorVariations, legendLabels, params, selectedFlows, size(selectedTotalValid,1));

    return; % 绘制完成后退出函数
    warning('输入的流量数据为空，无法绘图');
end

function printSelectedFlows(totalValid, pathValid, totalValidCosts, pathValidCosts)
    fprintf('\n========= 选择的流量方案 =========\n');
    fprintf('满足全部约束的流量方案:\n');
    for i = 1:size(totalValid, 1)
        fprintf('方案 %d: [', i); fprintf(' %.2f', totalValid(i, :)); fprintf(' ]\n');
        pathCosts = totalValidCosts{i};
        for j = 1:size(pathCosts, 1)
            fprintf('  路径 %d 成本: [时间: %.2f, 金钱: %.2f]\n', j, pathCosts(j, 1), pathCosts(j, 2));
        end
    end
    fprintf('\n只满足路径约束的流量方案:\n');
    for i = 1:size(pathValid, 1)
        fprintf('方案 %d: [', i); fprintf(' %.2f', pathValid(i, :)); fprintf(' ]\n');
        pathCosts = pathValidCosts{i};
        for j = 1:size(pathCosts, 1)
            fprintf('  路径 %d 成本: [时间: %.2f, 金钱: %.2f]\n', j, pathCosts(j, 1), pathCosts(j, 2));
        end
    end
    fprintf('==================================\n\n');
end

function plotPathComparison(allPathCosts, colorVariations, legendLabels, params, selectedFlows, validCount)
    % 绘制路径成本对比图
    
    % 创建图形
    fig = figure('Name', '路径成本对比图', 'NumberTitle', 'off', ...
          'Position', params.FigurePosition);
    configureFigure(fig, params);

    % 固定区域标记样式
    markerValid = 'd';  % 满足全部约束
    markerPathOnly = 'o'; % 只满足路径约束

    % 创建不同的线形样式
    lineStyles = {'-', '-', '-', '-', '-', '-', '-', '-'};
    
    % 创建隐藏的图例句柄
    q = size(colorVariations, 1);
    
    % 绘制所有点和连接线
    flowLegendHandles = [];
    flowLegendLabels = {};
    
    for i = 1:q
        costs = allPathCosts{i};
        if ~isempty(costs)
            % 根据区域选择标记样式
            if i <= validCount
                marker = markerValid;
            else
                marker = markerPathOnly;
            end
            % 确保lineStyles数组足够长
            lineIdx = mod(i-1, length(lineStyles)) + 1;
            
            % 先画散点 - 使用不同流量方案的不同标记
            h = scatter(costs(:,1), costs(:,2), 60, colorVariations(i,:), marker, ...
                'filled', 'MarkerEdgeColor', 'none');
            hold on;
            
            % 获取当前流量方案
            flow = selectedFlows(i, :);
            
            % 找出流量大于0的路径索引
            positiveFlowIdx = find(flow > 0);
            
            % 只对流量大于0的路径进行连线
            if ~isempty(positiveFlowIdx)
                % 提取流量大于0的路径成本
                positiveCosts = costs(positiveFlowIdx, :);
                
                % 按金钱成本排序
                [~, sortIdx] = sort(positiveCosts(:,2));
                sortedPositiveCosts = positiveCosts(sortIdx,:);
                
                % 画连接线 - 使用不同流量方案的不同线形
                plot(sortedPositiveCosts(:,1), sortedPositiveCosts(:,2), lineStyles{lineIdx}, 'Color', colorVariations(i,:), 'LineWidth', 1.5);
            end
            
            % 记录流量图例句柄
            flowLegendHandles = [flowLegendHandles; h];
            if i <= length(legendLabels)
                flowLegendLabels = [flowLegendLabels; legendLabels{i}];
            end
            
            % 添加路径编号标签
            for j = 1:size(costs, 1)
                text(costs(j,1), costs(j,2)-1, sprintf('P%d', j), ...
                    'FontSize', 8, 'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'bottom', 'Color', [0.3 0.3 0.3]);
            end
        end
    end
    
    % 配置坐标轴和标签
    xlabel('Time Cost', 'FontName', params.FontName, 'FontSize', params.FontSize);
    ylabel('Money Cost', 'FontName', params.FontName, 'FontSize', params.FontSize);
    % title('路径成本分析 (所有路径点按流量方案连线)', 'FontName', params.FontName, 'FontSize', params.FontSize+1);
    
    if params.ShowGrid
        grid on;
        grid minor;
    end
    
    % 添加图例
    legend(flowLegendHandles, flowLegendLabels, ...
        'Location', 'northeast', ...
        'FontName', params.FontName, 'FontSize', params.FontSize-1, ...
        'EdgeColor', [0.7, 0.7, 0.7], ...
        'Box', 'on', ...
        'Interpreter', 'latex');
    
    % 保存图形
    saveFigure(fig);
end

function boundary = calculateBoundary(allPathCosts)
    % 计算边界
    % 合并所有路径成本
    allCosts = [];
    for i = 1:length(allPathCosts)
        if ~isempty(allPathCosts{i})
            allCosts = [allCosts; allPathCosts{i}];
        end
    end
    
    if isempty(allCosts)
        % 如果没有成本数据，返回默认边界
        boundary.leftX = [0; 10];
        boundary.leftY = [10; 0];
        boundary.rightX = [20; 10];
        boundary.rightY = [0; 10];
        return;
    end
    
    % 确定数据边界
    minX = min(allCosts(:,1)) * 0.8;
    maxX = max(allCosts(:,1)) * 1.2;
    minY = min(allCosts(:,2)) * 0.8;
    maxY = max(allCosts(:,2)) * 1.2;
    
    % 创建边界框
    boundary.leftX = [minX; maxX];
    boundary.leftY = [maxY; minY];
    boundary.rightX = [maxX; minX];
    boundary.rightY = [minY; maxY];
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
        fileName = [outputDir 'path_costs_comparison_' datestr(now, 'yyyymmdd_HHMMSS') '.pdf'];
    else
        % 确保zeta和subset_index是数值
        if ~isnumeric(zeta)
            zeta = str2double(char(zeta));
        end
        if ~isnumeric(subset_index)
            subset_index = str2double(char(subset_index));
        end
        
        % 使用zeta和subset_index生成名称
        fileName = sprintf('%spath_costs_comparison_zeta%d_subset%d.pdf', outputDir, zeta, subset_index);
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