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

    % 固定种子
    rng(4);
    
    % 随机选择流量
    if ~isempty(standardTotalValidFlow) && ~isempty(standardTotalPathValidFlow)
        % 从满足所有约束的流量中随机选择
        validIndices = randperm(size(standardTotalValidFlow, 1), min(numFlows, size(standardTotalValidFlow, 1)));
        selectedTotalValid = standardTotalValidFlow(validIndices, :);
        
        % 从只满足路径约束的流量中随机选择
        pathValidIndices = randperm(size(standardTotalPathValidFlow, 1), min(numFlows, size(standardTotalPathValidFlow, 1)));
        selectedPathValid = standardTotalPathValidFlow(pathValidIndices, :);
        
        % 组合选择的流量
        selectedFlows = [selectedTotalValid; selectedPathValid];
        numSelected = size(selectedFlows, 1);
        
        % 打印选择的流量方案
        fprintf('\n========= 选择的流量方案 =========\n');
        fprintf('满足全部约束的流量方案:\n');
        for i = 1:size(selectedTotalValid, 1)
            fprintf('方案 %d: [', i);
            fprintf(' %.2f', selectedTotalValid(i, :));
            fprintf(' ]\n');
        end
        
        fprintf('\n只满足路径约束的流量方案:\n');
        for i = 1:size(selectedPathValid, 1)
            fprintf('方案 %d: [', i);
            fprintf(' %.2f', selectedPathValid(i, :));
            fprintf(' ]\n');
        end
        fprintf('==================================\n\n');
        
        % 计算每个选择的流量对应的所有路径成本
        allPathCosts = cell(numSelected, 1);
        
        % 计算每个流量对应的所有路径成本
        for i = 1:numSelected
            flow = selectedFlows(i, :);
            pathCosts = calculateAllPathCosts(flow, fullRelationMatrix);
            allPathCosts{i} = pathCosts;
        end
        
        % 创建两组不同的颜色方案，与plotThreeRegionsPathCosts.m保持一致
        % 第一组: 蓝色系 (满足所有约束，对应plotThreeRegionsPathCosts中的区域3)
        % 第二组: 绿色系 (只满足路径约束，对应plotThreeRegionsPathCosts中的区域1)
        colorGroup1 = [
            0.0, 0.4, 0.8;  % 深蓝，与plotThreeRegionsPathCosts中的区域3一致
            0.0, 0.4, 0.8;  % 相同深蓝
            0.0, 0.4, 0.8   % 相同深蓝
        ];
        
        colorGroup2 = [
            0.0, 0.6, 0.3;  % 绿色，与plotThreeRegionsPathCosts中的区域1一致
            0.0, 0.6, 0.3;  % 相同绿色
            0.0, 0.6, 0.3   % 相同绿色
        ];
        
        % 确保颜色足够
        while size(colorGroup1, 1) < numFlows
            colorGroup1 = [colorGroup1; colorGroup1(end,:) * 0.9];
        end
        while size(colorGroup2, 1) < numFlows
            colorGroup2 = [colorGroup2; colorGroup2(end,:) * 0.9];
        end
        
        % 合并颜色
        colorVariations = [
            colorGroup1(1:numFlows, :); 
            colorGroup2(1:numFlows, :)
        ];
        
        % 创建图例标签
        legendLabels = cell(numSelected, 1);
        for i = 1:numFlows
            if i <= size(selectedTotalValid, 1)
                legendLabels{i} = sprintf('$BS_0^{\\zeta}$ %d', i); % 使用latex格式
            end
        end
        for i = 1:numFlows
            if i <= size(selectedPathValid, 1)
                legendLabels{i+numFlows} = sprintf('$S_0^{\\zeta}$ %d', i); % 使用latex格式
            end
        end
        
        % 绘制对比图
        plotPathComparison(allPathCosts, colorVariations, legendLabels, params, selectedFlows);
    else
        warning('输入的流量数据为空，无法绘图');
    end
end

function plotPathComparison(allPathCosts, colorVariations, legendLabels, params, selectedFlows)
    % 绘制路径成本对比图
    
    % 创建图形
    fig = figure('Name', '路径成本对比图', 'NumberTitle', 'off', ...
          'Position', params.FigurePosition);
    configureFigure(fig, params);

    % 创建不同的标记样式
    markerStyles = {'o', 's', 'd', '^', 'v', '>', 'p', 'h'};  % 为不同流量方案使用不同标记
    
    % 创建不同的线形样式
    lineStyles = {'-', '--', ':', '-.', '-', '--', ':', '-.'};  % 为不同流量方案使用不同线形
    
    % 创建隐藏的图例句柄
    q = size(colorVariations, 1);
    
    % 绘制所有点和连接线
    flowLegendHandles = [];
    flowLegendLabels = {};
    
    for i = 1:q
        costs = allPathCosts{i};
        if ~isempty(costs)
            % 确保markerStyles数组足够长
            markerIdx = mod(i-1, length(markerStyles)) + 1;
            % 确保lineStyles数组足够长
            lineIdx = mod(i-1, length(lineStyles)) + 1;
            
            % 先画散点 - 使用不同流量方案的不同标记
            h = scatter(costs(:,1), costs(:,2), 60, colorVariations(i,:), markerStyles{markerIdx}, ...
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