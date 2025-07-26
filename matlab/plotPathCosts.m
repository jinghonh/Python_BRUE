function plotPathCosts(totalValidFlow, relationMatrix, selectedIndices, varargin)
    % Plot time and money costs for each path with scientific style
    % 
    % Input parameters:
    %   totalValidFlow  - All feasible flow matrix (M x n)
    %   relationMatrix  - Relation matrix
    %   selectedIndices - Selected flow vector indices
    %   varargin        - Optional name-value pairs:
    %                      'FigurePosition': [100, 100, 800, 600] (default)
    %                      'SavePath': 'results/' (default)
    %                      'FontName': 'Arial' (default)
    %                      'FontSize': 10 (default)
    %                      'ShowGrid': true (default)
    
    % Parse optional parameters
    params = parseInputParameters(varargin{:});
    
    % Define constants
    money = [20, 15, 1, 0, 0, 0, 0, 1];
    freeFlowTime = [18,22.5,12,24,2.4,6,24,12];
    maxCapacity = [3600,3600,1800,1800,1800,1800,1800,1800];
    money = money * relationMatrix';  % 1 x n
    
    % Process data for plotting
    [allPathCosts, allPathTimeCosts, allPathMoneyCosts, colorVariations, legendLabels] = ...
        preparePathCostsData(totalValidFlow, selectedIndices, relationMatrix, money, freeFlowTime, maxCapacity);
    
    % Calculate boundary for feasible region
    boundary = calculateFeasibleRegionBoundary(allPathTimeCosts, allPathMoneyCosts);
    
    % Plot time-money cost relationship
    plotTimeMoneyCostRelationship(allPathCosts, boundary, colorVariations, legendLabels, params);
    
    % Plot costs with upper limit
    plotPathCostsWithUpperLimit(allPathCosts, boundary, colorVariations, params);
end

function params = parseInputParameters(varargin)
    % Parse optional input parameters with defaults
    p = inputParser;
    
    % Define default values
    defaultFigPosition = [100, 100, 800, 600];
    defaultSavePath = 'results/';
    defaultFontName = 'Arial';
    defaultFontSize = 10;
    defaultShowGrid = true;
    
    % Add parameters
    addParameter(p, 'FigurePosition', defaultFigPosition);
    addParameter(p, 'SavePath', defaultSavePath);
    addParameter(p, 'FontName', defaultFontName);
    addParameter(p, 'FontSize', defaultFontSize);
    addParameter(p, 'ShowGrid', defaultShowGrid);
    
    % Parse inputs
    parse(p, varargin{:});
    
    % Get results
    params = p.Results;
end

function [allPathCosts, allPathTimeCosts, allPathMoneyCosts, colorVariations, legendLabels] = ...
        preparePathCostsData(totalValidFlow, selectedIndices, relationMatrix, money, freeFlowTime, maxCapacity)
    % Process and prepare path costs data for plotting
    q = length(selectedIndices);
    legendLabels = cell(q, 1);
    
    % Store all path time and money costs
    allPathTimeCosts = [];
    allPathMoneyCosts = [];
    allFlowVectorIndices = [];
    allPathCosts = cell(q, 1);
    
    % Calculate costs for each selected flow vector
    for i = 1:q
        % Get current flow vector
        currentFlow = totalValidFlow(selectedIndices(i), :);
        
        % Calculate time cost
        x = currentFlow * relationMatrix;
        RT = calculateRealTime(x, relationMatrix, freeFlowTime, maxCapacity);
        
        % Get time and money costs for each path
        pathTimeCosts = RT(1, :);  % 1 x n
        pathMoneyCosts = money;    % 1 x n
        
        % Remove paths with zero cost
        validPaths = pathTimeCosts > 0;
        pathMoneyCosts = pathMoneyCosts(validPaths);
        pathTimeCosts = pathTimeCosts(validPaths);
        
        % Store all path data
        allPathTimeCosts = [allPathTimeCosts; pathTimeCosts'];
        allPathMoneyCosts = [allPathMoneyCosts; pathMoneyCosts'];
        allFlowVectorIndices = [allFlowVectorIndices; i*ones(length(pathTimeCosts), 1)];
        
        % Combine and sort costs by money cost
        costs = [pathTimeCosts', pathMoneyCosts'];
        [~, sortIdx] = sort(pathMoneyCosts);
        costs = costs(sortIdx, :);
        allPathCosts{i} = costs;
        
        % Create legend labels
        legendLabels{i} = sprintf('Flow Vector %d', i);
    end
    
    % Create elegant color scheme
    colorVariations = createColorScheme(q);
end

function colorVariations = createColorScheme(numColors)
    % Create an elegant color scheme based on the number of flow vectors
    baseColor = [0.2 0.4 0.8]; % Base blue color
    colorVariations = zeros(numColors, 3);
    
    for i = 1:numColors
        % Create slight variations around the base color for each flow vector
        colorShift = (i-1)/(numColors-1+eps) * 0.4; % Vary within 0.4 range
        colorVariations(i,:) = min(1, baseColor + [-0.15+colorShift, colorShift, 0.1-colorShift]);
    end
end

function boundary = calculateFeasibleRegionBoundary(allPathTimeCosts, allPathMoneyCosts)
    % Calculate the boundary of the feasible region based on path costs
    
    % Find unique money cost values
    uniqueMoneyValues = unique(allPathMoneyCosts);
    leftBoundaryX = [];
    leftBoundaryY = [];
    rightBoundaryX = [];
    rightBoundaryY = [];
    
    % For each money cost, find minimum and maximum time costs
    for i = 1:length(uniqueMoneyValues)
        currMoney = uniqueMoneyValues(i);
        sameMoneyIdx = abs(allPathMoneyCosts - currMoney) < 0.001;
        
        if sum(sameMoneyIdx) > 0
            timesForMoney = allPathTimeCosts(sameMoneyIdx);
            minTimeForMoney = min(timesForMoney);
            maxTimeForMoney = max(timesForMoney);
            
            % Add to boundary arrays
            leftBoundaryX = [leftBoundaryX; minTimeForMoney];
            leftBoundaryY = [leftBoundaryY; currMoney];
            rightBoundaryX = [rightBoundaryX; maxTimeForMoney];
            rightBoundaryY = [rightBoundaryY; currMoney];
        end
    end
    
    % Sort boundary points by money cost
    [leftBoundaryY, sortIdx] = sort(leftBoundaryY);
    leftBoundaryX = leftBoundaryX(sortIdx);
    
    [rightBoundaryY, sortIdx] = sort(rightBoundaryY);
    rightBoundaryX = rightBoundaryX(sortIdx);
    
    % Create boundary structure
    boundary = struct(...
        'leftX', leftBoundaryX, ...
        'leftY', leftBoundaryY, ...
        'rightX', rightBoundaryX, ...
        'rightY', rightBoundaryY ...
    );
end

function plotTimeMoneyCostRelationship(allPathCosts, boundary, colorVariations, legendLabels, params)
    % Plot time-money cost relationship figure
    
    % Create figure
    fig = figure('Name', 'Path Time-Money Cost Relationship', 'NumberTitle', 'off', ...
          'Position', params.FigurePosition);
    configureFigure(fig, params);
    
    % Combine boundaries to form closed region
    boundaryX = [boundary.leftX; flipud(boundary.rightX)];
    boundaryY = [boundary.leftY; flipud(boundary.rightY)];
    
    % Create hidden legend handles
    q = size(colorVariations, 1);
    h_legend = zeros(q+2, 1);
    
    % Create gradient fill
    h_legend(1) = patch('XData', boundaryX, 'YData', boundaryY, ...
          'FaceColor', [0.9 0.95 1], ... % Very light blue
          'EdgeColor', 'none', ... 
          'LineWidth', 0.1, ...
          'FaceAlpha', 1);
          
    % Draw boundary lines
    h_legend(2) = plot(boundary.leftX, boundary.leftY, '-', 'Color', [0.4 0.5 0.8], 'LineWidth', 0.1);
    plot(boundary.rightX, boundary.rightY, '-', 'Color', [0.4 0.5 0.8], 'LineWidth', 0.1);

    % Plot all points
    for i = 1:q
        costs = allPathCosts{i};
        if ~isempty(costs)
            % Plot with elegant lines
            plot(costs(:,1), costs(:,2), '-', 'Color', [colorVariations(i,:), 0.7], 'LineWidth', 1.2);
            
            % Plot elegant markers
            h_legend(i+2) = scatter(costs(:,1), costs(:,2), 30, colorVariations(i,:), 'o', ...
                'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.8);
        end
    end
    
    % Configure axes and labels
    configureAxes(params);
    
    % Add legend
    legend([h_legend(1), h_legend(3)], {'Feasible Region', 'Selected Paths'}, ...
        'Location', 'northeast', ...
        'FontName', params.FontName, 'FontSize', params.FontSize-1, ...
        'EdgeColor', [0.7, 0.7, 0.7], ...
        'Box', 'on');
    
    % Save figure
    saveFigure(fig, [params.SavePath 'path_time_money_relationship_']);
end

function plotPathCostsWithUpperLimit(allPathCosts, boundary, colorVariations, params)
    % Plot path costs with upper limit, differentiating feasible and infeasible vectors
    
    % Calculate cost upper limit - using midpoints of boundaries
    upperLimitX = calculateUpperLimit(boundary.leftX, boundary.rightX);
    upperLimitY = boundary.leftY;  % Money cost is the same
    
    % Calculate equilibrium line (T_eqm) - position varies based on money cost
    eqmLimitX = calculateEquilibriumLine(boundary.leftX, upperLimitX, upperLimitY);
    [T, M] = calculateTandM(eqmLimitX, upperLimitY);
    
    % Create figure
    fig = figure('Name', 'Path Costs with Upper Limit', 'NumberTitle', 'off', ...
          'Position', params.FigurePosition);
    configureFigure(fig, params);
    
    % Draw feasible region
    [h_legend, boundaryX, boundaryY] = drawFeasibleRegion(boundary);
    
    % Check each path flow vector for feasibility
    [feasiblePaths, infeasiblePaths, feasibleCosts, feasibleColors] = ...
        checkPathFeasibility(allPathCosts, upperLimitX, upperLimitY);
    
    % Draw paths and regions
    [h_legend, h_infeasible] = drawPathsAndRegions(h_legend, allPathCosts, feasiblePaths, infeasiblePaths, ...
        feasibleCosts, feasibleColors, colorVariations);
    
    % Draw equilibrium line (T_eqm) - use a distinctive purple color and thicker dash-dot pattern
    h_legend(3) = plot(eqmLimitX, upperLimitY, '-.', 'Color', [0.5 0.0 0.8], 'LineWidth', 2.5);
    
    % Draw upper limit line (T_max)
    h_legend(7) = plot(upperLimitX, upperLimitY, '-.', 'Color', [0.8 0.2 0.2], 'LineWidth', 2);
    % h_legend(7) = plot(T, M, '-.', 'Color', [0.8 0.2 0.2], 'LineWidth', 2);
    
    % Configure axes
    configureAxes(params);
    
    % Add legend
    addLegendForUpperLimitPlot(h_legend, feasiblePaths, infeasiblePaths);
    
    % Save figure
    saveFigure(fig, [params.SavePath 'path_costs_with_upper_limit_']);
end

function upperLimitX = calculateUpperLimit(leftBoundaryX, rightBoundaryX)
    % Calculate the upper limit X values (time costs) as midpoints between boundaries
    upperLimitX = (leftBoundaryX + rightBoundaryX) / 2;
end

function eqmLimitX = calculateEquilibriumLine(leftBoundaryX, upperLimitX, upperLimitY)
    % Calculate equilibrium line position between left boundary and upper limit
    % The width between eqmLimitX and upperLimitX increases as money cost decreases
    
    % Get normalized money costs
    moneyMax = max(upperLimitY);
    moneyMin = min(upperLimitY);
    moneyRange = moneyMax - moneyMin;
    
    if moneyRange > 0
        % Calculate weights: lower money costs get higher weights
        % These weights determine distance from upperLimitX (T_max)
        % When money cost is low, weight is high, distance from T_max is large
        weights = 0.7 + 0.1 * ((upperLimitY - moneyMin) / moneyRange);
    else
        weights = 0.5 * ones(size(upperLimitY));
    end
    
    % Calculate equilibrium line position
    % With low money cost (small weight), eqmLimitX will be closer to leftBoundaryX
    % With high money cost (large weight), eqmLimitX will be closer to upperLimitX
    % eqmLimitX = upperLimitX - 50./(10+upperLimitY);
    eqmLimitX = leftBoundaryX + (upperLimitX - leftBoundaryX) .* weights;
end


function [T, M] = calculateTandM(eqmLimitX, upperLimitY)
    % 计算分段线性函数 T = k*M + b + 50/(10+M)
    % eqmLimitX, upperLimitY: 分段端点
    T = [];
    M = [];
    for i = 1:length(eqmLimitX)-1
        % 当前区间的M
        M_segment = linspace(upperLimitY(i), upperLimitY(i+1), 50);
        % 计算斜率和截距
        k = (eqmLimitX(i+1) - eqmLimitX(i)) / (upperLimitY(i+1) - upperLimitY(i));
        b = eqmLimitX(i) - k * upperLimitY(i);
        % 计算T
        T_segment = k * M_segment + b + 50 ./ (10 + M_segment);
        % 拼接
        T = [T, T_segment];
        M = [M, M_segment];
    end
end

function [h_legend, boundaryX, boundaryY] = drawFeasibleRegion(boundary)
    % Draw the feasible region on the current axes
    
    % Combine boundaries to form closed region
    boundaryX = [boundary.leftX; flipud(boundary.rightX)];
    boundaryY = [boundary.leftY; flipud(boundary.rightY)];
    
    % Create legend handles
    h_legend = zeros(6, 1);
    
    % Draw feasible region
    h_legend(1) = patch('XData', boundaryX, 'YData', boundaryY, ...
          'FaceColor', [0.9 0.95 1], ... % Very light blue
          'EdgeColor', 'none', ... 
          'FaceAlpha', 1);
    
    % Draw boundary lines
    h_legend(2) = plot(boundary.leftX, boundary.leftY, '-', 'Color', [0.4 0.5 0.8], 'LineWidth', 0.1);
    plot(boundary.rightX, boundary.rightY, '-', 'Color', [0.4 0.5 0.8], 'LineWidth', 0.1);
end

function [feasiblePaths, infeasiblePaths, feasibleCosts, feasibleColors] = ...
    checkPathFeasibility(allPathCosts, eqmLimitX, upperLimitY)
    % Check which path flow vectors are feasible based on upper limit
    
    % Preallocate result arrays
    feasible_flags = cellfun(@(costs) ...
        ~isempty(costs) && ...
        all(costs(:,1) <= interp1(upperLimitY, eqmLimitX, costs(:,2), 'nearest', 'extrap')), ...
        allPathCosts);

    % Get feasible and infeasible paths using logical indexing
    feasiblePaths = find(feasible_flags);
    infeasiblePaths = find(~feasible_flags);
    
    % Extract feasible costs and colors
    feasibleCosts = allPathCosts(feasiblePaths);
    feasibleColors = cell(length(feasiblePaths), 1);
    
    % If there are feasible paths, prepare colors
    if ~isempty(feasiblePaths)
        for i = 1:length(feasiblePaths)
            feasibleColors{i} = [0.2, 0.6, 0.8]; % Default blue color
        end
    end
end

function [h_legend, h_infeasible] = drawPathsAndRegions(h_legend, allPathCosts, ...
    feasiblePaths, infeasiblePaths, feasibleCosts, feasibleColors, colorVariations)
    % Draw paths and regions for feasible and infeasible flow vectors
    
    h_infeasible = [];
    
    % Draw infeasible flow vectors (gray)
    infeasibleColor = [0.9, 0.9, 0.9];
    for i = 1:length(infeasiblePaths)
        idx = infeasiblePaths(i);
        costs = allPathCosts{idx};
        
        plot(costs(:,1), costs(:,2), '-', 'Color', [infeasibleColor, 0.7], 'LineWidth', 1.2);
        h_tmp = scatter(costs(:,1), costs(:,2), 30, infeasibleColor, 'o', 'filled', ...
            'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.3);
        
        if isempty(h_infeasible)
            h_infeasible = h_tmp;
            h_legend(5) = h_tmp;
        end
    end
    
    % Draw feasible flow vectors and region
    if ~isempty(feasibleCosts)
        % Draw feasible region if possible
        [feasibleRegion, h_feasible_region] = calculateAndDrawFeasibleRegion(feasibleCosts);
        
        % Update legend if feasible region was drawn
        if ~isempty(h_feasible_region)
            h_legend(6) = h_feasible_region;
        end
        
        % Draw each feasible vector
        for i = 1:length(feasibleCosts)
            costs = feasibleCosts{i};
            
            % Use the original color for the corresponding index if available
            if i <= length(feasiblePaths) && feasiblePaths(i) <= size(colorVariations, 1)
                color = colorVariations(feasiblePaths(i), :);
            else
                color = feasibleColors{i};
            end
            
            % Draw lines and points
            plot(costs(:,1), costs(:,2), '-', 'Color', [color, 0.7], 'LineWidth', 1.2);
            h_legend(4) = scatter(costs(:,1), costs(:,2), 30, color, 'o', 'filled', ...
                'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.8);
        end
    end
end

function [feasibleRegion, h_feasible_region] = calculateAndDrawFeasibleRegion(feasibleCosts)
    % Calculate and draw the region for feasible flow vectors
    
    feasibleRegion = [];
    h_feasible_region = [];
    
    try
        % Collect all feasible points
        allFeasibleX = [];
        allFeasibleY = [];
        
        for i = 1:length(feasibleCosts)
            costs = feasibleCosts{i};
            allFeasibleX = [allFeasibleX; costs(:,1)];
            allFeasibleY = [allFeasibleY; costs(:,2)];
        end
        
        if isempty(allFeasibleX)
            return;
        end
        
        % Collect all feasible points
        feasiblePoints = [allFeasibleX, allFeasibleY];
        
        % Group points by money cost
        uniqueMoneyValues = unique(round(feasiblePoints(:,2)*100)/100); % Round to two decimal places
        
        % Find min and max time costs for each money cost
        leftBoundX = [];
        leftBoundY = [];
        rightBoundX = [];
        rightBoundY = [];
        
        for i = 1:length(uniqueMoneyValues)
            currMoney = uniqueMoneyValues(i);
            sameMoneyIdx = abs(feasiblePoints(:,2) - currMoney) < 0.01;
            
            if sum(sameMoneyIdx) > 0
                timesForMoney = feasiblePoints(sameMoneyIdx, 1);
                minTimeForMoney = min(timesForMoney);
                maxTimeForMoney = max(timesForMoney);
                
                leftBoundX = [leftBoundX; minTimeForMoney];
                leftBoundY = [leftBoundY; currMoney];
                rightBoundX = [rightBoundX; maxTimeForMoney];
                rightBoundY = [rightBoundY; currMoney];
            end
        end
        
        % Sort boundary points
        [leftBoundY, sortIdx] = sort(leftBoundY);
        leftBoundX = leftBoundX(sortIdx);
        leftBoundX = ones(length(leftBoundY),1)*(min(leftBoundX)-1);
        
        [rightBoundY, sortIdx] = sort(rightBoundY);
        rightBoundX = rightBoundX(sortIdx);
        % Form closed region
        boundX = [leftBoundX; flipud(rightBoundX)];
        boundY = [leftBoundY; flipud(rightBoundY)];
        
        % Draw feasible region - use a more distinct light green color
        h_feasible_region = fill(boundX, boundY, [0.7, 0.95, 0.7], ...
            'FaceAlpha', 0.5, 'EdgeColor', [0.4, 0.6, 0.4], 'LineWidth', 0.1);
        
        feasibleRegion = [boundX, boundY];
    catch e
        fprintf('Unable to draw feasible flow region: %s\n', e.message);
    end
end

function configureFigure(fig, params)
    % Configure figure appearance
    set(fig, 'Color', 'white');
    set(gca, 'FontName', params.FontName, 'FontSize', params.FontSize, ...
        'Box', 'on', 'LineWidth', 1);
    hold on;
end

function configureAxes(params)
    % Configure axes appearance
    if params.ShowGrid
        grid on;
        grid minor;
        set(gca, 'GridAlpha', 0.1, 'MinorGridAlpha', 0.05, 'Layer', 'top');
    end
    
    % Add labels
    xlabel('Time Cost', 'FontSize', params.FontSize+2, 'FontWeight', 'bold');
    ylabel('Money Cost', 'FontSize', params.FontSize+2, 'FontWeight', 'bold');
end

function addLegendForUpperLimitPlot(h_legend, feasiblePaths, infeasiblePaths)
    % Add legend to upper limit plot with simplified labels
    
    % Create arrays for legend handles and labels
    legendHandles = [];
    legendItems = {};
    
    % 1. Feasible Region (always show first)
    if h_legend(1) ~= 0
        legendHandles(end+1) = h_legend(1);
        legendItems{end+1} = 'Feasible Region';
    end
    
    % 2. Feasible Flow Region (show second if available)
    if ~isempty(feasiblePaths) && length(h_legend) >= 6 && h_legend(6) ~= 0
        legendHandles(end+1) = h_legend(6);
        legendItems{end+1} = 'Feasible Flow Region';
    end
    
    % 3. T_eqm
    if h_legend(3) ~= 0
        legendHandles(end+1) = h_legend(3);
        legendItems{end+1} = 'T_{eqm}';
    end
    
    % 4. T_max
    if length(h_legend) >= 7 && h_legend(7) ~= 0
        legendHandles(end+1) = h_legend(7);
        legendItems{end+1} = 'T_{max}';
    end
    
    legend(legendHandles, legendItems, ...
        'Location', 'best', ...
        'FontName', 'Arial', 'FontSize', 9, ...
        'EdgeColor', [0.7, 0.7, 0.7], ...
        'Box', 'on');
end

% function printFeasibilityResults(feasiblePaths, infeasiblePaths)
%     % Print feasibility results
%     fprintf('可行流量方案索引: %s\n', mat2str(feasiblePaths));
%     fprintf('不可行流量方案索引: %s\n', mat2str(infeasiblePaths));
% end

function saveFigure(fig, baseFilename)
    % Save figure as PDF with formatted filename including zeta and subset_index
    
    % 获取当前的zeta和subset_index值
    zeta = getappdata(0, 'current_zeta');
    subset_index = getappdata(0, 'current_subset_index');
    
    % 如果找不到参数，使用默认保存方式
    if isempty(zeta) || isempty(subset_index)
        % 保存为PNG格式，带时间戳
        figFile = sprintf('%s%s.png', baseFilename, datestr(now, 'yyyymmdd_HHMMSS'));
        print(fig, figFile, '-dpng', '-r300');
    else
        % 保存为PDF格式，使用参数化文件名，不包含时间戳
        if contains(baseFilename, 'path_time_money_relationship')
            baseOutputFile = sprintf('%spath_time_money_zeta%d_subset%d', ...
                     regexprep(baseFilename, 'path_time_money_relationship_.*', ''), ...
                     zeta, subset_index);
        elseif contains(baseFilename, 'path_costs_with_upper_limit')
            baseOutputFile = sprintf('%spath_costs_upper_limit_zeta%d_subset%d', ...
                     regexprep(baseFilename, 'path_costs_with_upper_limit_.*', ''), ...
                     zeta, subset_index);
        else
            baseOutputFile = sprintf('%s_zeta%d_subset%d', baseFilename, zeta, subset_index);
        end
        
        % 生成PDF和FIG文件名
        pdfFile = [baseOutputFile '.pdf'];
        figFile = [baseOutputFile '.fig'];
        
        % 确保目录存在
        [directory,~,~] = fileparts(pdfFile);
        if ~isempty(directory) && ~exist(directory, 'dir')
            mkdir(directory);
        end
        
        % 设置PDF大小与图窗大小一致
        set(fig, 'Units', 'inches');
        figPos = get(fig, 'Position');
        
        % 设置纸张大小与图窗大小一致
        set(fig, 'PaperPositionMode', 'auto');
        set(fig, 'PaperUnits', 'inches');
        set(fig, 'PaperSize', [figPos(3), figPos(4)]);
        
        % 保存为PDF格式，使用矢量渲染引擎
        print(fig, pdfFile, '-dpdf', '-painters');
        fprintf('PDF图像已保存为: %s (大小: %.1f × %.1f 英寸)\n', pdfFile, figPos(3), figPos(4));
        
        % 保存为FIG格式，便于将来编辑
        savefig(fig, figFile);
        fprintf('MATLAB图形已保存为: %s\n', figFile);
    end
    
    hold off;
end

function RT = calculateRealTime(x, M, freeFlowTime, maxCapacity)
    % Calculate real travel time
    index = find(sum(M)~=0);
    time = calculateTravelTime(x(:,index), freeFlowTime(index), maxCapacity(index));
    pathTime = time * M(:,index)';
    RT = pathTime + 15 * (1 - exp(-0.02 * pathTime));
end

function t = calculateTravelTime(x, t0, C)
    % Calculate travel time
    t = t0 .* (1 + 0.15 .* (x ./ C).^4);
end 