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
    parse(p, varargin);
    
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
    
    % Draw upper limit line
    h_legend(3) = plot(upperLimitX, upperLimitY, '--', 'Color', [0.8 0.2 0.2], 'LineWidth', 2);
    
    % Configure axes
    configureAxes(params);
    
    % Add legend
    addLegendForUpperLimitPlot(h_legend, feasiblePaths, infeasiblePaths);
    
    % Print feasibility results
    printFeasibilityResults(feasiblePaths, infeasiblePaths);
    
    % Save figure
    saveFigure(fig, [params.SavePath 'path_costs_with_upper_limit_']);
end

function upperLimitX = calculateUpperLimit(leftBoundaryX, rightBoundaryX)
    % Calculate the upper limit X values (time costs) as midpoints between boundaries
    upperLimitX = (leftBoundaryX + rightBoundaryX) / 2;
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
    checkPathFeasibility(allPathCosts, upperLimitX, upperLimitY)
    % Check which path flow vectors are feasible based on upper limit
    
    % Preallocate result arrays
    feasible_flags = cellfun(@(costs) ...
        ~isempty(costs) && ...
        all(costs(:,1) <= interp1(upperLimitY, upperLimitX, costs(:,2), 'nearest', 'extrap')), ...
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
        
        [rightBoundY, sortIdx] = sort(rightBoundY);
        rightBoundX = rightBoundX(sortIdx);
        
        % Form closed region
        boundX = [leftBoundX; flipud(rightBoundX)];
        boundY = [leftBoundY; flipud(rightBoundY)];
        
        % Draw feasible region
        h_feasible_region = fill(boundX, boundY, [0.8, 1, 0.8], ...
            'FaceAlpha', 0.5, 'EdgeColor', [0.4, 0.5, 0.8], 'LineWidth', 0.1);
        
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
    % Add legend to upper limit plot
    legendItems = {'Feasible Region', 'Boundary', 'Cost Upper Limit'};
    legendHandles = [h_legend(1), h_legend(2), h_legend(3)];
    
    if ~isempty(feasiblePaths)
        legendItems{end+1} = 'Feasible Flow Vectors';
        legendHandles(end+1) = h_legend(4);
        
        if length(h_legend) >= 6 && h_legend(6) ~= 0
            legendItems{end+1} = 'Feasible Flow Region';
            legendHandles(end+1) = h_legend(6);
        end
    end
    
    if ~isempty(infeasiblePaths)
        legendItems{end+1} = 'Infeasible Flow Vectors';
        legendHandles(end+1) = h_legend(5);
    end
    
    legend(legendHandles, legendItems, ...
        'Location', 'best', ...
        'FontName', 'Arial', 'FontSize', 9, ...
        'EdgeColor', [0.7, 0.7, 0.7], ...
        'Box', 'on');
end

function printFeasibilityResults(feasiblePaths, infeasiblePaths)
    % Print feasibility results
    fprintf('可行流量方案索引: %s\n', mat2str(feasiblePaths));
    fprintf('不可行流量方案索引: %s\n', mat2str(infeasiblePaths));
end

function saveFigure(fig, baseFilename)
    % Save figure with timestamp
    figFile = sprintf('%s%s.png', baseFilename, datestr(now, 'yyyymmdd_HHMMSS'));
    print(figFile, '-dpng', '-r300');
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