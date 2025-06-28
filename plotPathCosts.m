function plotPathCosts(totalValidFlow, relationMatrix, selectedIndices)
    % Plot time and money costs for each path with scientific style
    % Input parameters:
    %   totalValidFlow  - All feasible flow matrix (M x n)
    %   relationMatrix  - Relation matrix
    %   selectedIndices - Selected flow vector indices
    
    % Define constants
    money = [20, 15, 1, 0, 0, 0, 0, 1];
    freeFlowTime = [18,22.5,12,24,2.4,6,24,12];
    maxCapacity = [3600,3600,1800,1800,1800,1800,1800,1800];
    money = money * relationMatrix';  % 1 x n
    
    % Process data for plotting
    q = length(selectedIndices);
    colors = parula(q);  % Scientific colormap
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
    
    % ----------------- Figure 1: Time-Money Cost Relationship -----------------
    % Create separate figure for time-money cost relationship
    fig1 = figure('Name', 'Path Time-Money Cost Relationship', 'NumberTitle', 'off', 'Position', [100, 100, 800, 600]);
    set(fig1, 'Color', 'white');
    set(gca, 'FontName', 'Arial', 'FontSize', 10, 'Box', 'on', 'LineWidth', 1);
    
    % Create a more elegant color scheme
    % Use a subdued blue gradient for better visualization
    baseColor = [0.2 0.4 0.8]; % Base blue color
    colorVariations = [];
    
    for i = 1:q
        % Create slight variations around the base color for each flow vector
        colorShift = (i-1)/(q-1+eps) * 0.4; % Vary within 0.4 range
        colorVariations(i,:) = min(1, baseColor + [-0.15+colorShift, colorShift, 0.1-colorShift]);
    end
    
    hold on;
    
    % 改进的边界计算逻辑 - 使用金钱成本固定的特性
    % 用所有散点数据进行计算
    allX = allPathTimeCosts;
    allY = allPathMoneyCosts;
    
    % 按金钱成本分组，为每个金钱成本找到时间的最小值和最大值
    % 这种方法比凸包更适合这种特定的数据结构
    
    % 找到所有唯一的金钱成本值
    uniqueMoneyValues = unique(allY);
    leftBoundaryX = [];
    leftBoundaryY = [];
    rightBoundaryX = [];
    rightBoundaryY = [];
    
    % 为每个金钱成本找到对应的最小和最大时间成本
    for i = 1:length(uniqueMoneyValues)
        currMoney = uniqueMoneyValues(i);
        % 找到具有相同金钱成本的所有点
        sameMoneyIdx = abs(allY - currMoney) < 0.001;
        
        if sum(sameMoneyIdx) > 0
            timesForMoney = allX(sameMoneyIdx);
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
    
    % 合并上下边界以形成封闭区域
    boundaryX = [leftBoundaryX; flipud(rightBoundaryX)];
    boundaryY = [leftBoundaryY; flipud(rightBoundaryY)];
    
    % 先创建隐藏的图例句柄
    h_legend = zeros(q+2, 1);
    
    % 创建渐变填充
    h_legend(1) = patch('XData', boundaryX, 'YData', boundaryY, ...
          'FaceColor', [0.9 0.95 1], ... % 非常浅的蓝色
          'EdgeColor', 'none', ... 
          'LineWidth', 1.5, ...
          'FaceAlpha', 1);
          
    % 绘制边界线
    h_legend(2) = plot(leftBoundaryX, leftBoundaryY, '-', 'Color', [0.4 0.5 0.8], 'LineWidth', 1.5);
    plot(rightBoundaryX, rightBoundaryY, '-', 'Color', [0.4 0.5 0.8], 'LineWidth', 1.5);

    % 使用所有点计算边界后绘制散点
    for i = 1:q
        costs = allPathCosts{i};
        if ~isempty(costs)
            % Plot with thinner, more elegant lines
            plot(costs(:,1), costs(:,2), '-', 'Color', [colorVariations(i,:), 0.7], 'LineWidth', 1.2);
            
            % Use smaller, more elegant markers
            h_legend(i+2) = scatter(costs(:,1), costs(:,2), 30, colorVariations(i,:), 'o', 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.8);
        end
    end
    
    % 保存边界以便在第二个图中使用
    saveBoundaryX = boundaryX;
    saveBoundaryY = boundaryY;
    
    % Add subtle grid lines for scientific style
    grid on;
    grid minor;
    set(gca, 'GridAlpha', 0.1, 'MinorGridAlpha', 0.05, 'Layer', 'top');
    
    % Figure properties
    xlabel('Time Cost', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Money Cost', 'FontSize', 12, 'FontWeight', 'bold');
    % title('Time-Money Cost Relationship for Each Path', 'FontSize', 14, 'FontWeight', 'bold');
    
    % 添加科研风格的图例
    % 只显示概括性图例
    legend([h_legend(1), h_legend(2), h_legend(3)], {'Feasible Region', 'Boundary', 'Selected Paths'}, ...
        'Location', 'best', ...
        'FontName', 'Arial', 'FontSize', 9, ...
        'EdgeColor', [0.7, 0.7, 0.7], ...
        'Box', 'on');

    
    figFile = sprintf('results/path_time_money_relationship_%s.png', datestr(now, 'yyyymmdd_HHMMSS'));
    print(figFile, '-dpng', '-r300');
    hold off;
    
    % 调用新函数以生成带有上限线的新图
    plotPathCostsWithUpperLimit(allPathCosts, leftBoundaryX, leftBoundaryY, rightBoundaryX, rightBoundaryY, colorVariations);
    
    % ----------------- Figure 2: Path Cost Scatter Plot -----------------
    % Create separate figure for path cost scatter plot
    % fig2 = figure('Name', 'Path Cost Scatter Plot', 'NumberTitle', 'off', 'Position', [100, 100, 600, 500]);
    % set(fig2, 'Color', 'white');
    % set(gca, 'FontName', 'Arial', 'FontSize', 10, 'Box', 'on', 'LineWidth', 1);
    % 
    % hold on;
    % for i = 1:q
    %     indices = allFlowVectorIndices == i;
    %     scatter(allPathTimeCosts(indices), allPathMoneyCosts(indices), 70, colors(i,:), 'filled', 'MarkerEdgeColor', 'black', 'MarkerEdgeAlpha', 0.3);
    % end
    % 
    % % 为第二个图添加边界线
    % if exist('saveBoundaryX', 'var') && exist('saveBoundaryY', 'var') && ~isempty(saveBoundaryX) && ~isempty(saveBoundaryY)
    %     patch('XData', saveBoundaryX, 'YData', saveBoundaryY, ...
    %           'FaceColor', [0.95 0.95 1], ... % 非常浅的蓝色
    %           'EdgeColor', [0.5 0.5 0.8], ... % 中等蓝色边缘
    %           'LineWidth', 1.2, ...
    %           'FaceAlpha', 0.3);  % 更低的透明度
    % 
    %     % 绘制边界线
    %     plot(saveBoundaryX, saveBoundaryY, '-', 'Color', [0.5 0.5 0.8], 'LineWidth', 1.2);
    % end
    % 
    % % Figure properties
    % xlabel('Time Cost', 'FontSize', 12, 'FontWeight', 'bold');
    % ylabel('Money Cost', 'FontSize', 12, 'FontWeight', 'bold');
    % title('Path Cost Scatter Plot', 'FontSize', 14, 'FontWeight', 'bold');
    % grid on;
    % grid minor;
    % set(gca, 'GridAlpha', 0.15, 'MinorGridAlpha', 0.1, 'Layer', 'top');
    % 
    % % No legend as requested
    % figFile = sprintf('results/path_cost_scatter_%s.png', datestr(now, 'yyyymmdd_HHMMSS'));
    % print(figFile, '-dpng', '-r300');
    % hold off;
    % 
    % % ----------------- Figure 3: Time Cost Distribution -----------------
    % % Create separate figure for time cost distribution
    % fig3 = figure('Name', 'Path Time Cost Distribution', 'NumberTitle', 'off', 'Position', [100, 100, 600, 500]);
    % set(fig3, 'Color', 'white');
    % set(gca, 'FontName', 'Arial', 'FontSize', 10, 'Box', 'on', 'LineWidth', 1);
    % 
    % % Optimize x-axis by using custom box plot with more space
    % positions = 1:q;
    % boxplot(allPathTimeCosts, allFlowVectorIndices, 'Positions', positions, 'Width', 0.6, 'Symbol', 'k+');
    % 
    % % Customize box plot colors
    % h = findobj(gca, 'Tag', 'Box');
    % for j = 1:length(h)
    %     patch(get(h(j), 'XData'), get(h(j), 'YData'), colors(q-j+1,:), 'FaceAlpha', 0.7);
    % end
    % 
    % % Add custom scatter points for better visualization
    % hold on;
    % for i = 1:q
    %     indices = allFlowVectorIndices == i;
    %     % Add jitter to x positions for better visibility
    %     x_jitter = positions(i) + (rand(sum(indices),1)-0.5)*0.3;
    %     scatter(x_jitter, allPathTimeCosts(indices), 20, colors(i,:), 'filled', 'MarkerEdgeColor', 'black', 'MarkerEdgeAlpha', 0.3);
    % end
    % 
    % % Figure properties
    % title('Path Time Cost Distribution by Flow Vector', 'FontSize', 14, 'FontWeight', 'bold');
    % xlabel('Flow Vector', 'FontSize', 12, 'FontWeight', 'bold');
    % ylabel('Time Cost', 'FontSize', 12, 'FontWeight', 'bold');
    % 
    % % Fix x-axis overcrowding
    % xticks(positions);
    % xticklabels(cellfun(@(x) sprintf('FV%d', str2double(regexp(x, '\d+', 'match'))), legendLabels, 'UniformOutput', false));
    % xtickangle(0);
    % 
    % grid on;
    % set(gca, 'GridAlpha', 0.15, 'Layer', 'top');
    % 
    % % No legend as requested
    % 
    % figFile = sprintf('results/path_time_distribution_%s.png', datestr(now, 'yyyymmdd_HHMMSS'));
    % print(figFile, '-dpng', '-r300');
    hold off;
end

function plotPathCostsWithUpperLimit(allPathCosts, leftBoundaryX, leftBoundaryY, rightBoundaryX, rightBoundaryY, colorVariations)
    % 创建一个新图，在可行区域内添加成本上限线，并区分可行与不可行流量分配方案
    % 输入参数:
    %   allPathCosts      - 路径成本数据单元格数组
    %   leftBoundaryX     - 左边界X坐标
    %   leftBoundaryY     - 左边界Y坐标
    %   rightBoundaryX    - 右边界X坐标
    %   rightBoundaryY    - 右边界Y坐标
    %   colorVariations   - 颜色变化数组
    
    % 计算可行成本上限 - 使用边界的中点
    upperLimitX = [];
    upperLimitY = [];
    
    for i = 1:length(leftBoundaryY)
        % 对于每一个金钱成本，计算时间成本上限
        midTimeForMoney = (leftBoundaryX(i) + rightBoundaryX(i)) / 2;
        upperLimitX = [upperLimitX; midTimeForMoney];
        upperLimitY = [upperLimitY; leftBoundaryY(i)];
    end
    
    % 创建新图
    fig_new = figure('Name', 'Path Costs with Upper Limit', 'NumberTitle', 'off', 'Position', [100, 100, 800, 600]);
    set(fig_new, 'Color', 'white');
    set(gca, 'FontName', 'Arial', 'FontSize', 10, 'Box', 'on', 'LineWidth', 1);
    
    hold on;
    
    % 合并边界以形成封闭区域
    boundaryX = [leftBoundaryX; flipud(rightBoundaryX)];
    boundaryY = [leftBoundaryY; flipud(rightBoundaryY)];
    
    % 创建图例句柄
    h_legend = zeros(5, 1);
    
    % 绘制可行区域
    h_legend(1) = patch('XData', boundaryX, 'YData', boundaryY, ...
          'FaceColor', [0.9 0.95 1], ... % 非常浅的蓝色
          'EdgeColor', 'none', ... 
          'FaceAlpha', 1);
    
    % 绘制边界线
    h_legend(2) = plot(leftBoundaryX, leftBoundaryY, '-', 'Color', [0.4 0.5 0.8], 'LineWidth', 1.5);
    plot(rightBoundaryX, rightBoundaryY, '-', 'Color', [0.4 0.5 0.8], 'LineWidth', 1.5);
    
    % 检查每个路径流量方案是否可行
    feasiblePaths = [];
    infeasiblePaths = [];
    
    % 用于存储不可行方案的图例句柄
    h_infeasible = [];
    
    for i = 1:length(allPathCosts)
        costs = allPathCosts{i};
        if ~isempty(costs)
            % 检查每个点是否在上限线下
            isPointFeasible = zeros(size(costs, 1), 1);
            
            for j = 1:size(costs, 1)
                currPoint = costs(j, :);  % [时间成本, 金钱成本]
                
                % 找到最接近的金钱成本点
                [~, idx] = min(abs(upperLimitY - currPoint(2)));
                
                % 检查时间成本是否低于或等于上限
                isPointFeasible(j) = currPoint(1) <= upperLimitX(idx);
            end
            
            % 如果所有点都可行，则整个方案可行
            if ~all(isPointFeasible)
                % 不可行方案（）
                infeasiblePaths = [infeasiblePaths, i];
                infeasibleColor = [0.9, 0.9, 0.9]; 
                
                plot(costs(:,1), costs(:,2), '-', 'Color', [infeasibleColor, 0.7], 'LineWidth', 1.2);
                h_tmp = scatter(costs(:,1), costs(:,2), 30, infeasibleColor, 'o', 'filled', ...
                    'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.3);
                
                if isempty(h_infeasible)
                    h_infeasible = h_tmp;
                    h_legend(5) = h_tmp;
                end
            else
                feasiblePaths = [feasiblePaths, i];
                
                % 绘制可行方案（蓝色）
                plot(costs(:,1), costs(:,2), '-', 'Color', [colorVariations(i,:), 0.7], 'LineWidth', 1.2);
                h_legend(4) = scatter(costs(:,1), costs(:,2), 30, colorVariations(i,:), 'o', 'filled', ...
                    'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.8);
            end
        end
    end
    % 绘制可行成本上限折线
    h_legend(3) = plot(upperLimitX, upperLimitY, '--', 'Color', [0.8 0.2 0.2], 'LineWidth', 2);
    
    % 添加网格线
    grid on;
    grid minor;
    set(gca, 'GridAlpha', 0.1, 'MinorGridAlpha', 0.05, 'Layer', 'top');
    
    % 图形属性
    xlabel('Time Cost', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Money Cost', 'FontSize', 12, 'FontWeight', 'bold');
    
    % 添加图例
    legendItems = {'Feasible Region', 'Boundary', 'Cost Upper Limit'};
    legendHandles = [h_legend(1), h_legend(2), h_legend(3)];
    
    if ~isempty(feasiblePaths)
        legendItems{end+1} = 'Feasible Flow Vectors';
        legendHandles(end+1) = h_legend(4);
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
    
    % 添加可行和不可行方案的信息
    fprintf('可行流量方案索引: %s\n', mat2str(feasiblePaths));
    fprintf('不可行流量方案索引: %s\n', mat2str(infeasiblePaths));
    
    % 保存图形
    figFile = sprintf('results/path_costs_with_upper_limit_%s.png', datestr(now, 'yyyymmdd_HHMMSS'));
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