function plotPathCosts(totalValidFlow, relationMatrix, selectedIndices)
    % Plot time and money costs for each path
    % Input parameters:
    %   totalValidFlow  - All feasible flow matrix (M x n)
    %   relationMatrix  - Relation matrix
    %   selectedIndices - Selected flow vector indices
    
    % Define constants
    money = [20, 15, 1, 0, 0, 0, 0, 1];
    freeFlowTime = [18,22.5,12,24,2.4,6,24,12];
    maxCapacity = [3600,3600,1800,1800,1800,1800,1800,1800];
    money = money * relationMatrix';  % 1 x n
    
    % Create figure
    figure('Name', 'Path Cost Analysis', 'NumberTitle', 'off');
    hold on;
    
    % Calculate and plot costs for each selected flow vector
    q = length(selectedIndices);
    colors = jet(q);  % Create different colors for each flow vector
    legendLabels = cell(q, 1);
    
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
        pathMoneyCosts(pathTimeCosts == 0) = [];
        pathTimeCosts(pathTimeCosts == 0) = [];
        
        % Combine and sort costs by money cost
        costs = [pathTimeCosts', pathMoneyCosts'];
        [~, sortIdx] = sort(pathMoneyCosts);
        costs = costs(sortIdx, :);
        
        % Plot sorted points and lines
        h = plot(costs(:,1), costs(:,2), 'o-', 'Color', colors(i,:), ...
            'LineWidth', 1.5, 'MarkerSize', 6);
        
        % Create legend labels
        legendLabels{i} = sprintf('Flow Vector %d', i);
    end
    
    % Set figure properties
    xlabel('Time Cost');
    ylabel('Money Cost');
    title('Time-Money Cost Relationship for Each Path (Sorted by Money Cost)');
    grid on;
    legend(legendLabels, 'Location', 'best');
    hold off;
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