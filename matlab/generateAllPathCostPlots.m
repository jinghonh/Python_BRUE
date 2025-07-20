function generateAllPathCostPlots(rangeMin, rangeMax, varargin)
    % 一键生成所有zeta和subset_index组合的路径成本图
    %
    % 输入参数:
    %   rangeMin     - 每个维度的最小值数组 [d1_min, d2_min, ...]
    %   rangeMax     - 每个维度的最大值数组 [d1_max, d2_max, ...]
    %   varargin     - 可选参数:
    %                   'SavePath': 图像保存路径 (默认: 'results/')
    %                   'FigurePosition': 图像位置和大小 (默认: [100, 100, 800, 600])
    %                   'FontName': 字体名称 (默认: 'Arial')
    %                   'FontSize': 字体大小 (默认: 10)
    %                   'ShowGrid': 是否显示网格 (默认: true)

    % 解析可选参数
    p = inputParser;
    defaultSavePath = 'results/pdf_outputs/';
    defaultFigPosition = [100, 100, 800, 600];
    defaultFontName = 'Arial';
    defaultFontSize = 10;
    defaultShowGrid = true;
    
    addParameter(p, 'SavePath', defaultSavePath);
    addParameter(p, 'FigurePosition', defaultFigPosition);
    addParameter(p, 'FontName', defaultFontName);
    addParameter(p, 'FontSize', defaultFontSize);
    addParameter(p, 'ShowGrid', defaultShowGrid);
    
    parse(p, varargin{:});
    params = p.Results;

    % 确保结果目录存在
    if ~exist(params.SavePath, 'dir')
        mkdir(params.SavePath);
    end
    
    % 定义zeta和subset_index的所有可能组合
    combinations = [
        15, 0;
        15, 1;
        31, 0;
        31, 1;
        31, 2
    ];
    
    % 遍历所有组合并生成图像
    fprintf('开始生成所有组合的路径成本图...\n');
    
    for i = 1:size(combinations, 1)
        zeta = combinations(i, 1);
        subset_index = combinations(i, 2);
        
        fprintf('处理组合: zeta=%d, subset_index=%d\n', zeta, subset_index);
        
        % 调用analyzeTrafficNetwork生成相应的图
        analyzeTrafficNetwork(zeta, rangeMin, rangeMax, subset_index, ...
            'SavePath', params.SavePath, ...
            'FigurePosition', params.FigurePosition, ...
            'FontName', params.FontName, ...
            'FontSize', params.FontSize, ...
            'ShowGrid', params.ShowGrid, ...
            'ZetaValue', zeta, ...
            'SubsetIndex', subset_index);
            
        fprintf('完成组合: zeta=%d, subset_index=%d\n', zeta, subset_index);
    end
    
    % 显示生成的PDF文件列表
    pdfFiles = dir(fullfile(params.SavePath, '*.pdf'));
    if ~isempty(pdfFiles)
        fprintf('\n已生成的PDF文件:\n');
        for j = 1:length(pdfFiles)
            fprintf('  - %s\n', pdfFiles(j).name);
        end
    else
        fprintf('\n警告：未找到生成的PDF文件。请检查保存路径：%s\n', params.SavePath);
    end
    
    fprintf('所有图像生成完毕!\n');
end 