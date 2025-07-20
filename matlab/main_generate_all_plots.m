% =========================================================================
% main_generate_all_plots.m
% 一键生成所有zeta和subset_index组合的路径成本PDF图
% =========================================================================

% 清除工作区和关闭所有图窗口
clear;
close all;

% 设置初始化参数
rangeMin = [0, 0, 0, 0, 0]; % 每个维度的最小值
rangeMax = [6000, 6000, 6000, 6000, 6000]; % 每个维度的最大值

% 设置输出目录
outputDir = 'results/pdf_outputs/';

% 确保输出目录存在
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
    fprintf('已创建PDF输出目录: %s\n', outputDir);
end

% 设置图像参数
figPosition = [100, 100, 800, 600];
fontName = 'Arial';
fontSize = 10;
showGrid = true;

% 调用生成函数
fprintf('开始批量生成所有zeta和subset_index组合的PDF路径成本图...\n');
tic;
generateAllPathCostPlots(rangeMin, rangeMax, ...
    'SavePath', outputDir, ...
    'FigurePosition', figPosition, ...
    'FontName', fontName, ...
    'FontSize', fontSize, ...
    'ShowGrid', showGrid);
elapsedTime = toc;

fprintf('\n批量生成完成!\n');
fprintf('总耗时: %.2f 秒\n', elapsedTime);
fprintf('输出目录: %s\n', outputDir);

% 列出生成的文件
pdfFiles = dir(fullfile(outputDir, '*.pdf'));
if isempty(pdfFiles)
    fprintf('\n警告: 未找到PDF文件! 请检查错误消息。\n');
else
    fprintf('\n生成了 %d 个PDF文件:\n', length(pdfFiles));
    for i = 1:length(pdfFiles)
        fprintf('  - %s (%s)\n', pdfFiles(i).name, formatFileSize(pdfFiles(i).bytes));
    end
end

% 辅助函数：格式化文件大小
function sizeStr = formatFileSize(bytes)
    if bytes < 1024
        sizeStr = sprintf('%d B', bytes);
    elseif bytes < 1024^2
        sizeStr = sprintf('%.1f KB', bytes/1024);
    else
        sizeStr = sprintf('%.1f MB', bytes/(1024^2));
    end
end 