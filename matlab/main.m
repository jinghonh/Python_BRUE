% 清理工作空间
clear;
clc;
close all;

% 解析命令行参数
args = parseInputArgs();

% 设置参数
zeta = args.zeta;          % 可选值：8, 15 或 31
subset_index = args.subset;   % zeta=15时可选0,1；zeta=31时可选0,1,2
forceRecalculate = args.force; % 是否强制重新计算，不使用缓存
useConfigFile = args.useConfig; % 是否使用配置文件管理搜索范围

% 验证输入参数
if ~ismember(zeta, [8, 15, 31])
    error('zeta必须为8、15或31');
end

if zeta == 8 && ~ismember(subset_index, 0)
    error('当zeta=8时，subset_index必须为0');
end

if zeta == 15 && ~ismember(subset_index, [0, 1])
    error('当zeta=15时，subset_index必须为0或1');
end

if zeta == 31 && ~ismember(subset_index, [0, 1, 2])
    error('当zeta=31时，subset_index必须为0、1或2');
end

% 加载搜索范围（从配置文件或使用默认值）
if useConfigFile
    % 从配置文件加载搜索范围
    config = loadSaveConfig('load', zeta, subset_index);
    rangeMin = config.rangeMin;
    rangeMax = config.rangeMax;
else
    % 使用硬编码的范围
    if zeta == 8
        switch subset_index
            case 0  % zeta 8 0
                rangeMin = [0, 0, 0];
                rangeMax = [10000, 10000, 10000];
        end
    elseif zeta == 15
        switch subset_index
            case 0  % zeta 15 0
                rangeMin = [3201, 3151, 1201];
                rangeMax = [5002, 4800, 2400];
            case 1  % zeta 15 1
                rangeMin = [4200, 3000, 0, 0, 0];
                rangeMax = [4900, 4200, 1610, 1610, 2000];
        end
    elseif zeta == 31
        switch subset_index
            case 0  % zeta 31 0
                rangeMin = [1588, 2800, 1198];
                rangeMax = [4903, 5823, 2840];
            case 1  % zeta 31 1
                rangeMin = [500, 1200, 0, 0, 0];
                rangeMax = [4899, 6001, 2507, 2503, 2717];
            case 2  % zeta 31 2
                rangeMin = [4000, 3000, 0, 0, 902, 0];
                rangeMax = [5007, 4201, 3909, 3909, 1353, 1601];
        end
    end
end

% 显示当前参数设置
fprintf('当前参数设置：\n');
fprintf('zeta = %d\n', zeta);
fprintf('subset_index = %d\n', subset_index);
fprintf('rangeMin = [%s]\n', num2str(rangeMin));
fprintf('rangeMax = [%s]\n', num2str(rangeMax));
fprintf('强制重新计算 = %d\n', forceRecalculate);
fprintf('使用配置文件 = %d\n', useConfigFile);

% 创建结果目录（如果不存在）
if ~exist('results', 'dir')
    mkdir('results');
    fprintf('创建结果目录: results\n');
end

% 运行分析
tic
analyzeTrafficNetwork(zeta, rangeMin, rangeMax, subset_index, 'ForceRecalculate', forceRecalculate, 'UseConfigFile', useConfigFile); 
toc

function args = parseInputArgs()
    % 使用inputParser解析MATLAB工作区中的变量作为参数
    % 
    % 返回:
    %   args - 包含所有解析参数的结构体
    %     .zeta - zeta值（默认=8）
    %     .subset - subset索引（默认=0）
    %     .force - 是否强制重新计算（默认=false）
    %     .useConfig - 是否使用配置文件（默认=true）
    
    % 初始化默认参数
    defaultZeta = 8;
    defaultSubset = 0;
    defaultForce = false;
    defaultUseConfig = true;
    
    % 从工作区获取参数（如果存在）
    try
        if evalin('base', 'exist(''zeta'', ''var'')')
            zeta = evalin('base', 'zeta');
        else
            zeta = defaultZeta;
        end
        
        if evalin('base', 'exist(''subset_index'', ''var'')')
            subset = evalin('base', 'subset_index');
        else
            subset = defaultSubset;
        end
        
        if evalin('base', 'exist(''forceRecalculate'', ''var'')')
            force = evalin('base', 'forceRecalculate');
        else
            force = defaultForce;
        end
        
        if evalin('base', 'exist(''useConfigFile'', ''var'')')
            useConfig = evalin('base', 'useConfigFile');
        else
            useConfig = defaultUseConfig;
        end
    catch
        % 如果出错，使用默认值
        zeta = defaultZeta;
        subset = defaultSubset;
        force = defaultForce;
        useConfig = defaultUseConfig;
    end
    
    % 允许从命令行输入覆盖
    % 尝试解析命令行参数（如果有的话）
    args = struct('zeta', zeta, 'subset', subset, 'force', force, 'useConfig', useConfig);
    
    % 检查是否有命令行参数（通过变量ans传递）
    try
        if evalin('base', 'exist(''ans'', ''var'')')
            cmdArgs = evalin('base', 'ans');
            
            % 如果ans是字符串数组，说明有命令行参数
            if iscell(cmdArgs) && ~isempty(cmdArgs)
                i = 1;
                while i <= length(cmdArgs)
                    if ischar(cmdArgs{i}) || isstring(cmdArgs{i})
                        switch lower(cmdArgs{i})
                            case '-zeta'
                                if i+1 <= length(cmdArgs) && isnumeric(str2double(cmdArgs{i+1})) && ~isnan(str2double(cmdArgs{i+1}))
                                    args.zeta = str2double(cmdArgs{i+1});
                                    i = i + 2;
                                else
                                    warning('无效的zeta参数，使用默认值%d', args.zeta);
                                    i = i + 1;
                                end
                            case '-subset'
                                if i+1 <= length(cmdArgs) && isnumeric(str2double(cmdArgs{i+1})) && ~isnan(str2double(cmdArgs{i+1}))
                                    args.subset = str2double(cmdArgs{i+1});
                                    i = i + 2;
                                else
                                    warning('无效的subset参数，使用默认值%d', args.subset);
                                    i = i + 1;
                                end
                            case '-force'
                                args.force = true;
                                i = i + 1;
                            case '-noconfig'
                                args.useConfig = false;
                                i = i + 1;
                            otherwise
                                warning('未知参数: %s', cmdArgs{i});
                                i = i + 1;
                        end
                    else
                        warning('忽略非字符串参数');
                        i = i + 1;
                    end
                end
            end
        end
    catch
        % 无法解析命令行参数，使用之前设置的值
    end
    
    % 验证参数
    if ~ismember(args.zeta, [8, 15, 31])
        warning('zeta值%d不在允许范围内，将重置为默认值%d', args.zeta, defaultZeta);
        args.zeta = defaultZeta;
    end
    
    % 显示参数使用说明
    fprintf('参数使用说明:\n');
    fprintf('  -zeta <值>   : 设置zeta参数 (当前值: %d)，可选值：[8, 15, 31]\n', args.zeta);
    fprintf('  -subset <值> : 设置subset索引 (当前值: %d)\n', args.subset);
    fprintf('  -force       : 强制重新计算，不使用缓存 (当前值: %d)\n', args.force);
    fprintf('  -noconfig    : 不使用配置文件管理搜索范围 (当前值: %d)\n', ~args.useConfig);
end