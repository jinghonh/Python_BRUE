function [config, updated] = loadSaveConfig(mode, varargin)
    % 配置文件管理函数 - 用于加载和保存搜索范围配置
    %
    % 用法:
    %   config = loadSaveConfig('load', zeta, subset_index)
    %   updated = loadSaveConfig('save', zeta, subset_index, rangeMin, rangeMax, validFlow)
    %
    % 输入:
    %   mode - 'load'或'save'，指定操作模式
    %   zeta - zeta参数值
    %   subset_index - 子集索引
    %   rangeMin - (仅保存模式) 最小范围向量
    %   rangeMax - (仅保存模式) 最大范围向量
    %   validFlow - (仅保存模式) 有效流量，用于自适应更新搜索范围
    %
    % 输出:
    %   config - (加载模式) 包含rangeMin和rangeMax的配置结构体
    %   updated - (保存模式) 布尔值，指示配置是否已更新
    
    % 配置文件路径
    configDir = 'config';
    if ~exist(configDir, 'dir')
        mkdir(configDir);
    end
    
    % 确保配置文件夹存在
    configFile = fullfile(configDir, 'search_ranges.mat');
    
    % 默认配置
    defaultConfig = struct();
    
    % zeta=8 默认范围
    defaultConfig.zeta8.subset0.rangeMin = [0, 0, 0];
    defaultConfig.zeta8.subset0.rangeMax = [10000, 10000, 10000];
    
    % zeta=15 默认范围
    defaultConfig.zeta15.subset0.rangeMin = [3201, 3151, 1201];
    defaultConfig.zeta15.subset0.rangeMax = [5002, 4800, 2400];
    defaultConfig.zeta15.subset1.rangeMin = [4200, 3000, 0, 0, 0];
    defaultConfig.zeta15.subset1.rangeMax = [4900, 4200, 1610, 1610, 2000];
    
    % zeta=31 默认范围
    defaultConfig.zeta31.subset0.rangeMin = [1588, 2800, 1198];
    defaultConfig.zeta31.subset0.rangeMax = [4903, 5823, 2840];
    defaultConfig.zeta31.subset1.rangeMin = [500, 1200, 0, 0, 0];
    defaultConfig.zeta31.subset1.rangeMax = [4899, 6001, 2507, 2503, 2717];
    defaultConfig.zeta31.subset2.rangeMin = [4000, 3000, 0, 0, 902, 0];
    defaultConfig.zeta31.subset2.rangeMax = [5007, 4201, 3909, 3909, 1353, 1601];
    
    switch lower(mode)
        case 'load'
            % 加载模式
            if nargin < 3
                error('加载模式需要指定zeta和subset_index参数');
            end
            
            zeta = varargin{1};
            subset_index = varargin{2};
            
            % 验证参数
            if ~ismember(zeta, [8, 15, 31])
                error('zeta必须为8、15或31');
            end
            
            % 加载现有配置文件或使用默认值
            if exist(configFile, 'file')
                loadedConfig = load(configFile);
                config = loadedConfig.config;
            else
                config = defaultConfig;
            end
            
            % 获取指定zeta和subset的配置
            zetaField = ['zeta' num2str(zeta)];
            subsetField = ['subset' num2str(subset_index)];
            
            % 检查配置是否存在
            if isfield(config, zetaField) && isfield(config.(zetaField), subsetField)
                % 返回指定配置
                rangeMin = config.(zetaField).(subsetField).rangeMin;
                rangeMax = config.(zetaField).(subsetField).rangeMax;
                
                fprintf('从配置文件加载搜索范围: zeta=%d, subset=%d\n', zeta, subset_index);
                fprintf('rangeMin = [%s]\n', num2str(rangeMin));
                fprintf('rangeMax = [%s]\n', num2str(rangeMax));
            else
                % 使用默认配置
                if isfield(defaultConfig, zetaField) && isfield(defaultConfig.(zetaField), subsetField)
                    rangeMin = defaultConfig.(zetaField).(subsetField).rangeMin;
                    rangeMax = defaultConfig.(zetaField).(subsetField).rangeMax;
                    
                    % 更新配置文件
                    config.(zetaField).(subsetField).rangeMin = rangeMin;
                    config.(zetaField).(subsetField).rangeMax = rangeMax;
                    save(configFile, 'config');
                    
                    fprintf('使用默认搜索范围并保存到配置: zeta=%d, subset=%d\n', zeta, subset_index);
                else
                    error('没有找到匹配的默认配置: zeta=%d, subset=%d', zeta, subset_index);
                end
            end
            
            % 返回配置结构体
            config = struct('rangeMin', rangeMin, 'rangeMax', rangeMax);
            
        case 'save'
            % 保存模式
            if nargin < 6
                error('保存模式需要指定zeta、subset_index、rangeMin、rangeMax和validFlow参数');
            end
            
            zeta = varargin{1};
            subset_index = varargin{2};
            rangeMin = varargin{3};
            rangeMax = varargin{4};
            validFlow = varargin{5};
            
            % 加载现有配置文件或使用默认值
            if exist(configFile, 'file')
                loadedConfig = load(configFile);
                config = loadedConfig.config;
            else
                config = defaultConfig;
            end
            
            % 获取指定zeta和subset的字段名
            zetaField = ['zeta' num2str(zeta)];
            subsetField = ['subset' num2str(subset_index)];
            
            % 如果有有效样本，优化搜索范围
            if ~isempty(validFlow)
                % 计算有效样本的统计信息
                minFlow = min(validFlow);
                maxFlow = max(validFlow);
                
                % 扩展边界，给予一定余量（10%）
                range = maxFlow - minFlow;
                minFlow = max(0, minFlow - range * 0.1);
                maxFlow = maxFlow + range * 0.1;
                
                % 更新配置
                if ~isfield(config, zetaField)
                    config.(zetaField) = struct();
                end
                
                % 更新配置或创建新配置
                config.(zetaField).(subsetField).rangeMin = minFlow;
                config.(zetaField).(subsetField).rangeMax = maxFlow;
                config.(zetaField).(subsetField).lastUpdate = datestr(now);
                
                % 保存配置
                save(configFile, 'config');
                updated = true;
                
                fprintf('已更新配置文件的搜索范围: zeta=%d, subset=%d\n', zeta, subset_index);
                fprintf('新rangeMin = [%s]\n', num2str(minFlow));
                fprintf('新rangeMax = [%s]\n', num2str(maxFlow));
            else
                % 没有有效样本，保存当前范围
                if ~isfield(config, zetaField)
                    config.(zetaField) = struct();
                end
                
                % 更新配置
                config.(zetaField).(subsetField).rangeMin = rangeMin;
                config.(zetaField).(subsetField).rangeMax = rangeMax;
                config.(zetaField).(subsetField).lastUpdate = datestr(now);
                
                % 保存配置
                save(configFile, 'config');
                updated = true;
                
                fprintf('已保存当前搜索范围到配置文件: zeta=%d, subset=%d\n', zeta, subset_index);
            end
            
        otherwise
            error('未知的模式: %s，应为"load"或"save"', mode);
    end
end 