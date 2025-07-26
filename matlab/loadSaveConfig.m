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
    
    % 配置文件路径 - 使用当前目录的相对路径
    currentDir = pwd;
    configDir = fullfile(currentDir, 'config');
    
    % 初始化返回值
    updated = false;
    
    % 确保配置文件夹存在
    try
        if ~exist(configDir, 'dir')
            fprintf('正在创建配置目录: %s\n', configDir);
            [status, msg] = mkdir(configDir);
            if ~status
                error('创建配置目录失败: %s', msg);
            end
        end
        configFile = fullfile(configDir, 'search_ranges.mat');
    catch e
        % 如果无法创建config目录，则使用cache目录作为备选
        warning('无法创建config目录: %s，将使用cache目录代替', sprintf('%s', e.message));
        configDir = fullfile(currentDir, 'cache');
        if ~exist(configDir, 'dir')
            [status, msg] = mkdir(configDir);
            if ~status
                error('创建缓存目录失败: %s', msg);
            end
        end
        configFile = fullfile(configDir, 'search_ranges.mat');
    end
    
    % 默认配置
    defaultConfig = struct();
    
    % zeta <= 9 默认范围
    defaultConfig.zeta_low = struct();
    defaultConfig.zeta_low.subset0.rangeMin = [0, 0, 0];
    defaultConfig.zeta_low.subset0.rangeMax = [10000, 10000, 10000];
    
    % 10 <= zeta <= 26 默认范围
    defaultConfig.zeta_mid = struct();
    defaultConfig.zeta_mid.subset0.rangeMin = [3201, 3151, 1201];
    defaultConfig.zeta_mid.subset0.rangeMax = [5002, 4800, 2400];
    defaultConfig.zeta_mid.subset1.rangeMin = [4200, 3000, 0, 0, 0];
    defaultConfig.zeta_mid.subset1.rangeMax = [4900, 4200, 1610, 1610, 2000];
    
    % zeta > 26 默认范围
    defaultConfig.zeta_high = struct();
    defaultConfig.zeta_high.subset0.rangeMin = [1588, 2800, 1198];
    defaultConfig.zeta_high.subset0.rangeMax = [4903, 5823, 2840];
    defaultConfig.zeta_high.subset1.rangeMin = [500, 1200, 0, 0, 0];
    defaultConfig.zeta_high.subset1.rangeMax = [4899, 6001, 2507, 2503, 2717];
    defaultConfig.zeta_high.subset2.rangeMin = [4000, 3000, 0, 0, 902, 0];
    defaultConfig.zeta_high.subset2.rangeMax = [5007, 4201, 3909, 3909, 1353, 1601];
    
    switch lower(mode)
        case 'load'
            % 加载模式
            if nargin < 3
                error('加载模式需要指定zeta和subset_index参数');
            end
            
            zeta = varargin{1};
            subset_index = varargin{2};
            
            % 验证参数
            if zeta <= 0
                error('zeta必须为正数');
            end
            
            % 根据zeta值范围验证subset_index
            if zeta <= 9 && ~ismember(subset_index, 0)
                error('当zeta <= 9时，subset_index必须为0');
            elseif zeta >= 10 && zeta <= 26 && ~ismember(subset_index, [0, 1])
                error('当10 <= zeta <= 26时，subset_index必须为0或1');
            elseif zeta > 26 && ~ismember(subset_index, [0, 1, 2])
                error('当zeta > 26时，subset_index必须为0、1或2');
            end
            
            % 加载现有配置文件或使用默认值
            try
                if exist(configFile, 'file')
                    loadedConfig = load(configFile);
                    config = loadedConfig.config;
                else
                    config = defaultConfig;
                end
            catch e
                warning('加载配置文件失败: %s，将使用默认配置', sprintf('%s', e.message));
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
                if zeta <= 9
                    zetaCategory = 'zeta_low';
                elseif zeta >= 10 && zeta <= 26
                    zetaCategory = 'zeta_mid';
                else
                    zetaCategory = 'zeta_high';
                end
                
                if isfield(defaultConfig, zetaCategory) && isfield(defaultConfig.(zetaCategory), subsetField)
                    rangeMin = defaultConfig.(zetaCategory).(subsetField).rangeMin;
                    rangeMax = defaultConfig.(zetaCategory).(subsetField).rangeMax;
                    
                    % 尝试更新配置文件
                    try
                        % 确保zeta字段存在
                        if ~isfield(config, zetaField)
                            config.(zetaField) = struct();
                        end
                        
                        % 更新配置并保存
                        config.(zetaField).(subsetField).rangeMin = rangeMin;
                        config.(zetaField).(subsetField).rangeMax = rangeMax;
                        save(configFile, 'config');
                        fprintf('使用默认搜索范围并保存到配置: zeta=%d, subset=%d\n', zeta, subset_index);
                    catch e
                        warning('保存默认配置失败: %s，将只使用内存中的配置', sprintf('%s', e.message));
                    end
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
            try
                if exist(configFile, 'file')
                    loadedConfig = load(configFile);
                    config = loadedConfig.config;
                else
                    config = defaultConfig;
                    fprintf('未找到配置文件，将创建新配置\n');
                end
            catch e
                warning('加载配置文件失败: %s，将使用默认配置', sprintf('%s', e.message));
                config = defaultConfig;
            end
            
            % 获取指定zeta和subset的字段名
            zetaField = ['zeta' num2str(zeta)];
            subsetField = ['subset' num2str(subset_index)];
            
            % 确保zeta字段存在
            if ~isfield(config, zetaField)
                config.(zetaField) = struct();
            end
            
            % 如果有有效样本，优化搜索范围
            if ~isempty(validFlow)
                % 计算有效样本的统计信息
                minFlow = min(validFlow);
                maxFlow = max(validFlow);
                
                % 扩展边界，给予一定余量（10%）
                range = maxFlow - minFlow;
                minFlow = max(0, minFlow - range * 0.1);
                maxFlow = maxFlow + range * 0.1;
                
                % 更新配置或创建新配置
                config.(zetaField).(subsetField).rangeMin = minFlow;
                config.(zetaField).(subsetField).rangeMax = maxFlow;
                config.(zetaField).(subsetField).lastUpdate = datestr(now);
                
                % 尝试保存配置
                try
                    % 保存配置
                    fprintf('正在保存配置到: %s\n', configFile);
                    save(configFile, 'config');
                    updated = true;
                    fprintf('已更新配置文件的搜索范围: zeta=%d, subset=%d\n', zeta, subset_index);
                    fprintf('新rangeMin = [%s]\n', num2str(minFlow));
                    fprintf('新rangeMax = [%s]\n', num2str(maxFlow));
                catch e
                    warning('保存配置文件失败: %s\n尝试使用临时文件...', sprintf('%s', e.message));
                    % 尝试保存到临时文件
                    try
                        tempFile = fullfile(tempdir, 'search_ranges_temp.mat');
                        save(tempFile, 'config');
                        warning('已将配置保存到临时文件: %s', tempFile);
                    catch e2
                        warning('保存到临时文件也失败: %s', sprintf('%s', e2.message));
                    end
                end
            else
                % 没有有效样本，保存当前范围
                config.(zetaField).(subsetField).rangeMin = rangeMin;
                config.(zetaField).(subsetField).rangeMax = rangeMax;
                config.(zetaField).(subsetField).lastUpdate = datestr(now);
                
                % 尝试保存配置
                try
                    % 保存配置
                    save(configFile, 'config');
                    updated = true;
                    fprintf('已保存当前搜索范围到配置文件: zeta=%d, subset=%d\n', zeta, subset_index);
                catch e
                    warning('保存配置文件失败: %s', sprintf('%s', e.message));
                end
            end
            
        otherwise
            error('未知的模式: %s，应为"load"或"save"', mode);
    end
end 