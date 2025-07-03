% 测试脚本，运行 main.m 并捕获错误
try
    disp('开始运行 main.m...');
    main;
    disp('main.m 运行完成。');
catch e
    disp('运行 main.m 时出错:');
    disp(e.message);
    disp('错误堆栈:');
    disp(getReport(e));
end 