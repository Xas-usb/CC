%% USV 轨迹对比绘图脚本
% 对比期望路径 (Trajectory) 与实际路径 (Plant)

figure('Color', 'w', 'Name', 'USV 轨迹跟踪对比图');
hold on; grid on;

% 从工作区读取 Simulink 记录数据
% 注意：确保 Simulink 模型中 To Workspace 模块命名与此一致
x_ref = Ref_Data.Data(:,1);  
y_ref = Ref_Data.Data(:,2);
x_act = Act_Data.Data(:,1);
y_act = Act_Data.Data(:,2);

% 绘制期望路径 (红色虚线)
plot(x_ref, y_ref, 'r--', 'LineWidth', 2, 'DisplayName', '期望轨迹 (Desired)');

% 绘制实际路径 (蓝色实线)
plot(x_act, y_act, 'b-', 'LineWidth', 1.5, 'DisplayName', '实际轨迹 (Actual)');

% 标注起点与终点
plot(x_act(1), y_act(1), 'go', 'MarkerSize', 10, 'LineWidth', 2, ...
    'DisplayName', '实际起点');
plot(x_ref(1), y_ref(1), 'rx', 'MarkerSize', 10, 'LineWidth', 2, ...
    'DisplayName', '期望起点');

% 图形修饰
xlabel('X (m)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Y (m)', 'FontSize', 12, 'FontWeight', 'bold');
title('欠驱动无人艇轨迹跟踪对比图', 'FontSize', 14);
legend('Location', 'best', 'FontSize', 10);
axis equal;
set(gca, 'GridLineStyle', ':', 'GridColor', 'k', 'GridAlpha', 0.2);

% 末端跟踪误差
final_error = sqrt((x_ref(end)-x_act(end))^2 + (y_ref(end)-y_act(end))^2);
fprintf('仿真结束时的末端跟踪误差: %.4f 米\n', final_error);
