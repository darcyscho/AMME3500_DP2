% Extract time and data
time = out.nonlinearDissipation.Time;           % 52x1
raw_data = out.nonlinearDissipation.Data;       % 3x1x52

% Reshape to get 52x3
data = squeeze(raw_data)';  % size: 52x3

% Plot using stairs to replicate ZOH behavior
figure;
stairs(time, data(:,1), 'r', 'LineWidth', 1.5); hold on;
stairs(time, data(:,2), 'g', 'LineWidth', 1.5);
stairs(time, data(:,3), 'b', 'LineWidth', 1.5);
hold off;

xlabel('Time (s)');
ylabel('Nonlinear Dissipation (degrees)');
title('Nonlinear Dissipation per Joint');
legend('Joint 1', 'Joint 2', 'Joint 3');
grid on;
