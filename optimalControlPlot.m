% --- Extract control input timeseries ---
u_ts = out.control_inputs;

% Validate
if ~isa(u_ts, 'timeseries')
    error('"out.control_inputs" is not a timeseries object.');
end

% Extract time and data
t_u = u_ts.Time;
u = u_ts.Data;  % Nx3 matrix

% --- Plot control inputs ---
figure;
plot(t_u, u(:,1), 'm', 'LineWidth', 1.5); hold on;
plot(t_u, u(:,2), 'c', 'LineWidth', 1.5);
plot(t_u, u(:,3), 'k', 'LineWidth', 1.5);

% --- Labeling ---
xlabel('Time (s)');
ylabel('Control Input (Nm)');
title('Joint Control Inputs');
legend('u1 (hip)', 'u2 (knee)', 'u3 (ankle)');
grid on;
