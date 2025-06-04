% --- Extract actual joint angles ---
joint_ts = out.joint_angle;
if ~isa(joint_ts, 'timeseries')
    error('"out.joint_angle" is not a timeseries object.');
end
t = joint_ts.Time;
y = joint_ts.Data;  % Nx3 matrix

% --- Extract reference joint angles ---
ref_ts = out.reference;
if ~isa(ref_ts, 'timeseries')
    error('"out.reference" is not a timeseries object.');
end
t_ref = ref_ts.Time;
y_ref = ref_ts.Data;  % Nx3 matrix

% --- Plot both actual and reference joint angles ---
figure;
plot(t, y(:,1), 'r', 'LineWidth', 1.5); hold on;
plot(t, y(:,2), 'g', 'LineWidth', 1.5);
plot(t, y(:,3), 'b', 'LineWidth', 1.5);

plot(t_ref, y_ref(:,1), '--r', 'LineWidth', 1.5);
plot(t_ref, y_ref(:,2), '--g', 'LineWidth', 1.5);
plot(t_ref, y_ref(:,3), '--b', 'LineWidth', 1.5);

% --- Labeling ---
xlabel('Time (s)');
ylabel('Joint Angle (rad)');
title('Joint Angles vs Reference');
legend('q1', 'q2', 'q3', 'q1 ref', 'q2 ref', 'q3 ref');
grid on;
