% Time span for simulation
tspan = [0 10];  % seconds

% Equilibrium values
theta1_eq = 1.75;
tau1_eq = -18.53;

% Damping coefficient (adjust this value)
b1 = 1;  % Try values like 1.0 to 5.0

% Define ODE system: [x1 = theta1, x2 = dtheta1]
odefun = @(t, x) [
    x(2);  % dx1/dt = x2
    2.35*tau1_eq - 2.35*b1*x(2) - 1.33*x(1) + 45.88  % dx2/dt = ddtheta1
];

% Initial conditions: [theta1(0), dtheta1(0)]
x0 = [0; 0];  % Starting from rest, away from equilibrium

% Solve ODE
[t, x] = ode45(odefun, tspan, x0);

% Plot theta1(t)
figure;
plot(t, x(:,1), 'b', 'LineWidth', 2); hold on;
yline(theta1_eq, '--r', 'Equilibrium \theta_1', 'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'middle');
xlabel('Time (s)');
ylabel('\theta_1 (rad)');
title('System Trajectory Given Eqilibrium for Link 1');
grid on;
