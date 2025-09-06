clc; clearvars; close all;

% Define system of ODEs
f1 = @(y) y(2);
f2 = @(y) 1000 * (1 - y(1)^2) * y(2) - y(1);

% Initial conditions
y0 = [2; 0];
t_span = [0 3000]; % Time span

% Define function handle for ODE solver
van_der_pol = @(t, y) [f1(y); f2(y)];

% Solve using ode15s
[t, y] = ode15s(van_der_pol, t_span, y0);

% **Plot Results**
figure;
subplot(1,2,1);
plot(t, y(:,1), 'b', 'LineWidth', 1.2);
xlabel('Time'); ylabel('y_1 (Position)');
title('ode15s: y_1 vs Time');
grid on;

subplot(1,2,2);
plot(t, y(:,2), 'r', 'LineWidth', 1.2);
xlabel('Time'); ylabel('y_2 (Velocity)');
title('ode15s: y_2 vs Time');
grid on;

disp(y)