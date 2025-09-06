clc; clearvars; close all;

% Define system of ODEs
f1 = @(y) y(2);
f2 = @(y) 1000 * (1 - y(1)^2) * y(2) - y(1);

% Initial conditions
y0 = [2; 0];
h = 0.001; % Small step size for stability
t_end = 3000;
t_steps = 0:h:t_end;
n = length(t_steps);

% Initialize solution array
y_fwd = zeros(2, n);
y_fwd(:, 1) = y0;

% **Forward Euler Method**
for i = 1:n-1
    y_fwd(:, i+1) = y_fwd(:, i) + h * [f1(y_fwd(:, i)); f2(y_fwd(:, i))];
end

% **Plot Results**
figure;
subplot(1,2,1);
plot(t_steps, y_fwd(1,:), 'b');
xlabel('Time'); ylabel('y_1 (Position)');
title('Forward Euler: y_1 vs Time');
grid on;

subplot(1,2,2);
plot(t_steps, y_fwd(2,:), 'r');
xlabel('Time'); ylabel('y_2 (Velocity)');
title('Forward Euler: y_2 vs Time');
grid on;

disp(y_fwd)