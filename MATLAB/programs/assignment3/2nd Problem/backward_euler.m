clc; clearvars; close all;

% Define system of ODEs
f1 = @(y) y(2);
f2 = @(y) 1000 * (1 - y(1)^2) * y(2) - y(1);

% Initial conditions
y0 = [2; 0];
h = 0.01; % Larger step size possible due to stability
t_end = 3000;
t_steps = 0:h:t_end;
n = length(t_steps);

% Initialize solution array
y_bwd = zeros(2, n);
y_bwd(:, 1) = y0;

%**Backward Euler Method (Fixed-Point Approximation)**
for i = 1:n-1
    % Initial guess using Forward Euler
    y_guess = y_bwd(:, i) + h * [f1(y_bwd(:, i)); f2(y_bwd(:, i))];

    % Iterative correction (fixed-point iteration)
    for j = 1:5
        y_guess = y_bwd(:, i) + h * [f1(y_guess); f2(y_guess)];
    end

    y_bwd(:, i+1) = y_guess;
end

% **Plot Results**
figure;
subplot(1,2,1);
plot(t_steps, y_bwd(1,:), 'b', 'LineWidth', 1.2);
xlabel('Time'); ylabel('y_1 (Position)');
title('Backward Euler: y_1 vs Time');
grid on;

subplot(1,2,2);
plot(t_steps, y_bwd(2,:), 'r', 'LineWidth', 1.2);
xlabel('Time'); ylabel('y_2 (Velocity)');
title('Backward Euler: y_2 vs Time');
grid on;
