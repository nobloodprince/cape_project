clc; clearvars;

% Function definitions
f1 = @(y) y(2);
f2 = @(y) 1000 * (1 - y(1)^2) * y(2) - y(1);

% Initial values
y = [2, 0];
h = 1; % Step size
t_end = 300; % Final time
t_values = 0:h:t_end; % Time array

% Initialize solution arrays for plotting
y1_values = zeros(1, length(t_values));
y2_values = zeros(1, length(t_values));
y1_values(1) = y(1);
y2_values(1) = y(2);

fprintf('Before start:\n');
disp(y)

% Runge-Kutta 4th Order Method loop
for i = 2:length(t_values)
    t = t_values(i-1);
    fprintf('Iteration at t = %d\n', t);

    k1 = h * [f1(y), f2(y)];
    k2 = h * [f1(y + k1/2), f2(y + k1/2)];
    k3 = h * [f1(y + k2/2), f2(y + k2/2)];
    k4 = h * [f1(y + k3), f2(y + k3)];

    y = y + (k1 + 2*k2 + 2*k3 + k4) / 6;

    % Store values for plotting
    y1_values(i) = y(1);
    y2_values(i) = y(2);

    disp(y)
end

% Plot the results
figure;
plot(t_values, y1_values, 'b', 'LineWidth', 1.5); hold on;
plot(t_values, y2_values, 'r', 'LineWidth', 1.5);
xlabel('Time (t)');
ylabel('Solution');
title('Solution of the Van der Pol System using RK4');
legend('y_1(t)', 'y_2(t)');
grid on;
