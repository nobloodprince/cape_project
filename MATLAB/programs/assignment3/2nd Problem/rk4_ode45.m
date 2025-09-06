clc, clearvars
ode_sys = @(t, y) [y(2);
    1000*(1-y(1)^2)*y(2)-y(1)];
tspan = [0,3000];
y0 = [2 0];

[t, y] = ode45(ode_sys, tspan, y0);

figure;
plot(t, y(:, 1), 'r*')
hold on;
plot(t, y(:, 2), 'b+')
xlabel('Time (t)');
ylabel('Solution');

disp(y)