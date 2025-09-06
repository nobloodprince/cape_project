clear; clc; close all;

% Parameters
L = 1; Nx = 21;
x = linspace(0, L, Nx);
times = [1, 5, 10, 50, 100];
alpha_vals = [1, 10, 100];

for alpha = alpha_vals
    sol = pdepe(0, @(x,t,T,dTdx) heat_eqn(x,t,T,dTdx,alpha), @initial_cond, @boundary_cond, x, times);

    % Print results
    for i = 1:length(times)
        fprintf('pdepe: α = %d, Time = %d s, Temperature = %s\n', alpha, times(i), mat2str(round(sol(i,:),1)));
    end

    figure; hold on;
    for i = 1:length(times)
        plot(x, sol(i,:), 'o-', 'DisplayName', ['t = ', num2str(times(i)), ' s']);
    end
    title(['pdepe: \alpha = ', num2str(alpha)]);
    xlabel('Position (m)'); ylabel('Temperature (K)');
    legend; grid on;
    hold off;
end

% PDE function for pdepe
function [c, f, s] = heat_eqn(~,~,T,dTdx,alpha)
    c = 1 / alpha;
    f = dTdx;
    s = 0;
end

% Initial Condition for pdepe
function T0 = initial_cond(~)
    T0 = 350;
end

% Boundary Conditions for pdepe
function [pl, ql, pr, qr] = boundary_cond(xl, Tl, xr, Tr, t)
    pl = Tl - 300;
    ql = 0;
    pr = Tr - 400;
    qr = 0;
end
