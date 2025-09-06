clc; clear; close all;
tic;
% Given parameters
T0 = 100;       % Temperature at x = 0
TL = 30;        % Temperature at x = L
T_inf = 30;     % Ambient temperature
L = 2;          % Length of the fin
beta = 1.5;     % Coefficient

N = 10;         % Number of grid points
dx = L / (N - 1); % Step size
x = linspace(0, L, N); % Discretized x domain

% Coefficients for finite difference scheme
A = zeros(N, N);
b = zeros(N, 1);

% Boundary conditions
A(1,1) = 1;
b(1) = T0; % T(0) = T0

A(N,N) = 1;
b(N) = TL; % T(L) = TL

% Finite Difference Method for interior points
for i = 2:N-1
    A(i, i-1) = 1/dx^2;
    A(i, i)   = -2/dx^2 - beta;
    A(i, i+1) = 1/dx^2;
    b(i) = -beta * T_inf;
end

% Solve the linear system
T = A \ b;
t = toc;
disp(table(x', T, 'VariableNames', {'x', 'T'}));

% Plot the results
figure;
plot(x, T, '-o', 'LineWidth', 1.5);
xlabel('x (length of fin)');
ylabel('Temperature T(x)');
title('Temperature Distribution Along the Fin');
grid on;
