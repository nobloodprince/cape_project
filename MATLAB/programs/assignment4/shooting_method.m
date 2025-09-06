clc; clear; close all;

tic;

% constants
T0 = 100;       
TL = 30;        
T_inf = 30;     
L = 2;          
beta = 1.5;     

% Define ODE as a system of first-order equations
odefun = @(x, y) [y(2); beta * (y(1) - T_inf)];

% Define shooting function to match boundary condition at x = L
shootingFunc = @(s) ode45(odefun, [0, L], [T0; s]); % Solve ODE
residualFunc = @(s) shootingFunc(s).y(1, end) - TL; % Residual at x = L

% Solve for correct initial slope using fzero
s_correct = fzero(residualFunc, 0);

% Solve ODE with correct slope
[x, y] = ode45(odefun, linspace(0, L, 100), [T0; s_correct]);

tt = toc;


% Plot the results
figure;
plot(x, y(:,1), 'r*', 'LineWidth', 1.5);
xlabel('x (length of fin)');
ylabel('Temperature T(x)');
title('Temperature Distribution Using Shooting Method');
grid on;

% Display results
disp('Temperature distribution along the fin:');
disp(table(x, y(:,1), 'VariableNames', {'x', 'T'}));
