clc; clear; close all;

tic;

%constants
T0 = 100;       
TL = 30;        
T_inf = 30;     
L = 2;          
beta = 1.5;     

% Define the system of first-order ODEs
odefun = @(x, y) [y(2); beta * (y(1) - T_inf)];

% Define boundary conditions
bcfun = @(ya, yb) [ya(1) - T0; yb(1) - TL];

% Initial mesh points
xmesh = linspace(0, L, 10); 

% Initial guess for the solution as a function handle
y_guess = @(x) [T0 + (TL - T0) * x / L; 0]; % Linear guess for T(x) and 0 slope

% Generate initial solution structure using bvpinit
solinit = bvpinit(xmesh, y_guess);

% Solve the boundary value problem
sol = bvp4c(odefun, bcfun, solinit);

% Extract refined solution
x = linspace(0, L, 100); % Refined x-mesh
y = deval(sol, x); % Evaluate the solution

ttt = toc;
% Plot the results
figure;
plot(x, y(1,:), 'go', 'LineWidth', 1.5);
xlabel('x (length of fin)');
ylabel('Temperature T(x)');
title('Temperature Distribution Using bvp4c');
grid on;

% Display results
disp('Temperature distribution along the fin:');
disp(table(x', y(1,:)', 'VariableNames', {'x', 'T'}));
