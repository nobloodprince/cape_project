clc; clear; close all;

% Given Data
D = 0.154; % Pipe diameter (m)
rho = 1000; % Water density (kg/m^3)
f = 0.005; % Fanning friction factor
L = [100, 300, 1200, 300, 1200, 1200, 1200]; % Pipe lengths (m)
P0 = 1601325; % Initial pressure (Pa)
P5 = 101325;  % Outlet pressure (Pa)

% Compute k_ij for each pipe
k = (2 * f * rho .* L) ./ (pi^2 * D^5);

% Define nonlinear system of equations
fun = @(q) [
    P0 - k(1)*q(1)^2 - k(2)*q(2)^2 - k(3)*q(3)^2; % Node 1 pressure balance
    k(2)*q(2)^2 - k(4)*q(4)^2; % Node 2 pressure balance
    k(3)*q(3)^2 - k(5)*q(5)^2; % Node 3 pressure balance
    k(4)*q(4)^2 - k(6)*q(6)^2; % Node 4 pressure balance
    k(5)*q(5)^2 - k(7)*q(7)^2; % Node 5 pressure balance
    P5 - (P0 - (k(1)*q(1)^2 + k(2)*q(2)^2 + k(3)*q(3)^2 + k(4)*q(4)^2 + k(5)*q(5)^2 + k(6)*q(6)^2 + k(7)*q(7)^2)); % Final outlet pressure check
    
    q(1) - q(2) - q(3); % Flow continuity at node 1
    q(2) - q(4); % Flow continuity at node 2
    q(3) - q(5); % Flow continuity at node 3
    q(4) - q(6); % Flow continuity at node 4
    q(5) - q(7); % Flow continuity at node 5
];

% Initial guess for flow rates (m^3/s)
q0 = ones(7,1) * 0.1; % Reasonable initial guess

% Solve the nonlinear system
options = optimoptions('fsolve', 'Display', 'iter', 'FunctionTolerance', 1e-6);
[q_sol, fval, exitflag] = fsolve(fun, q0, options);

% Display solved flow rates
fprintf('Flow Rates (m^3/s):\n');
for i = 1:7
    fprintf('q%d = %.6f\n', i, q_sol(i));
end

% Compute Pressures at Nodes
P = zeros(6,1);
P(1) = P0 - k(1)*q_sol(1)^2;
P(2) = P(1) - k(2)*q_sol(2)^2;
P(3) = P(1) - k(3)*q_sol(3)^2;
P(4) = P(2) - k(4)*q_sol(4)^2;
P(5) = P(3) - k(5)*q_sol(5)^2;
P(6) = P(5) - k(7)*q_sol(7)^2; % Should match P5

% Display pressures
fprintf('\nPressures (Pa):\n');
for i = 1:6
    fprintf('P%d = %.2f\n', i, P(i));
end

% Check if P6 matches P5
fprintf('\nCheck: P6 = %.2f, P5 (Expected) = %.2f\n', P(6), P5);
