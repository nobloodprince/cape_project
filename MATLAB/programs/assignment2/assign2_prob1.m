clc; clear; close all;

% Given Data
D = 0.154; % Pipe diameter (m)
rho = 1000; % Water density (kg/m^3)
f = 0.005; % Fanning friction factor
L = [100, 300, 1200, 300, 1200]; % Pipe lengths (m)
P0 = 15 * 1e5; % Pressure at pump exit (Pa)

% Compute k_ij for each pipe
k = (2 * f * rho .* L) ./ (pi^2 * D^5);

% Define the system of equations
fun = @(q) [
    P0 - k(1)*q(1)^2 - k(3)*q(3)^2 - k(2)*q(2)^2 - k(4)*q(4)^2 - k(5)*q(5)^2; % Node 1 pressure balance
    k(2)*q(2)^2 - k(3)*q(3)^2 - k(5)*q(5)^2; % Node 2 pressure balance
    k(3)*q(3)^2 - k(4)*q(4)^2; % Node 3 pressure balance
    k(4)*q(4)^2 - k(5)*q(5)^2; % Node 4 pressure balance
    q(1) - q(2) - q(3); % Flow continuity at node 1
    q(3) - q(4) - q(5); % Flow continuity at node 3
];

% Initial guess for flow rates (m^3/s)
q0 = ones(5,1); %

% Solve the nonlinear system
options = optimoptions('fsolve','Display','iter');
[q_sol, fval, exitflag] = fsolve(fun, q0, options);

%k values
fprintf('k values: \n k1 = %.2f\n k2 = %.2f\n k3 = %.2f\n k4 = %.2f\n k5 = %.2f\n \n', k(1), k(2), k(3), k(4), k(5));

% Display results
disp('Solved Flow Rates (m^3/s):');

%display
fprintf('Solution: \n q1 = %.6f\n q2 = %.6f\n q3 = %.6f\n q4 = %.6f\n q5 = %.6f\n \n', q_sol(1), q_sol(2), q_sol(3), q_sol(4), q_sol(5));

% Compute Pressures at Nodes
P = zeros(4,1);
P(1) = P0 - k(1)*q_sol(1)^2;
P(2) = P(1) - k(2)*q_sol(2)^2;
P(3) = P(2) - k(3)*q_sol(3)^2;
P(4) = P(3) - k(4)*q_sol(4)^2;

disp('Solved Pressures (Pa):');
fprintf('Solution: \n P1 = %.6f\n P2 = %.6f\n P3 = %.6f\n P4 = %.6f\n', P(1), P(2), P(3), P(4));
