clc; clear; close all;

%% Given Data
D = 0.154; % Pipe diameter (m)
rho = 1000; % Water density (kg/m^3)
f = 0.005; % Fanning friction factor
L = [100, 300, 1200, 300, 1200]; % Pipe lengths (m)
P0 = 15 * 1e5; % Pressure at pump exit (Pa)

%% Compute k_ij for each pipe
k = (2 * f * rho .* L) ./ (pi^2 * D^5);

%% Modified System of Equations for Blocked Pipe (q_4 = 0)
fun_blocked = @(q) [
    P0 - k(1)*q(1)^2 - k(3)*q(3)^2 - k(2)*q(2)^2 - k(5)*q(4)^2; % Node 1 pressure balance
    k(2)*q(2)^2 - k(3)*q(3)^2 - k(5)*q(4)^2; % Node 2 pressure balance
    k(3)*q(3)^2 - k(5)*q(4)^2; % Node 3 pressure balance
    q(1) - q(2) - q(3); % Flow continuity at node 1
    q(3) - q(4); % Flow continuity at node 3 (q4 = 0)
];

%% Initial Guess for Flow Rates
q0_blocked = ones(4,1); % Now only four unknowns

%% Solve the Nonlinear System for Blocked Pipeline
options = optimoptions('fsolve', 'Display', 'iter');
[q_sol_blocked, fval_blocked, exitflag_blocked] = fsolve(fun_blocked, q0_blocked, options);

% Compute Pressures at Nodes for Blocked Condition
P_blocked = zeros(4,1);
P_blocked(1) = P0 - k(1)*q_sol_blocked(1)^2;
P_blocked(2) = P_blocked(1) - k(2)*q_sol_blocked(2)^2;
P_blocked(3) = P_blocked(2) - k(3)*q_sol_blocked(3)^2;
P_blocked(4) = P_blocked(3); % Since q4 = 0, P3 = P4

% Display Results for Blocked Pipeline Case
fprintf('--- Blocked Pipeline Results ---\n');
fprintf('Flow Rates (m^3/s) with Pipe 2-4 Blocked:\n');
for i = 1:length(q_sol_blocked)
    fprintf('q_%d = %.6f\n', i, q_sol_blocked(i));
end

fprintf('\nPressures (Pa) with Pipe 2-4 Blocked:\n');
for i = 1:length(P_blocked)
    fprintf('P_%d = %.2f Pa\n', i, P_blocked(i));
end


