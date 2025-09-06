clc, clearvars

% Given constants
D = 0.154; % Diameter (m)
P0 = 1601325; % Pressure at node 0 (Pa)
P5 = 101325; % Atmospheric pressure (Pa)

% Pipe lengths (m)
L01 = 100;
L12 = 300;
L23 = 300;
L45 = 300;
L13 = 1200;
L24 = 1200;
L34 = 1200;

% Other constants
f = 0.005;
p = 1000; % Density (kg/m^3)

% Pressure drop function (using absolute value for flow direction)
delta_P = @(L, qij) ((2*f*p*L)/((pi^2)*(D^5)))*(abs(qij))^2;

% Function defining the system of equations
F = @(v) [
    v(8)  - (P0 + delta_P(L01, v(1)));  % P1
    v(9)  - (v(8) + delta_P(L12, v(2)));  % P2
    v(10) - (v(9) + delta_P(L23, v(4)));  % P3
    v(10) - (v(8) + delta_P(L13, v(3)));  % P3 from a different path
    v(11) - (v(10) + delta_P(L34, v(6))); % P4
    v(11) - (v(9) + delta_P(L24, v(5)));  % P4 from a different path
    P5 - (v(11) + delta_P(L45, v(7))); % P5 = 101325 Pa

    % Flow conservation at nodes
    v(3) + v(2) - v(1);  % Node 1
    v(4) + v(5) - v(2);  % Node 2
    v(4) + v(3) - v(6);  % Node 3
    v(5) + v(6) - v(7)   % Node 4
];

% Initial Guess
vo = [1; 1; 1; 1; 1; 1; 1; 1.5e6; 1.3e6; 1.2e6; 1.1e6];

% Solve the system
solution = fsolve(F, vo);

% Display results
disp('Solution [q01, q12, q13, q23, q24, q34, q45, P1, P2, P3, P4]:');
disp(solution);
