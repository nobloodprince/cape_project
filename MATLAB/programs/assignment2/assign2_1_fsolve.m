 %contants
 clc, clearvars
 F = 1;
 V = 1;
 ko = 36*10e6;
 delta_H = 6500;
 E = 12000;
 p_Cp = 500;
 Tf = 298;
 CAf = 10;
 UA = 150;
 Tjo = 298;
 pj_Cj = 600;
 Fj = 1.25;
 Vj = 0.25;
 R = 8.314;

 % Define function for fsolve
 tic;
func = @(x) [
    F*(CAf - x(1)) - ko * exp(-E / (8.314 * x(2))) * x(1) * V;
    p_Cp * F * (Tf - x(2)) + (delta_H) * V * ko * exp(-E / (8.314 * x(2))) * x(1) - UA * (x(2) - x(3))
    pj_Cj * Fj * (Tjo - x(3)) + UA * (x(2) - x(3))
];

% Initial guess for [CA, T, Tj]
x0 = [5; 350; 800];

% Solve using fsolve
%options = optimoptions('fsolve', 'Display', 'iter'); % Show iterations
x = fsolve(func, x0);

t2 = toc;
% Display solution
%disp('Solution (Steady-state values for CA, T, Tj):');
%disp(x);

fprintf('FSOLVE FUNCTION\n');
fprintf('Solution: CA = %.6f, T = %.6f, Tj = %.6f\n', x(1), x(2), x(3));
fprintf('time taken(sec): %.6f\n', t2);
fprintf('Converged in %d iterations\n', false);
fprintf('\n');