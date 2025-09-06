clearvars, clc
%Constants
R = 0.08206;  % L·atm·mol⁻¹·K⁻¹
Tc = 407.5;   % K
Pc = 111.3;   % atm

%Given Conditions
T = 250 + 273.15;  % Convert °C to K
P = 10;       % atm

% Calculate a and b
a = (27 * R^2 * Tc^2) / (64 * Pc);
b = (R * Tc) / (8 * Pc);

% function definition as f(v) = 0
f = @(v) (P + a/v^2) * (v - b) - R * T;

% Derivative of f(v) for Newton's method
f_prime = @(v) ((-2) * a * (v - b)/ v^3 + (P + a/v^2));

% Fixed-point iteration
function [v, iterations] = fixed_point_iteration(f, tol, max_iter, R, T, P, a, b)
    g = @(v) (R * T) ./ (P + a ./ v.^2) + b;
    v = 1.0;  % Initial guess
    for i = 1:max_iter
        v_new = g(v);
        if abs(v_new - v) < tol
            iterations = i;
            return;
        end
        v = v_new;
    end
    error('Fixed-point iteration did not converge');
end

% Bisection method
function [v, iterations] = bisection_method(f, tol, max_iter)
    v_low = 0.01;
    v_high = 10;  % Initial bracket (needs adjustment if necessary)
    if f(v_low) * f(v_high) > 0
        error('Bisection method requires a sign change in the interval');
    end

    for i = 1:max_iter
        v_mid = (v_low + v_high) / 2;
        if abs(f(v_mid)) < tol
            v = v_mid;
            iterations = i;
            return;
        elseif f(v_low) * f(v_mid) < 0
            v_high = v_mid;
        else
            v_low = v_mid;
        end
    end
    error('Bisection method did not converge');
end

% Newton's method
function [v, iterations] = newtons_method(f, f_prime, tol, max_iter)
    v = 1.0;  % Initial guess
    for i = 1:max_iter
        v_new = v - f(v) / f_prime(v);
        if abs(v_new - v) < tol
            iterations = i;
            return;
        end
        v = v_new;
    end
    error('Newton''s method did not converge');
end

% MATLAB fzero
function v = matlab_fzero(f)
    v = fzero(f, 1.0);  % Initial guess
end

% Compare methods
function compare_methods()
    % Constants
    R = 0.08206;  % L·atm·mol⁻¹·K⁻¹
    Tc = 407.5;   % K
    Pc = 111.3;   % atm
    T = 250 + 273.15;  % Convert °C to K
    P = 10;       % atm

    % Calculate a and b
    a = (27 * R^2 * Tc^2) / (64 * Pc);
    b = (R * Tc) / (8 * Pc);

    % Van der Waals equation as f(v) = 0
    f = @(v) (P + a ./ v.^2) .* (v - b) - R * T;

    % Derivative of f(v) for Newton's method
    f_prime = @(v) -2 * a * (v - b) ./ v.^3 + (P + a ./ v.^2);

    tol = 1e-6;
    max_iter = 1000;

    methods = {
        'Fixed-point Iteration', @() fixed_point_iteration(f, tol, max_iter, R, T, P, a, b);
        'Bisection Method', @() bisection_method(f, tol, max_iter);
        'Newton''s Method', @() newtons_method(f, f_prime, tol, max_iter);
        'MATLAB fzero', @() matlab_fzero(f)
    };

    for i = 1:size(methods, 1)
        method_name = methods{i, 1};
        method_func = methods{i, 2};

        try
            tic;
            if strcmp(method_name, 'MATLAB fzero')
                v = method_func();
                iterations = NaN;
            else
                [v, iterations] = method_func();
            end
            time_taken = toc;

            fprintf('\n%s:\n', method_name);
            fprintf('  Solution: %.6f\n', v);
            fprintf('  Iterations: %d\n', iterations);
            fprintf('  Time (s): %.6f\n', time_taken);
        catch ME
            fprintf('\n%s:\n', method_name);
            fprintf('  Error: %s\n', ME.message);
        end
    end
end

% Run comparison
compare_methods();
