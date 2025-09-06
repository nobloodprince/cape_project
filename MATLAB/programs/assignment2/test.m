clc; clear; close all;
syms x y z

% Define functions symbolically
F_sym = [ x^2 + y^2 + z^2 - 9;
          x * exp(y) + z - 3;
          x + y + z - 1];

% Compute Jacobian matrix
J_sym = jacobian(F_sym, [x, y, z]);

% Convert symbolic functions to numerical functions
F_func = matlabFunction(F_sym, 'Vars', {x, y, z});
J_func = matlabFunction(J_sym, 'Vars', {x, y, z});

% Initial guess
X = [1; 1; 1];  
tol = 1e-6;  
max_iter = 50;  

for k = 1:max_iter
    % Evaluate function and Jacobian numerically
    F_val = F_func(X(1), X(2), X(3));
    J_val = J_func(X(1), X(2), X(3));
    
    % Solve for update step: J(X) * dX = -F(X)
    dX = -J_val \ F_val;
    
    % Update solution
    X = X + dX;
    
    % Check convergence
    if norm(dX, inf) < tol
        fprintf('Converged in %d iterations\n', k);
        break;
    end
end

% Display results
fprintf('Solution: x = %.6f, y = %.6f, z = %.6f\n', X(1), X(2), X(3));
