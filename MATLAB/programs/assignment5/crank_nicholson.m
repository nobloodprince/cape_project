clear; clc; close all;

% Parameters
L = 1; Nx = 21; dx = L / (Nx - 1);
x = linspace(0, L, Nx);
dt = 0.1; times = [1, 5, 10, 50, 100];
alpha_vals = [1, 10, 100];

for alpha = alpha_vals
    r = alpha * dt / (2 * dx^2);

    % Construct A and B Matrices
    A = diag((1 + 2 * r) * ones(Nx-2,1)) + diag(-r * ones(Nx-3,1), 1) + diag(-r * ones(Nx-3,1), -1);
    B = diag((1 - 2 * r) * ones(Nx-2,1)) + diag(r * ones(Nx-3,1), 1) + diag(r * ones(Nx-3,1), -1);

    % Initial and boundary conditions
    T = 350 * ones(Nx-2, 1);

    for n = 1:max(times)/dt
        % Right-hand side
        RHS = B * T;
        RHS(1) = RHS(1) + r * (300 + 300); % Adjust for left boundary
        RHS(end) = RHS(end) + r * (400 + 400); % Adjust for right boundary

        % Solve system
        T = A \ RHS;

        % Store and print results at specified times
        if ismember(n * dt, times)
            fprintf('Crank-Nicholson: α = %d, Time = %.1f s, Temperature = %s\n', alpha, n * dt, mat2str(round([300; T; 400]',1)));
        end
    end

    % Plot results
    figure; plot(x, [300; T; 400], 'o-'); title(['Crank-Nicholson: \alpha = ', num2str(alpha)]);
    xlabel('Position (m)'); ylabel('Temperature (K)'); grid on;
end
