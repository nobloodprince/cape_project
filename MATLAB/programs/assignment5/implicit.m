clear; clc; close all;

% Parameters
L = 1; Nx = 21; dx = L / (Nx - 1);
x = linspace(0, L, Nx);
dt = 0.001; times = [1, 5, 10, 50, 100];
alpha_vals = [1, 10, 100];

for alpha = alpha_vals
    r = alpha * dt / dx^2;
    A = (1 + 2*r) * eye(Nx-2) - diag(r * ones(Nx-3, 1), 1) - diag(r * ones(Nx-3, 1), -1);
    T = 350 * ones(Nx-2, 1);

    for n = 1:max(times)/dt
        T = A \ (T + r * [300; zeros(Nx-4,1); 400]);

        if ismember(n * dt, times)
            fprintf('Implicit: α = %d, Time = %.1f s, Temperature = %s\n', alpha, n * dt, mat2str(round([300; T; 400]',1)));
        end
    end

    figure; plot(x, [300; T; 400], 'o-'); title(['Implicit: \alpha = ', num2str(alpha)]);
    xlabel('Position (m)'); ylabel('Temperature (K)'); grid on;
end
