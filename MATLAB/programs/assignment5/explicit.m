function explicit_method
    clear; clc; close all;

    % Problem setup
    L = 1;        % Length of rod
    Nx = 50;      % Number of spatial points
    x = linspace(0, L, Nx);
    dx = x(2) - x(1);

    alpha_vals = [1, 10, 100]; % Diffusivity values
    t_final = 10;

    for alpha = alpha_vals
        dt = 0.4 * (dx^2 / alpha);  % CFL condition for stability
        Nt = ceil(t_final / dt);
        times = linspace(0, t_final, Nt);

        % Initial condition
        T = ones(Nx, 1) * 350;
        T(1) = 300;  % Dirichlet BC at x = 0
        T(end) = 400; % Dirichlet BC at x = L

        % Storage for visualization
        T_store = zeros(Nt, Nx);
        T_store(1, :) = T';

        % Explicit method loop
        r = alpha * dt / dx^2;
        fprintf('\n===== Results for α = %d =====\n', alpha);
        fprintf('%-10s %-10s\n', 'Time (s)', 'Temperature at different positions');

        for t_idx = 2:Nt
            T_new = T; % Copy previous temperature
            for i = 2:Nx-1
                T_new(i) = T(i) + r * (T(i+1) - 2*T(i) + T(i-1));
            end
            T = T_new;
            T_store(t_idx, :) = T';

            % Print vertical table for selected time steps
            if mod(t_idx, round(Nt/5)) == 0
                fprintf('\nTime = %.2f s\n', times(t_idx));
                fprintf('%-10s %-10s\n', 'Position (m)', 'Temperature (K)');
                for i = 1:Nx
                    fprintf('%-10.2f %-10.2f\n', x(i), T(i));
                end
            end
        end

        % Plot results
        figure; hold on;
        for t_idx = round(linspace(1, Nt, 5))
            plot(x, T_store(t_idx, :), 'DisplayName', sprintf('t = %.1f s', times(t_idx)));
        end
        hold off; grid on;
        xlabel('Position (m)'); ylabel('Temperature (K)');
        title(['Explicit Method: \alpha = ', num2str(alpha)]);
        legend show;
    end
end
