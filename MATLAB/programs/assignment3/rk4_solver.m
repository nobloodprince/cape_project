clc; clear; close all;

% Time settings
tspan = [0 50]; % Time range
h = 0.1; % Step size
t = tspan(1):h:tspan(2);
n = length(t);

% Steady-state values
steady_states = [0.59 405 390; 1.12 355 340; 2.51 320 310];

% Loop through each steady state
for i = 1:3
    y0 = steady_states(i, :)'; % Convert to column vector
    y = zeros(3, n); % Preallocate (3 states × n time points)
    y(:, 1) = y0; % Set initial condition

    for j = 1:n-1
        k1 = h * reactor_odes(t(j), y(:, j));
        k2 = h * reactor_odes(t(j) + h/2, y(:, j) + k1/2);
        k3 = h * reactor_odes(t(j) + h/2, y(:, j) + k2/2);
        k4 = h * reactor_odes(t(j) + h, y(:, j) + k3);
        
        y(:, j+1) = y(:, j) + (k1 + 2*k2 + 2*k3 + k4) / 6;
    end

    % Plot results
    figure;
    plot(t, y(1, :), 'r', t, y(2, :), 'b', t, y(3, :), 'g');
    legend('C_A', 'T', 'T_j');
    title(['RK4 Solution for Steady State ', num2str(i)]);
    xlabel('Time');
    ylabel('State Variables');
end

function dydt = reactor_odes(~, y)
    % Define reactor ODEs (Replace with actual equations)
    CA = y(1);
    T = y(2);
    Tj = y(3);
    
    % Example ODEs (Modify based on system)
    dCA_dt = -0.1 * CA;
    dT_dt = 0.05 * (Tj - T);
    dTj_dt = 0.02 * (T - Tj);
    
    dydt = [dCA_dt; dT_dt; dTj_dt]; % Ensure column vector output
end
