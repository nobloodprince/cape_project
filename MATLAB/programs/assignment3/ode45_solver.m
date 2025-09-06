clc; clear; close all;

% Time settings
tspan = [0 50];

% Steady-state values
steady_states = [0.59 405 390; 1.12 355 340; 2.51 320 310];

% Loop through each steady state
for i = 1:3
    y0 = steady_states(i, :); % Initial condition

    % Solve using ode45
    [t, y] = ode45(@reactor_odes, tspan, y0);

    % Plot results
    figure;
    plot(t, y(:, 1), 'r', t, y(:, 2), 'b', t, y(:, 3), 'g');
    legend('C_A', 'T', 'T_j');
    title(['ODE45 Solution for Steady State ', num2str(i)]);
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
    
    dydt = [dCA_dt; dT_dt; dTj_dt];
end
