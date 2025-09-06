clc; clear; close all;

% Given parameters
F = 1;  V = 1;  k0 = 36e6;  UA = 150;
dH = -6500;  E = 12000;  Caf = 10;  
Tj0 = 298;  Fj = 1.25;  Vj = 0.25;
Cp = 500;  rho = 500;  Cpj = 600;  rhoj = 600;
Tf = 298; 

% Steady states from previous analysis
steady_states = [
    0.59, 405, 390;  % Steady state 1
    1.12, 355, 340;  % Steady state 2
    2.51, 320, 310   % Steady state 3
];

% Time settings
tspan = [0 50];
dt = 0.1;  % Step size for RK4

for i = 1:3
    % Select a steady state
    CA_ss = steady_states(i,1);
    T_ss = steady_states(i,2);
    Tj_ss = steady_states(i,3);

    perturbations = [0.01, 0.05, 0.25]; % 1%, 5%, 25%
    
    figure;
    hold on;
    
    for j = 1:length(perturbations)
        % Perturb the initial condition
        perturb = perturbations(j);
        y0 = [CA_ss * (1 + perturb), T_ss * (1 + perturb), Tj_ss * (1 + perturb)];
        
        % Solve using RK4
        [t, Y_rk4] = RK4(@(t, y) cstr_odes(t, y, F, V, k0, E, UA, dH, Caf, Tj0, Fj, Vj, Cp, rho, Cpj, rhoj, Tf), tspan, y0, dt);
        
        % Solve using ode45 for comparison
        [t_ode45, Y_ode45] = ode45(@(t, y) cstr_odes(t, y, F, V, k0, E, UA, dH, Caf, Tj0, Fj, Vj, Cp, rho, Cpj, rhoj, Tf), tspan, y0);
        
        % Plot results
        subplot(3,1,1);
        plot(t, Y_rk4(:,1), 'DisplayName', sprintf('Perturb %.1f%%', perturb*100));
        title(sprintf('Concentration C_A vs Time (Steady State %d)', i)); xlabel('Time'); ylabel('C_A');
        
        subplot(3,1,2);
        plot(t, Y_rk4(:,2), 'DisplayName', sprintf('Perturb %.1f%%', perturb*100));
        title('Reactor Temperature vs Time'); xlabel('Time'); ylabel('T');
        
        subplot(3,1,3);
        plot(t, Y_rk4(:,3), 'DisplayName', sprintf('Perturb %.1f%%', perturb*100));
        title('Jacket Temperature vs Time'); xlabel('Time'); ylabel('T_j');
    end
    
    legend;
end

% Function for ODEs
function dydt = cstr_odes(~, y, F, V, k0, E, UA, dH, Caf, Tj0, Fj, Vj, Cp, rho, Cpj, rhoj, Tf)
    CA = y(1); T = y(2); Tj = y(3);
    k = k0 * exp(-E / (1.987 * T)); 
    r = k * CA;

    dCA_dt = (F / V) * (Caf - CA) - r;
    dT_dt = (F / V) * (Tf - T) + (-dH / (rho * Cp)) * r - (UA / (rho * Cp * V)) * (T - Tj);
    dTj_dt = (Fj / Vj) * (Tj0 - Tj) + (UA / (rhoj * Cpj * Vj)) * (T - Tj);
    
    dydt = [dCA_dt; dT_dt; dTj_dt];
end

% RK4 implementation
function [t, Y] = RK4(odefun, tspan, y0, dt)
    t = tspan(1):dt:tspan(2);
    N = length(t);
    Y = zeros(N, length(y0));
    Y(1, :) = y0;
    
    for i = 1:N-1
        k1 = dt * odefun(t(i), Y(i, :)')';
        k2 = dt * odefun(t(i) + dt/2, (Y(i, :) + k1/2)')';
        k3 = dt * odefun(t(i) + dt/2, (Y(i, :) + k2/2)')';
        k4 = dt * odefun(t(i) + dt, (Y(i, :) + k3)')';
        Y(i+1, :) = Y(i, :) + (k1 + 2*k2 + 2*k3 + k4) / 6;
    end
end
