clc; clear; close all;

% Given Parameters
F = 1;          % m^3/h
V = 1;          % m^3
k0 = 36e6;      % h^-1
E = 12000;      % kcal/kmol
UA = 150;       % kcal/°C h
Tj0 = 298;      % K
C_Af = 10;      % kmol/m^3
Cp = 500;       % kcal/m^3°C
rho = 500;      % kg/m^3
Fj = 1.25;      % m^3/h
Vj = 0.25;      % m^3
rho_j = 600;    % kg/m^3
Cp_j = 600;     % kcal/m^3°C
H = 6500;       % kcal/kmol
T_f = 298;      % K

% Steady States (Replace with your actual values)
steady_states = [
    0.59, 405, 390;  % Unstable high temp
    1.12, 355, 340;  % Stable intermediate temp
    2.51, 320, 310   % Stable low temp
];

perturbations = [0.01, 0.05, 0.25];  % 1%, 5%, 25% perturbations

for i = 1:3 % Loop through steady states
    C_A_ss = steady_states(i,1);
    T_ss = steady_states(i,2);
    Tj_ss = steady_states(i,3);
    
    figure;
    hold on;
    title(['Steady State ', num2str(i)]);
    xlabel('Time (h)');
    ylabel('State Variables');
    
    for p = perturbations % Loop through perturbations
        init_cond = [(1+p)*C_A_ss, (1+p)*T_ss, (1+p)*Tj_ss];
        tspan = [0 50];
        
        % Solve using ode45
        [t, y] = ode45(@(t, y) cstr_ode(t, y, F, V, k0, E, UA, Tj0, C_Af, Cp, rho, Fj, Vj, rho_j, Cp_j, H, T_f), tspan, init_cond);
        
        % Plot results
        plot(t, y(:,1), 'DisplayName', ['C_A, perturbed ', num2str(p*100), '%']);
        plot(t, y(:,2), 'DisplayName', ['T, perturbed ', num2str(p*100), '%']);
        plot(t, y(:,3), 'DisplayName', ['T_j, perturbed ', num2str(p*100), '%']);
    end
    legend;
    hold off;
end

% Function defining the ODE system
function dydt = cstr_ode(~, y, F, V, k0, E, UA, Tj0, C_Af, Cp, rho, Fj, Vj, rho_j, Cp_j, H, T_f)
    C_A = y(1);
    T = y(2);
    T_j = y(3);
    
    k = k0 * exp(-E / (1.987 * T)); % Arrhenius equation
    r = k * C_A; % Reaction rate
    
    dC_A_dt = (F/V) * (C_Af - C_A) - r;
    dT_dt = (F/V) * (T_f - T) + (-H/Cp) * r - (UA/(rho * Cp * V)) * (T - T_j);
    dTj_dt = (Fj/Vj) * (Tj0 - T_j) + (UA/(rho_j * Cp_j * Vj)) * (T - T_j);
    
    dydt = [dC_A_dt; dT_dt; dTj_dt];
end
