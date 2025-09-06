function reactor_dynamics()
    % Given steady states
    steady_states = [
        0.59, 405, 390;
        1.12, 355, 340;
        2.51, 320, 310
    ];
    
    tspan = [0 50]; % Time span
    dt = 0.1; % Step size for RK4
    
    for i = 1:3
        CA0 = steady_states(i, 1);
        T0 = steady_states(i, 2);
        Tj0 = steady_states(i, 3);
        y0 = [CA0, T0, Tj0];
        
        % Solve using RK4
        [t_rk4, y_rk4] = rk4(@reactor_odes, tspan, y0, dt);
        
        % Solve using ode45
        [t_ode45, y_ode45] = ode45(@reactor_odes, tspan, y0);
        
        % Plot results
        figure;
        subplot(3,1,1);
        plot(t_rk4, y_rk4(:,1), 'r', t_ode45, y_ode45(:,1), 'b--');
        legend('RK4', 'ode45');
        title(['C_A vs Time for Steady State ', num2str(i)]);
        xlabel('Time'); ylabel('C_A');
        
        subplot(3,1,2);
        plot(t_rk4, y_rk4(:,2), 'r', t_ode45, y_ode45(:,2), 'b--');
        legend('RK4', 'ode45');
        title(['T vs Time for Steady State ', num2str(i)]);
        xlabel('Time'); ylabel('T');
        
        subplot(3,1,3);
        plot(t_rk4, y_rk4(:,3), 'r', t_ode45, y_ode45(:,3), 'b--');
        legend('RK4', 'ode45');
        title(['T_j vs Time for Steady State ', num2str(i)]);
        xlabel('Time'); ylabel('T_j');
    end
end

function dydt = reactor_odes(~, y)
    CA = y(1); T = y(2); Tj = y(3);
    % Define your ODEs here based on reactor model equations
    k = 0.1; % Example rate constant
    dydt = [-k * CA;
            -k * (T - 300);
            -k * (Tj - 300)];
end

function [t, y] = rk4(odefun, tspan, y0, dt)
    t = tspan(1):dt:tspan(2);
    y = zeros(length(t), length(y0));
    y(1,:) = y0;
    
    for i = 1:length(t)-1
        k1 = dt * odefun(t(i), y(i,:)');
        k2 = dt * odefun(t(i) + dt/2, y(i,:)' + k1/2);
        k3 = dt * odefun(t(i) + dt/2, y(i,:)' + k2/2);
        k4 = dt * odefun(t(i) + dt, y(i,:)' + k3);
        y(i+1,:) = y(i,:) + (k1' + 2*k2' + 2*k3' + k4')/6;
    end
end
