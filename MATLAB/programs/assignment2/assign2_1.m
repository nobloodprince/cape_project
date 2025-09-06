 %contants
 clc
 F = 1;
 V = 1;
 ko = 36*10e6;
 delta_H = 6500;
 E = 12000;
 p_Cp = 500;
 Tf = 298;
 CAf = 10;
 UA = 150;
 Tjo = 298;
 pj_Cj = 600;
 Fj = 1.25;
 Vj = 0.25;
 R = 8.314;

%rate of reaction
%r = ko*(exp(-E/(R*T)))*CA;

tic;

%func1 =  F*CAf - F*CA - (ko*(exp(-E/(R*T)))*CA)*V;
%func2 = p*Cp*F*(Tf-T)+(delta_H)*V*r - UA*(T-Tj);
%func3 = pj_Cj*Fj*(Tjo-Tj)+UA*(T-Tj);

syms CA T Tj
% Define functions symbolically
F_sym = [ F*CAf - F*CA - (ko*(exp(-E/(R*T)))*CA)*V;
          p_Cp*F*(Tf-T)+(delta_H)*V*(ko*(exp(-E/(R*T)))*CA) - UA*(T-Tj);
          pj_Cj*Fj*(Tjo-Tj)+UA*(T-Tj)];

% Compute Jacobian matrix
J_sym = jacobian(F_sym, [CA, T, Tj]);

% Convert symbolic functions to numerical functions
F_func = matlabFunction(F_sym, 'Vars', {CA, T, Tj});
J_func = matlabFunction(J_sym, 'Vars', {CA, T, Tj});

% Initial guess

X = [5; 350; 800];  
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
        %fprintf('Converged in %d iterations\n', k);
        break;
    end
end

t1 = toc;

% Display results
fprintf('Newton Rapshon Method\n');
fprintf('Solution: CA = %.6f, T = %.6f, Tj = %.6f\n', X(1), X(2), X(3));
fprintf('time taken(sec): %.6f\n', t1);
fprintf('Converged in %d iterations\n', k);
fprintf('\n');


