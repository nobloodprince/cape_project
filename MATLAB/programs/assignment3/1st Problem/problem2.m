clc, clearvars
% Function definition
f1 = @(y) y(2);
f2 = @(y) 1000*(1 - y(1)^2) * y(2) - y(1);

% Initial value
y0 = [2 0];

% Step size
h = 1; % sec

y = y0;
fprintf('Before start:\n');
disp(y)

for t = 0:300
    fprintf('%0.0f iteration\n', t+1)  % Fixed missing newline

    k11 = h * f1(y);
    k12 = h * f1(y + [k11/2, k11/2]); % Fixed argument passing
    k13 = h * f1(y + [k12/2, k12/2]);
    k14 = h * f1(y + [k13, k13]);

    k21 = h * f2(y);
    k22 = h * f2(y + [k11/2, k21/2]); % Ensuring correct argument update
    k23 = h * f2(y + [k12/2, k22/2]);
    k24 = h * f2(y + [k13, k23]);

    y(1) = y(1) + (k11 + 2*k12 + 2*k13 + k14) / 6;
    y(2) = y(2) + (k21 + 2*k22 + 2*k23 + k24) / 6;

    disp(y)
end
