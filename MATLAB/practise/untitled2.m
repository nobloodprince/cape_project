% Given values
T = 250 + 273; % Convert to Kelvin
P = 10; % Pressure
Tc = 407.5; % Critical Temperature
Pc = 111.3; % Critical Pressure
R = 0.08206; % Gas constant

% Calculate constants
a = (27 * R^2 * Tc^2) / (64 * Pc);
b = R * Tc / (8 * Pc);

% Define volume range, ensuring v > b
v_min = b + 0.01; % Slightly above b to avoid singularity
v_max = 20;
v = linspace(v_min, v_max, 1000); % More points for smoothness

% Define function
func = (P + (a ./ (v.^2))) ./ (v - b) - R * T;

% Plot the function
figure;
plot(v, func, 'b-', 'LineWidth', 2);
hold on;
yline(0, 'k--'); % Add horizontal line at y = 0 for reference
xlabel('Volume (v)');
ylabel('Function Value');
title('Variation of Function with Volume');
grid on;
