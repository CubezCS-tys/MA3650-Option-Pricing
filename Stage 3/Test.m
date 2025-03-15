clc; clear; close all;

% Parameters
S_values = linspace(80, 200, 200); % Range of spot prices
Min = 80;    % Minimum spot price over period
M = 400;     % Maximum spot price over period
r = 0.1;     % Risk-free rate (10%)
sigma = 0.4; % Volatility (40%)
T = 1;       % Time to expiry (1 year)

% Compute Analytical Lookback Call Price
lookback_call_prices = arrayfun(@(S) lookback_call(S, Min, r, sigma, T), S_values);

% Compute Numerical Crank-Nicholson Solution
Smax = 200;   % Maximum stock price in the numerical grid
dS = 0.05;       % Stock price step size
dt = 0.001;    % Time step


[callPrice, zGrid, matVal] = FloatingStrikeCallImp(100, Min, r, T, sigma, Smax, dS, dt);
figure; 
plot(zGrid, matVal(:,1), 'b-','LineWidth',1.5);
xlabel('z'); ylabel('u(z,0)'); title('Floating-Strike Lookback Call in z-space');
grid on;

% Convert the z-grid back to the stock price grid: S = Min * z.
SGrid = Min * zGrid;

% Compute the option values from the u(z,0) values:
VGrid = Min * matVal(:,1);

% Plot the option value against the stock price
figure;
plot(SGrid, VGrid, 'b-', 'LineWidth', 1.5);
% Plot Both Solutions
plot(S_values, lookback_call_prices, 'b', 'LineWidth', 1.5); hold on;
plot(SGrid, VGrid, 'b-', 'LineWidth', 1.5);
% plot(vetz * Min, matval(:,1), 'r--', 'LineWidth', 1.5);
xlabel('Stock Price (S)');
ylabel('Lookback Option Value');
title('Numerical vs Analytical Floating Strike Lookback Call');
legend('Analytical Solution', 'Numerical Solution (Crank-Nicholson)', 'Location', 'NorthWest');
grid on;
