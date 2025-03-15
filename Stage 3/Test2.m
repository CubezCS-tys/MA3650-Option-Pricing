clc; clear;close all;

% Define Parameters

Max = 120;    % Running minimum (historical min price)
S_values = linspace(0, Max, 200); % Range of spot prices for analytical pricing
r = 0.05;    % Risk-free rate
sigma = 0.2; % Volatility
T = 1;       % Time to expiry

% Compute Analytical Lookback Call Prices
lookback_put_prices = arrayfun(@(S) lookback_put(S, Max, r, sigma, T), S_values);

% Numerical Crank-Nicholson Parameters
Smax = 500;   % Maximum stock price in the numerical grid
dS = 0.5;    % Stock price step size
dt = 0.01;   % Time step

% Compute Numerical Solution for a range of S0 values
S0_values = linspace(0, Max, 100); % Range of initial stock prices
putPrices = zeros(size(S0_values));

for i = 1:length(S0_values)
    S0 = S0_values(i);
    [putPrice, zGrid, matVal] = FloatingStrikePutImp(S0, Max, r, T, sigma, Smax, dS, dt);
    putPrices(i) = putPrice;
end

% Convert z-space to Stock Price Grid
% SGrid = Min * zGrid;
% VGrid = Min * matVal(:,1); % Option value from numerical solution

% Plot 1: Solution in z-space (u(z,0))
figure;
plot(zGrid, matVal(:,1), 'g-', 'LineWidth', 1.5);
xlabel('z');
ylabel('u(z,0)');
title('Floating-Strike Lookback Put in z-space');
grid on;

% Plot 2: Compare Analytical vs Numerical Solution vs Stock Price
figure;
hold on;
% Analytical Solution (Red)
plot(S_values, lookback_put_prices, 'r-', 'LineWidth', 1.5);
% Numerical Solution (Blue)
%plot(SGrid, VGrid, 'b-', 'LineWidth', 1.5);
% Numerical Prices Over Range of S0 (Blue Dotted)
plot(S0_values, putPrices, 'b-', 'MarkerSize', 4);
xlabel('Stock Price S');
ylabel('Lookback Option Value');
title('Numerical vs Analytical Floating Strike Lookback Call');
legend('Analytical Solution', 'Numerical Solution (Crank-Nicholson)', 'Numerical Prices for Different S0', 'Location', 'NorthWest');
grid on;
hold off;

clc; clear; 

% Define Parameters
Max = 120;   % Running minimum (historical min price)
r = 0.05;   % Risk-free rate
sigma = 0.2; % Volatility
T = 1;      % Time to expiry

% Define stock price points for comparison
S_values = linspace(0, Max, 10); % Stock prices from 50 to 400

% Compute Analytical Lookback Call Prices
analytical_prices = arrayfun(@(S) lookback_put(S, Max, r, sigma, T), S_values);

% Numerical Crank-Nicholson Parameters
Smax = 500;   % Maximum stock price in the numerical grid
dS = 0.5;     % Stock price step size
dt = 0.01;    % Time step

% Compute Numerical Solution for a range of S0 values
numerical_prices = zeros(size(S_values));

for i = 1:length(S_values)
    S0 = S_values(i);
    [putPrice, ~, ~] = FloatingStrikePutImp(S0, Max, r, T, sigma, Smax, dS, dt);
    numerical_prices(i) = putPrice;
end

% Compute Absolute Error
absolute_error = abs(analytical_prices - numerical_prices);

% Create a Table
comparison_table = table(S_values', numerical_prices', analytical_prices', absolute_error', ...
    'VariableNames', {'Stock_Price', 'Numerical_Solution', 'Analytical_Solution', 'Absolute_Error'});

% Display Table
disp('Comparison of Numerical and Analytical Floating Strike Lookback Put Prices:');
disp(comparison_table);

% Save table to a CSV file
writetable(comparison_table, 'LookbackOption_Comparison.csv');
