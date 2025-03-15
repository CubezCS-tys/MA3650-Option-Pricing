% clc; clear; close all;
% 
% % Parameters
% S_values = linspace(80, 400, 200); % Range of spot prices for analytical pricing
% Min = 80;    % Running minimum
% M = 400;     % Maximum spot price over period (not used in PDE)
% r = 0.05;     % Risk-free rate (10%)
% sigma = 0.2; % Volatility (40%)
% T = 1;       % Time to expiry (1 year)
% 
% % Compute Analytical Lookback Call Price using a custom function "lookback_call"
% % (Ensure that you have a function "lookback_call" defined in your MATLAB path.)
% lookback_call_prices = arrayfun(@(S) lookback_call(S, Min, r, sigma, T), S_values);
% 
% % Compute Numerical Crank-Nicolson Solution using FloatingStrikeCallImp
% Smax = 500;   % Maximum stock price in the numerical grid
% dS = 0.05;    % Stock price step size
% dt = 0.001;   % Time step
% 
% % Call your numerical PDE solver (FloatingStrikeCallImp) to compute the solution.
% % It returns the option price, the z-grid, and the full solution matrix matVal.
% [callPrice, zGrid, matVal] = FloatingStrikeCallImp(100, Min, r, T, sigma, Smax, dS, dt);
% 
% % Plot 1: Plot the solution in z-space (u(z,0)) with green color.
% figure;
% plot(zGrid, matVal(:,1), 'g-', 'LineWidth', 1.5);
% xlabel('z');
% ylabel('u(z,0)');
% title('Floating-Strike Lookback Call in z-space');
% grid on;
% 
% % Convert the z-grid back to the stock price grid using S = Min * z.
% SGrid = Min * zGrid;
% 
% % Compute the option values from u(z,0) (recall V = m * u).
% VGrid = Min * matVal(:,1);
% 
% % Plot 2: Plot both Analytical and Numerical solutions versus the stock price.
% figure;
% hold on;
% % Plot Analytical solution with red solid line.
% plot(S_values, lookback_call_prices, 'r-', 'LineWidth', 1.5);
% % Plot Numerical solution with blue solid line.
% plot(SGrid, VGrid, 'b-', 'LineWidth', 1.5);
% xlabel('Stock Price (S)');
% ylabel('Lookback Option Value');
% title('Numerical vs Analytical Floating Strike Lookback Call');
% legend('Analytical Solution', 'Numerical Solution (Crank-Nicolson)', 'Location', 'NorthWest');
% grid on;
% hold off;

clc; clear;close all;

% Define Parameters

Min = 50;    % Running minimum (historical min price)
S_values = linspace(Min, 400, 200); % Range of spot prices for analytical pricing
r = 0.2;    % Risk-free rate
sigma = 0.1; % Volatility
T = 0.4;       % Time to expiry

% Compute Analytical Lookback Call Prices
lookback_call_prices = arrayfun(@(S) lookback_call(S, Min, r, sigma, T), S_values);

% Numerical Crank-Nicholson Parameters
Smax = 500;   % Maximum stock price in the numerical grid
dS = 0.5;    % Stock price step size
dt = 0.01;   % Time step

% Compute Numerical Solution for a range of S0 values
S0_values = linspace(Min, 400, 100); % Range of initial stock prices
callPrices = zeros(size(S0_values));

for i = 1:length(S0_values)
    S0 = S0_values(i);
    [callPrice, zGrid, matVal] = FloatingStrikeCallImp(S0, Min, r, T, sigma, Smax, dS, dt);
    callPrices(i) = callPrice;
end

% Convert z-space to Stock Price Grid
% SGrid = Min * zGrid;
% VGrid = Min * matVal(:,1); % Option value from numerical solution

% Plot 2: Compare Analytical vs Numerical Solution vs Stock Price
figure;
hold on;
% Analytical Solution (Red)
plot(S_values, lookback_call_prices, 'r-', 'LineWidth', 1.5);
% Numerical Solution (Blue)
%plot(SGrid, VGrid, 'b-', 'LineWidth', 1.5);
% Numerical Prices Over Range of S0 (Blue Dotted)
plot(S0_values, callPrices, 'b-', 'MarkerSize', 4);
xlabel('Stock Price S');
ylabel('Lookback Option Value');
title('Numerical vs Analytical Floating Strike Lookback Call');
legend('Analytical Solution', 'Numerical Solution (Crank-Nicholson)', 'Numerical Prices for Different S0', 'Location', 'NorthWest');
grid on;
hold off;

clc; clear; 

% Define Parameters
Min = 50;   % Running minimum (historical min price)
r = 0.05;   % Risk-free rate
sigma = 0.2; % Volatility
T = 1;      % Time to expiry

% Define stock price points for comparison
S_values = linspace(Min, 400, 10); % Stock prices from 50 to 400

% Compute Analytical Lookback Call Prices
analytical_prices = arrayfun(@(S) lookback_call(S, Min, r, sigma, T), S_values);

% Numerical Crank-Nicholson Parameters
Smax = 500;   % Maximum stock price in the numerical grid
dS = 0.5;     % Stock price step size
dt = 0.01;    % Time step

% Compute Numerical Solution for a range of S0 values
numerical_prices = zeros(size(S_values));

for i = 1:length(S_values)
    S0 = S_values(i);
    [callPrice, ~, ~] = FloatingStrikeCallImp(S0, Min, r, T, sigma, Smax, dS, dt);
    numerical_prices(i) = callPrice;
end

% Compute Absolute Error
absolute_error = abs(analytical_prices - numerical_prices);

% Create a Table
comparison_table = table(S_values', numerical_prices', analytical_prices', absolute_error', ...
    'VariableNames', {'Stock_Price', 'Numerical_Solution', 'Analytical_Solution', 'Absolute_Error'});

% Display Table
disp('Comparison of Numerical and Analytical Floating Strike Lookback Call Prices:');
disp(comparison_table);

% Save table to a CSV file
writetable(comparison_table, 'LookbackOption_Comparison.csv');



% clc; clear; close all;
% 
% % Parameters
% S0 = 100;  % Spot Price
% Min = 80;  % Minimum observed price
% Max = 200; % Maximum observed price
% r = 0.1;   % Risk-free rate
% T = 1;     % Time to expiry
% sigma = 0.4; % Volatility
% Smax = 500; % Maximum stock price in domain
% dS = 0.05;  % Stock price step size
% dt = 0.001; % Time step
% 
% % Compute Call and Put Prices
% [callPrice, putPrice, zGrid, matValCall, matValPut] = FloatingStrikeLookback(S0, Min, Max, r, T, sigma, Smax, dS, dt);
% 
% % Convert z-grid to Stock Price Grid
% SGrid = Min * zGrid;
% 
% % 3D Plot of the Call Option
% figure;
% surf(SGrid, linspace(0, T, size(matValCall, 2)), matValCall');
% xlabel('Stock Price S');
% ylabel('Time to Expiry T - t');
% zlabel('Call Option Value');
% title('3D Floating Strike Lookback Call Option');
% shading interp; colorbar; grid on; view(135, 30);
% 
% % 3D Plot of the Put Option
% figure;
% surf(SGrid, linspace(0, T, size(matValPut, 2)), matValPut');
% xlabel('Stock Price S');
% ylabel('Time to Expiry T - t');
% zlabel('Put Option Value');
% title('3D Floating Strike Lookback Put Option');
% shading interp; colorbar; grid on; view(135, 30);
% 
% % Print Final Prices
% fprintf('Floating Strike Lookback Call Price: %.4f\n', callPrice);
% fprintf('Floating Strike Lookback Put Price: %.4f\n', putPrice);

