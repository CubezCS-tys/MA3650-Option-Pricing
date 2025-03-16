%% Compare PDE and Monte Carlo Prices for a Partial Floating-Strike Lookback Call
% This script computes the option price for different initial stock prices S0
% using your PDE function and a Monte Carlo simulation (both pricing the partial
% lookback, i.e. payoff = max(S(T)-fixedMin,0)).

clear; clc; close all;

%% Parameters common to both methods
r = 0.05;
sigma = 0.2;
T = 1;
fixedMin = 80;      % locked-in minimum (same as "Min" for the PDE function)
Smax = 200;         % for the PDE grid
dS = 0.5;           % spatial grid spacing for PDE (in scaled space)
dtPDE = 0.01;      % time step for PDE
numPaths = 1000;  % Monte Carlo: number of simulated paths
numSteps = 10;    % Monte Carlo: number of time steps

%% Range of initial stock prices S0 (must be >= fixedMin)
S0_values = linspace(fixedMin, 200, 50);
PDE_prices = zeros(size(S0_values));
MC_prices = zeros(size(S0_values));

%% Loop over S0 values and compute prices
for i = 1:length(S0_values)
    S0 = S0_values(i);
    % PDE price from your function FloatingStrikeCallImp:
    % Syntax: [price, vetz, matval] = FloatingStrikeCallImp(S0, Min, r, T, sigma, Smax, dS, dt)
    [pricePDE, ~, ~] = FloatingStrikeCallImp(S0, fixedMin, r, T, sigma, Smax, dS, dtPDE);
    PDE_prices(i) = pricePDE;
    
    % Monte Carlo price using the helper function below:
    priceMC = MonteCarloPartialLookback(S0, fixedMin, r, sigma, T, numPaths, numSteps);
    MC_prices(i) = priceMC;

    lookback_call_prices = arrayfun(@(S) lookback_call(S, fixedMin, r, sigma, T), S0_values);

end

%% Plot the results
figure;
plot(S0_values, PDE_prices, 'b-o', 'LineWidth', 2);
hold on;
plot(S0_values, MC_prices, 'r-s', 'LineWidth', 2);
hold on;
plot(S0_values, lookback_call_prices, 'g-s', 'LineWidth', 2);
xlabel('Initial Stock Price, S_0');
ylabel('Option Price');
title('Partial Floating-Strike Lookback Call: PDE vs Monte Carlo');
legend('PDE','Monte Carlo','Location','NorthWest');
grid on;
