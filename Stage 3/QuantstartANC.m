clc; clear; close all;

% Parameters
S_values = linspace(1, 200, 200); % Range of spot prices
m = 80;    % Minimum spot price over period
M = 200;   % Maximum spot price over period
r = 0.1;   % Risk-free rate (10%)
v = 0.4;   % Volatility (40%)
T = 1;     % Time to expiry (1 year)

% Compute Lookback Call and Put prices
lookback_call_prices = arrayfun(@(S) lookback_call(S, m, r, v, T), S_values);
lookback_put_prices = arrayfun(@(S) lookback_put(S, M, r, v, T), S_values);

% Plotting
figure;
plot(S_values, lookback_call_prices, 'b', 'LineWidth', 1.5); hold on;
plot(S_values, lookback_put_prices, 'r', 'LineWidth', 1.5);
xlabel('Stock Price (S)');
ylabel('Lookback Option Value');
title('Lookback Call and Put Option Prices vs Spot Price');
legend('Lookback Call', 'Lookback Put', 'Location', 'NorthWest');
grid on;

% Functions

% Approximation to the cumulative normal distribution function
function cdf = norm_cdf(x)
    k = 1.0 / (1.0 + 0.2316419 * abs(x));
    k_sum = k * (0.319381530 + k * (-0.356563782 + k * (1.781477937 + k * (-1.821255978 + 1.330274429 * k))));
    cdf = 1 - (1 / sqrt(2 * pi)) * exp(-0.5 * x^2) * k_sum;
    
    if x < 0
        cdf = 1 - cdf;
    end
end

% Functions a_1, a_2, a_3
function result = a_1(S, H, r, v, T)
    result = (log(S / H) + (r + 0.5 * v^2) * T) / (v * sqrt(T));
end

function result = a_2(S, H, r, v, T)
    result = a_1(S, H, r, v, T) - v * sqrt(T);
end

function result = a_3(S, H, r, v, T)
    result = a_1(S, H, r, v, T) - (2 * r * sqrt(T) / v);
end

% Lookback Call Option Pricing
function price = lookback_call(S, m, r, v, T)
    a1 = a_1(S, m, r, v, T);
    a2 = a_2(S, m, r, v, T);
    a3 = a_3(S, m, r, v, T);
    
    term1 = S * norm_cdf(a1);
    term2 = m * exp(-r * T) * norm_cdf(a2);
    mult = (S * v^2) / (2 * r);
    term3 = norm_cdf(-a1) - exp(-r * T) * (m / S)^((2 * r) / v^2) * norm_cdf(-a3);
    
    price = term1 - term2 - mult * term3;
end

% Lookback Put Option Pricing
function price = lookback_put(S, M, r, v, T)
    a1 = a_1(S, M, r, v, T);
    a2 = a_2(S, M, r, v, T);
    a3 = a_3(S, M, r, v, T);
    
    term1 = -S * norm_cdf(-a1);
    term2 = M * exp(-r * T) * norm_cdf(-a2);
    mult = (S * v^2) / (2 * r);
    term3 = norm_cdf(a1) - exp(-r * T) * (M / S)^((2 * r) / v^2) * norm_cdf(a3);
    
    price = term1 + term2 + mult * term3;
end
