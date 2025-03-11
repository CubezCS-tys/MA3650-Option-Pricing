%% Floating-Strike Lookback Put Option: Price vs. Stock Price
clear; clc; close all;

% Parameters
r     = 0.1;    % risk-free rate
q     = 0;      % dividend yield
sigma = 0.4;    % volatility
T     = 1;   % time to maturity (years)
Smax  = 400;     % maximum observed stock price

% Define a range for S0 (current stock price) where S0 <= Smax.
S0_vec = linspace(10, Smax, 100);  % for example, from 30 to 50
PFL_vec = zeros(size(S0_vec));

% Compute option price for each S0
for i = 1:length(S0_vec)
    S0 = S0_vec(i);
    PFL_vec(i) = PFL(S0, Smax, r, q, sigma, T);
end

% Plot the option value vs. the stock price
figure;
plot(S0_vec, PFL_vec, 'LineWidth', 2);
xlabel('Stock Price, S_0');
ylabel('Floating-Strike Lookback Put Price');
title('Floating-Strike Lookback Put Option Price vs. Stock Price');
grid on;

%% Local Function Definitions

function pfl = PFL(S0, Smax, r, q, sigma, T)
    % Computes the floating-strike lookback put price using the given parameters.
    b1 = b_1(S0, Smax, r, q, sigma, T);
    b2 = b_2(S0, Smax, r, q, sigma, T);
    b3 = b_3(S0, Smax, r, q, sigma, T);
    Y2 = Y_2(S0, Smax, r, q, sigma, T);
    
    term1 = (Smax * exp(-r * T)) * ( normcdf(b1) - ((sigma^2) / (2*(r-q))) * exp(Y2) * normcdf(-b3) );
    term2 = (S0 * exp(-q * T)) * ((sigma^2) / (2*(r-q))) * normcdf(-b2) - (S0 * exp(-q * T)) * normcdf(b2);
    
    pfl = term1 + term2;
end

function b1 = b_1(S0, Smax, r, q, sigma, T)
    % Computes b1 based on the parameters.
    b1 = ( log(Smax/S0) + (-r + q + 0.5 * sigma^2) * T ) / ( sigma * sqrt(T) );
end

function b2 = b_2(S0, Smax, r, q, sigma, T)
    % Computes b2 based on b1.
    b1_val = b_1(S0, Smax, r, q, sigma, T);
    b2 = b1_val - sigma * sqrt(T);
end

function b3 = b_3(S0, Smax, r, q, sigma, T)
    % Computes b3 based on the parameters.
    b3 = ( log(Smax/S0) + (r - q - 0.5 * sigma^2) * T ) / ( sigma * sqrt(T) );
end

function Y2 = Y_2(S0, Smax, r, q, sigma, T)
    % Computes Y2 based on the parameters.
    Y2 = (2 * (r - q - 0.5 * sigma^2) * log(Smax/S0)) / (sigma^2);
end
