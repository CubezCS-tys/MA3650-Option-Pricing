clc; clear; close all;

% Parameters
S0 = 150;   % Initial stock price
K = 140;     % Strike price
r = 0.05;   % Risk-free rate
T = 0.5;      % Time to maturity
sigma = 0.2; % Volatility
Smax = 2 * max(S0, K) * exp(r * T); % Maximum asset price for CN method
dS = 1;     % Space step (increased for efficiency)
dt = 0.001;  % Time step (increased for efficiency)

% p values from 1 to 2 with step 0.1
p_values = 1:0.1:2;
num_p = length(p_values);
results = zeros(num_p, 4); % Store [p, BS Price, CN Price, Error]

for i = 1:num_p
    p = p_values(i);

    % Black-Scholes Pricing
    bs_price = BlackScholesPowerCall(S0, K, T, 0, r, sigma, p);

    % Crank-Nicholson Pricing
    cn_price = CrankNicholsonCall(S0, K, r, T, sigma, Smax, dS, dt, p);

    % Compute absolute error
    error = abs(bs_price - cn_price);

    % Store results
    results(i, :) = [p, bs_price, cn_price, error];
end

% Save results to CSV file
%writematrix('power_option_pricing_errors.csv', results);

% Display results
disp(array2table(results, 'VariableNames', {'p', 'BlackScholesPrice', 'CrankNicholsonPrice', 'AbsoluteError'}));




