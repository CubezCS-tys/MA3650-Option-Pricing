% MATLAB Code for Pricing a Power Option and Plotting the Price vs. Underlying Price

% Clear workspace and command window


% Parameters for the option and underlying asset
K = 100;        % Strike price
r = 0.05;       % Risk-free interest rate (annual)
sigma = 0.2;    % Volatility (annual)
T = 1;          % Time to maturity in years
n = 2;          % Power exponent (payoff is max(S_T^n - K, 0))
numPaths = 1e5; % Number of Monte Carlo simulation paths

% ---------------------------
% Price the Power Option
% ---------------------------
% Assume a given initial stock price S0
S0 = 100;

% Simulate terminal stock prices under the risk-neutral measure:
% S_T = S0 * exp((r - 0.5*sigma^2)*T + sigma*sqrt(T)*Z)
Z = randn(numPaths, 1);  % Standard normal random numbers
ST = S0 * exp((r - 0.5*sigma^2)*T + sigma*sqrt(T)*Z);

% Compute the payoff of a power call option: max(S_T^n - K, 0)
payoffs = max(ST.^n - K, 0);

% Discount the average payoff to obtain the option price
price = exp(-r*T) * mean(payoffs);

fprintf('Monte Carlo Price for the Power Call Option with S0 = %.2f: %.4f\n', S0, price);

% ---------------------------
% Plot Option Price vs. Initial Stock Price S0
% ---------------------------
S0_range = linspace(10, 150, 100);  % Range of initial stock prices
prices = zeros(size(S0_range));     % Preallocate array for option prices

for i = 1:length(S0_range)
    S0_i = S0_range(i);
    % Simulate terminal prices for each S0_i
    ST_i = S0_i * exp((r - 0.5*sigma^2)*T + sigma*sqrt(T)*randn(numPaths, 1));
    % Compute payoff and option price for this S0_i
    payoff_i = max(ST_i.^n - K, 0);
    prices(i) = exp(-r*T) * mean(payoff_i);
end

% Create the plot
figure;
plot(S0_range, prices, 'LineWidth', 2);
xlabel('Initial Stock Price, S_0');
ylabel('Power Option Price');
title('Power Call Option Price as a Function of S_0');
grid on;
