%% Monte Carlo Pricing of a Floating-Strike Lookback Call Option
% This script simulates asset paths under the Black–Scholes model,
% computes the running minimum for each path, and then prices the option.
% The payoff at maturity is: S(T) - min_{0<=t<=T} S(t)
% The option price is the discounted average payoff.

clear; clc; close all;

%% Parameters
S0       = 30;     % initial asset price
r        = 0.1;    % risk-free interest rate
sigma    = 0.4;     % volatility
T        = 1;       % time to maturity (years)
numPaths = 100000;  % number of simulated paths
numSteps = 1000;     % number of time steps (e.g., trading days in a year)
dt       = T/numSteps;

%% Simulate Asset Paths
% Pre-allocate the matrix for asset prices:
S = zeros(numSteps+1, numPaths);
S(1,:) = S0;

% Generate random increments (each column is one path)
dW = sqrt(dt) * randn(numSteps, numPaths);

% Simulate the asset paths using the exact solution for geometric Brownian motion:
% S(t+dt) = S(t) * exp((r - 0.5*sigma^2)*dt + sigma*sqrt(dt)*Z)
for i = 2:numSteps+1
    S(i,:) = S(i-1,:) .* exp((r - 0.5*sigma^2)*dt + sigma*dW(i-1,:));
end

%% Compute Payoffs
% For each path, determine the minimum asset price reached:
minS = min(S, [], 1);
% The final asset price for each path:
S_T = S(end, :);
% Floating-strike lookback call payoff:
payoff = S_T - minS;  % always nonnegative

%% Discount and Estimate Price
discountFactor = exp(-r*T);
priceMC = discountFactor * mean(payoff);
stdError  = discountFactor * std(payoff) / sqrt(numPaths);

fprintf('Monte Carlo Price: %.4f\n', priceMC);
fprintf('Standard Error: %.4f\n', stdError);

%% (Optional) Plot a few sample paths and histogram of payoffs
figure;
plot(0:dt:T, S(:,1:50));
xlabel('Time (years)');
ylabel('Asset Price');
title('Sample Simulated Asset Paths');

figure;
histogram(payoff,50);
xlabel('Payoff');
ylabel('Frequency');
title('Histogram of Payoffs');
