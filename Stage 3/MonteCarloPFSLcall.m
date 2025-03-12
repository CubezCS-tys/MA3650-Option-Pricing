%% Monte Carlo Pricing for a Partial Floating-Strike Lookback Call Option
% This script prices an option with payoff:
%    payoff = max(S(T) - fixedMin, 0)
% where fixedMin is the locked-in minimum value (as in your PDE).
%
% This is different from the full lookback option where the running minimum 
% is computed along each simulated path.

clear; clc; close all;

%% Parameters
S0       = 100;      % Initial asset price
fixedMin = 80;       % Fixed (locked-in) minimum, same as "Min" in your PDE
r        = 0.1;      % Risk-free interest rate
sigma    = 0.4;      % Volatility
T        = 1;        % Time to maturity (years)
numPaths = 100000;   % Number of simulated paths
numSteps = 1000;      % Number of time steps (e.g., trading days in a year)
dt       = T/numSteps;

%% Simulate Asset Paths
S = zeros(numSteps+1, numPaths);
S(1,:) = S0;
dW = sqrt(dt)*randn(numSteps, numPaths);

for i = 2:numSteps+1
    S(i,:) = S(i-1,:) .* exp((r - 0.5*sigma^2)*dt + sigma*dW(i-1,:));
end

%% Compute Payoffs Using Fixed Minimum
% Here we do NOT compute the running minimum dynamically.
% Instead, the payoff is calculated using the given fixed minimum.
S_T = S(end, :);
payoff = max(S_T - fixedMin, 0);

%% Discount and Estimate Price
discountFactor = exp(-r*T);
priceMC = discountFactor * mean(payoff);
stdError  = discountFactor * std(payoff) / sqrt(numPaths);

fprintf('Monte Carlo Price (Partial Lookback): %.4f\n', priceMC);
fprintf('Standard Error: %.4f\n', stdError);

%% (Optional) Plot Sample Paths and Payoff Histogram
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
