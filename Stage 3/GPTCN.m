clc; clear; close all;

% Parameters
Smax = 200;   % Max stock price
Smin = 50;    % Min stock price in domain
M = 1000;      % Number of stock price steps
N = 10000;     % Number of time steps
T = 1;        % Time to maturity
r = 0.1;      % Risk-free rate
sigma = 0.4;  % Volatility
dS = (Smax - Smin) / M;  % Stock price step
dt = T / N;   % Time step

% Stock price grid
S = linspace(Smin, Smax, M+1)'; 

% Terminal Condition: Payoff at Expiry
V = max(S - Smin, 0); 

% Coefficients for Crank-Nicholson
alpha = (0.25 * dt * (sigma^2 * (S.^2) / dS^2 - (r * S) / dS));
beta = -0.5 * dt * (sigma^2 * (S.^2) / dS^2 + r);
gamma = (0.25 * dt * (sigma^2 * (S.^2) / dS^2 + (r * S) / dS));

% Tridiagonal matrices
A = diag(1 - beta(2:M)) + diag(-alpha(3:M), -1) + diag(-gamma(2:M-1), 1);
B = diag(1 + beta(2:M)) + diag(alpha(3:M), -1) + diag(gamma(2:M-1), 1);

% Time stepping with Crank-Nicholson
for n = N:-1:1
    rhs = B * V(2:M);
    V(2:M) = A \ rhs; % Solve the tridiagonal system
end

% Plot results
figure;
plot(S, V, 'b', 'LineWidth', 1.5);
xlabel('Stock Price (S)');
ylabel('Lookback Option Value');
title('Crank-Nicholson Pricing of Floating Strike Lookback Call');
grid on;
legend('Numerical Solution');

