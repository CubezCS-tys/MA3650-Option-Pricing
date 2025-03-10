%% Floating-Strike Lookback Call Analytical Solution
% This script computes the analytical price of a floating-strike lookback call
% using the following formulas:
%
%   a_1 = [log(S0/Smin) + (r - q + 0.5*sigma^2)*T] / (sigma*sqrt(T))
%   a_2 = a_1 - sigma*sqrt(T)
%   a_3 = [log(S0/Smin) + (-r + q + 0.5*sigma^2)*T] / (sigma*sqrt(T))
%   Y_1 = - (2*(r - q - 0.5*sigma^2)*log(S0/Smin)) / (sigma^2)
%
% The analytical call price is given by:
%
%   CFL = S0*exp(-q*T)*normcdf(a_1)
%         - S0*exp(-q*T)*((sigma^2)/(2*(r-q)))*normcdf(-a_1)
%         - Smin*exp(-r*T)*[normcdf(a_2) - ((sigma^2)/(2*(r-q)))*exp(Y_1)*normcdf(-a_3)]
%
% The script then plots CFL as a function of S0 over a specified range.

clear; close all; clc;

%% Parameter Setup
% Fixed parameters
Smin = 80;      % Running minimum
r    = 0.1;     % Risk-free interest rate
q    = 0;       % Dividend yield
sigma = 0.4;    % Volatility
T    = 1;       % Time to maturity (years)

% Choose a range for the current stock price S0.
S0_min = 0;
S0_max = 450;
numS = 100;
S0_range = linspace(S0_min, S0_max, numS);

% Preallocate analytical call price
CFL_values = zeros(size(S0_range));

%% Compute the analytical solution for each S0
for i = 1:length(S0_range)
    S0 = S0_range(i);
    CFL_values(i) = CFL(S0, Smin, r, q, sigma, T);
end

%% Plot the Analytical Solution
figure;
plot(S0_range, CFL_values, 'b-', 'LineWidth', 1.5);
xlabel('Current Stock Price, S_0');
ylabel('Floating-Strike Lookback Call Price');
title('Analytical Floating-Strike Lookback Call Price');
grid on;

%% Local Functions

function a1 = a_1(S0, Smin, r, q, sigma, T)
    % Compute a_1 as defined in the analytical formula
    a1 = (log(S0/Smin) + (r - q + 0.5*sigma^2)*T) / (sigma*sqrt(T));
end

function a2 = a_2(S0, Smin, r, q, sigma, T)
    % Compute a_2 using a_1
    a1 = a_1(S0, Smin, r, q, sigma, T);
    a2 = a1 - sigma*sqrt(T);
end

function a3 = a_3(S0, Smin, r, q, sigma, T)
    % Compute a_3 directly
    a3 = (log(S0/Smin) + (-r + q + 0.5*sigma^2)*T) / (sigma*sqrt(T));
end

function Y1 = Y_1(S0, Smin, r, q, sigma, T)
    % Compute Y_1 as defined in the formula
    Y1 = -(2*(r - q - 0.5*sigma^2)*log(S0/Smin))/(sigma^2);
end

function cfl = CFL(S0, Smin, r, q, sigma, T)
    % Compute the floating-strike lookback call price analytically.
    % Note: This formula assumes r ~= q.
    a1 = a_1(S0, Smin, r, q, sigma, T);
    a2 = a_2(S0, Smin, r, q, sigma, T);
    a3 = a_3(S0, Smin, r, q, sigma, T);
    Y1 = Y_1(S0, Smin, r, q, sigma, T);
    
    term1 = S0 * exp(-q*T) * normcdf(a1);
    term2 = S0 * exp(-q*T) * ((sigma^2)/(2*(r-q))) * normcdf(-a1);
    term3 = Smin * exp(-r*T) * ( normcdf(a2) - ((sigma^2)/(2*(r-q)))*exp(Y1)*normcdf(-a3) );
    
    cfl = term1 - term2 - term3;
end
