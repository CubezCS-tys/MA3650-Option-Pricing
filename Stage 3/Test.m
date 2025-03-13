% demo_FullLookback2D_ADI.m
clear; clc; close all;

% Parameters
S0 = 100;      % initial asset price (and running minimum at t=0)
r = 0.1;
sigma = 0.4;
T = 1;         % time to maturity
Smax = 200;    % maximum S in grid
dS = 0.5;        % grid spacing in S
dm = 0.5;        % grid spacing in m (m in [0,S0])
dt = 0.001;    % time step for ADI

[price, S_grid, m_grid, V] = FullLookback2D_ADI(S0, r, sigma, T, Smax, dS, dm, dt);
fprintf('Full floating-strike lookback call price at t=0 (S=S0, m=S0): %.4f\n', price);

% Suppose we have V, S_grid, m_grid from the PDE solver at t=0.
% We want the slice where the running minimum m = S0.

S0 = 100;  % the same S0 you used in your PDE code

% 1) Find the index in m_grid closest to S0
[~, i0] = min(abs(m_grid - S0));

% 2) Extract that row from V
vals_at_m_S0 = V(i0, :);

% 3) Optionally, only plot where S >= S0 (because for S < S0 it's invalid)
validMask = (S_grid >= S0);

% 4) Plot
figure;
plot(S_grid(validMask), vals_at_m_S0(validMask), 'b-o', 'LineWidth', 2);
xlabel('Asset Price S');
ylabel('Option Value');
title('Lookback Option Value at t=0, for m = S_0');
grid on;
