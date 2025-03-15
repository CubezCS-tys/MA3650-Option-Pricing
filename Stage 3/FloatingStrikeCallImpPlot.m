% S0 = 50;
% Min = 30;
% r = 0.05;
% T = 1;
% sigma = 0.2;
% Smax = 500;
% dS = 0.5;
% dt = 0.001;
% 
S0 = 100;
Min = 80;
r = 0.1;
T = 1;
sigma = 0.4;
Smax = 500;
dS = 0.05;
dt = 0.001;

price = FloatingStrikeCallImp(S0, Min, r, T, sigma, Smax, dS, dt)
[callPrice, zGrid, matVal] = FloatingStrikeCallImp(S0, Min, r, T, sigma, Smax, dS, dt);
figure; 
plot(zGrid, matVal(:,1), 'b-','LineWidth',1.5);
xlabel('z'); ylabel('u(z,0)'); title('Floating-Strike Lookback Call in z-space');
grid on;

% Convert the z-grid back to the stock price grid: S = Min * z.
SGrid = Min * zGrid;

% Compute the option values from the u(z,0) values:
VGrid = Min * matVal(:,1);

% Plot the option value against the stock price
figure;
plot(SGrid, VGrid, 'b-', 'LineWidth', 1.5);
xlabel('Stock Price S');
ylabel('Option Value V(S, Min, 0)');
title('Floating-Strike Lookback Call Option Value vs. Stock Price');
grid on;


% clc; clear; close all;
% 
% % Define parameters
% S0    = 100;
% Min   = 80;
% r     = 0.05;
% T     = 1;
% sigma = 0.2;
% Smax  = 500;
% dS    = 0.005;
% dt    = 0.0001;
% 
% % Compute the Floating Strike Lookback Call
% [callPrice, zGrid, matVal] = FloatingStrikeCallImp(S0, Min, r, T, sigma, Smax, dS, dt);
% 
% % Convert the z-grid back to stock price grid: S = Min * z
% SGrid = Min * zGrid;
% VGrid = Min * matVal;
% 
% % Time grid (discretized from 0 to T)
% timeGrid = linspace(0, T, size(matVal, 2));
% 
% % ---------------------------
% % CHOOSE YOUR STARTING S HERE
% % ---------------------------
% Sstart = 80;  % or 80, 50, etc.
% i0 = find(SGrid >= Sstart, 1, 'first');  % index where S >= Sstart
% 
% % Slice the data
% SGridPlot   = SGrid(i0:end);
% matValPlot  = VGrid(i0:end, :);
% 
% % 3D Surface Plot
% figure;
% surf(SGridPlot, timeGrid, matValPlot');  % Transpose so time is on Y-axis
% xlabel('Stock Price S');
% ylabel('Time to Expiry T - t');
% zlabel('Option Value V(S, Min, t)');
% title('Floating-Strike Lookback Call (S \ge Starting Value)');
% colorbar;
% shading interp;
% grid on;
% view(135, 30);
% 
% 
