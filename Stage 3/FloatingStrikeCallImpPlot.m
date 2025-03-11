% S0 = 50;
% Min = 30;
% r = 0.05;
% T = 1;
% sigma = 0.2;
% Smax = 500;
% dS = 0.5;
% dt = 0.001;

S0 = 70;
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


%     
% S0     = 100;    % Current underlying price
% m      = 100;    % Running minimum (constant)
% r      = 0.05;   % 5% interest
% T      = 1.0;    % 1 year
% sigma  = 0.2;    % 20% volatility
% Smax   = 300;    % Max underlying for the grid
% dz     = 0.01;   % Step in z
% dt     = 0.001;   % Time step
% [callPrice, zGrid, matVal] = FloatingStrikeCallImp(S0, m, r, T, sigma, Smax, dz, dt);
% disp(callPrice)
% 
% figure; 
% plot(zGrid, matVal(:,1), 'b-','LineWidth',1.5);
% xlabel('z'); ylabel('u(z,0)'); title('Floating-Strike Lookback Call in z-space');
% grid on;



