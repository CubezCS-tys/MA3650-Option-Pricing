% S0 = 50;
% Min = 30;
% r = 0.05;
% T = 1;
% sigma = 0.2;
% Smax = 500;
% dS = 0.5;
% dt = 0.001;

% S0 = 100;
% Min = 80;
% r = 0.1;
% T = 1;
% sigma = 0.4;
% Smax = 500;
% dS = 0.05;
% dt = 0.001;
% 
% price = FloatingStrikeCallImp(S0, Min, r, T, sigma, Smax, dS, dt)
% [callPrice, zGrid, matVal] = FloatingStrikeCallImp(S0, Min, r, T, sigma, Smax, dS, dt);
% figure; 
% plot(zGrid, matVal(:,1), 'b-','LineWidth',1.5);
% xlabel('z'); ylabel('u(z,0)'); title('Floating-Strike Lookback Call in z-space');
% grid on;
% 
% % Convert the z-grid back to the stock price grid: S = Min * z.
% SGrid = Min * zGrid;
% 
% % Compute the option values from the u(z,0) values:
% VGrid = Min * matVal(:,1);
% 
% % Plot the option value against the stock price
% figure;
% plot(SGrid, VGrid, 'b-', 'LineWidth', 1.5);
% xlabel('Stock Price S');
% ylabel('Option Value V(S, Min, 0)');
% title('Floating-Strike Lookback Call Option Value vs. Stock Price');
% grid on;
% 

S0 = 80;
Max = 200;
r = 0.1;
T = 1;
sigma = 0.4;
Smax = 500;
dS = 0.05;
dt = 0.001;

price = FloatingStrikePutImp(S0, Max, r, T, sigma, Smax, dS, dt)
[callPrice, zGrid, matVal] = FloatingStrikePutImp(S0, Max, r, T, sigma, Smax, dS, dt);
figure; 
plot(zGrid, matVal(:,1), 'b-','LineWidth',1.5);
xlabel('z'); ylabel('u(z,0)'); title('Floating-Strike Lookback Put in z-space');
grid on;

% Convert the z-grid back to the stock price grid: S = Min * z.
SGrid = Max * zGrid;

% Compute the option values from the u(z,0) values:
VGrid = Max * matVal(:,1);

% Plot the option value against the stock price
figure;
plot(SGrid, VGrid, 'b-', 'LineWidth', 1.5);
xlabel('Stock Price S');
ylabel('Option Value V(S, Max, 0)');
title('Floating-Strike Lookback Put Option Value vs. Stock Price');
grid on;




