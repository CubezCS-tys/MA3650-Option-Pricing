S0 = 100;
Min = 100;
r = 0.05;
T = 1;
sigma = 0.2;
Smax = 500;
dS = 0.5;
dt = 0.001;

price = FloatingStrikeCallImp(S0, Min, r, T, sigma, Smax, dS, dt)
    
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