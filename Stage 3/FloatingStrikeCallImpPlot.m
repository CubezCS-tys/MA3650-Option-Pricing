% S0 = 50;
% Min = 50;
% r = 0.05;
% T = 1;
% sigma = 0.2;
% Smax = 500;
% dS = 0.5;
% dt = 0.001;
% 
% 
% price = FloatingStrikeCallImp(S0, Min, r, T, sigma, Smax, dS, dt)
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

% Example usage:

S0    = 100;     % current spot
M0    = 80;      % running minimum
r     = 0.05;    % risk-free rate
q     = 0.00;    % dividend yield
sigma = 0.20;    % volatility
T     = 1.0;     % time to maturity (in years)

Nx    = 400;     % number of x steps
Xmax  = 3.0;     % maximum x-value (S/M up to 3x or so)
Nt    = 400;     % number of time steps

price_est = LookbackFloatingCall_CN(S0, M0, r, q, sigma, T, Nx, Xmax, Nt);
fprintf('Estimated floating-strike call price: %.4f\n', price_est);
