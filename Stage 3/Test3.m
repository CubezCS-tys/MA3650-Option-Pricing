% Example parameters
r = 0.05;
sigma = 0.2;
T = 5.0;
Xmax = 1.0;      % large enough ratio
Nx = 500;        % spatial steps
Nt = 1000;        % time steps

[u0, Xgrid, U] = LookbackCallFloating_CN(r, sigma, T, Xmax, Nx, Nt);
% 'u0' is u(0,1).

% If your underlying starts at S0 and M(0)=S0 => ratio X=1,
% the floating‐strike call value is:
S0 = 100;  
callPrice = S0 * u0;

disp(['Floating‐strike call = ', num2str(callPrice)]);
