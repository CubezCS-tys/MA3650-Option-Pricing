% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % Example usage of the MC approach + compare
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % clear; clc;
% % 
% % r      = 0.05;
% % sigma  = 0.2;
% % T      = 1.0;
% % S0     = 100;
% % M0     = 80;    % Suppose the running min so far
% % Npaths = 5e5;   % number of simulations
% % Nsteps = 252;   % daily steps for 1 year
% % 
% % % We'll do a "shifted" approach to replicate the same scenario in MC:
% % % If the running min is M0 < S0 at time 0, we can simulate from M0 up to T
% % % but that is tricky to replicate exactly in a single-asset path from S0.
% % % For demonstration, let's just do standard MC from S0 and keep track of min.
% % 
% % % PDE solver call (still "explicit" version, but with more NT)
% % [Vgrid, Svec, Mvec] = LookbackFloatingCallPDE(r, sigma, T, Smax, NS, NM, NT);
% % 
% % % Suppose we want the PDE value at S0=100, M0=80
% % S0 = 100;
% % M0 = 80;
% % [~, iS0] = min(abs(Svec - S0));
% % [~, iM0] = min(abs(Mvec - M0));
% % pricePDE = Vgrid(iS0, iM0);
% % 
% % fprintf('PDE price at (S0=%.2f, M0=%.2f) = %.4f\n', Svec(iS0), Mvec(iM0), pricePDE);
% % 
% % % Now let's plot a slice for M0=80
% % valsPDE = Vgrid(:, iM0);  % This is dimension NS x 1
% % figure;
% % plot(Svec, valsPDE, 'LineWidth', 1.2);
% % xlabel('S'); ylabel('Option Value');
% % title('Lookback Floating-Strike Call (M=80 slice, t=0)');
% % grid on;
% % 
% % % -------------- Plot PDE vs Monte Carlo for varying S0 --------------
% % % For demonstration, fix M0=80, let S vary
% % iM0 = find(abs(Mvec - 80) < 1e-5, 1);
% % valsPDE = Vgrid(:, iM0);
% % 
% % figure;
% % plot(Svec, valsPDE, 'LineWidth', 1.2);
% % xlabel('S0'); ylabel('Value');
% % title('Floating-Strike Lookback Call via PDE (M=80, t=0)');
% % grid on;
% % 
% % % If you want, you can sample a few S0 points in Monte Carlo and overlay them.
% % hold on;
% % testS0s = [80, 90, 100, 110, 120];
% % mcVals = zeros(size(testS0s));
% % for k = 1:length(testS0s)
% %     mcVals(k) = LookbackFloatingCallMC(testS0s(k), r, sigma, T, 5e5, Nsteps);
% % end
% % plot(testS0s, mcVals, 'o');  % overlay MC points
% % legend('PDE','Monte Carlo Points','Location','Best');
% 
% 
% % -------------------------------------------------
% % Example driver code for the dimension-reduced PDE
% % -------------------------------------------------
% clear; clc;
% 
% % Parameters
% r     = 0.05;
% sigma = 0.2;
% T     = 1.0;
% Xmax  = 5.0;    % e.g. let X go up to 5
% NX    = 200;    % # X steps
% NT    = 5000;   % # time steps (explicit scheme => large # for stability)
% 
% [uGrid, Xvec] = LookbackFloatingCallPDE_1D(r, sigma, T, Xmax, NX, NT);
% 
% % The PDE solution at t=0 is in uGrid(:,1) [assuming 1-based indexing for times]
% u0 = uGrid(:,1);
% 
% % The option price for a given (S0, M0) is V(0,S0,M0)= M0 * u(0, X0),
% % where X0= S0 / M0.
% S0 = 100;
% M0 = 80;  % so X0=1.25
% X0 = S0/M0;
% % find index near X0
% [~, iX0] = min(abs(Xvec - X0));
% V0 = M0 * u0(iX0);
% 
% fprintf('Dimension-reduced PDE estimate at (S0=%.2f, M0=%.2f) ~ %.4f\n',...
%          S0, M0, V0);
% 
% % Plot the initial-time solution u(0,X) vs X
% figure;
% plot(Xvec, u0, 'LineWidth',1.2);
% xlabel('X = S / M'); ylabel('u(0, X) = V(0,S,M)/M');
% title('Dimension-Reduced PDE for Floating-Strike Lookback Call');
% grid on;


S0     = 100;
zmax   = 5;
K      = 80;
r      = 0.1;
T      = 1.0;
sigma  = 0.4;
Smax   = 400;   % max S used for discretization
ds     = 0.5;     % spacing in S
dt     = 0.01;  % time step

price = EuCallImpl3(S0, zmax, K, r, T, sigma, Smax, ds, dt);
fprintf('Computed price = %.4f\n', price);
