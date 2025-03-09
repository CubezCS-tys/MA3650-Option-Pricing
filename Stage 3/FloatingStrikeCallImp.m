% function [price, vetS, matval] = FloatingStrikeCallImp(S0,Min,r,T,sigma,Smax,dS,dt)
% % set up grid and adjust increments if necessary
% Sz = S0/Min;
% z = Smax/Min;
% M = round((z)/dS);
% dS = (z)/M;
% N = round(T/dt);
% dt = T/N;
% matval = zeros(M+1,N+1);
% vetz = linspace(0,z,M+1)';
% veti = 0:M;
% vetj = 0:N;
% % set up boundary conditions
% matval(:,N+1) = vetz;
% matval(1,:) = 0;
% matval(M+1,:) = max(z-1, 0);
% % set up the tridiagonal coefficients matrix
% a = 0.5*(r*dt*veti-sigma^2*dt*(veti.^2));
% b = 1+sigma^2*dt*(veti.^2)+r*dt;
% c = -0.5*(r*dt*veti+sigma^2*dt*(veti.^2));
% coeff = diag(a(3:M),-1) + diag(b(2:M)) + diag(c(2:M-1),1);
% [L,U] = lu(coeff);
% % solve the sequence of linear systems
% aux = zeros(M-1,1);
% for j=N:-1:1
%    aux(M-1) = - c(M) * matval(M+1,j);
%    matval(2:M,j) = U \ (L \ (matval(2:M,j+1) + aux));
% end
% % return price, possibly by linear interpolation outside the grid
% price = Min*interp1(vetz, matval(:,1), Sz);
% matval;
% 

function [price, vetS, matval] = FloatingStrikeCallImp(S0, Min, r, T, sigma, Smax, dS, dt)
% FloatingStrikeCallImp
% Prices a floating-strike lookback call option using an implicit finite-difference scheme.
%
% INPUTS:
%   S0    - Current underlying price
%   Min   - Running minimum of S (used to transform to z = S/Min)
%   r     - Risk-free interest rate
%   T     - Time to maturity
%   sigma - Volatility
%   Smax  - Maximum underlying price for the grid (used to set z_max)
%   dS    - Step size in z (initial guess)
%   dt    - Time step (initial guess)
%
% OUTPUTS:
%   price - Option price computed at time 0
%   vetS  - Grid in the original S-space (S = Min * z)
%   matval- PDE solution matrix in z-space

    % Set up grid and adjust increments if necessary
    Sz = S0 / Min;      % dimensionless S0
    z   = Smax / Min;   % maximum z
    M   = round(z / dS);
    dS  = z / M;
    N   = round(T / dt);
    dt  = T / N;
    matval = zeros(M+1, N+1);
    vetz = linspace(0, z, M+1)';  % z grid from 0 to z
    % veti = 0:M;   % Not used further
    % vetj = 0:N;   % Not used further

    % Set up boundary conditions
    % Terminal condition: payoff u(z,T) = max(z - 1, 0)
    matval(:, N+1) = max(vetz - 1, 0);
    % Lower boundary: at z = 0, u(0,t) = 0
    matval(1, :) = 0;
    % Upper boundary: at z = z_max, set to max(z - 1, 0)
    matval(M+1, :) = max(z - 1, 0);

    % Set up the tridiagonal coefficients matrix using z values
    % Note: using vetz instead of the index vector (veti)
    a = 0.5 * (r*dt*vetz - sigma^2*dt*(vetz.^2));
    b = 1 + sigma^2*dt*(vetz.^2) + r*dt;
    c = -0.5 * (r*dt*vetz + sigma^2*dt*(vetz.^2));

    % Construct matrix for interior nodes: indices 2:M (exclude boundaries)
    coeff = diag(a(3:M), -1) + diag(b(2:M)) + diag(c(2:M-1), 1);
    [L, U] = lu(coeff);

    % Solve the sequence of linear systems (backward time-stepping)
    for j = N:-1:1
        aux = zeros(M-1, 1);  % Reset aux for each time step
        % Right boundary contribution from node i = M+1 affecting interior node i = M:
        aux(M-1) = - c(M) * matval(M+1, j);
        % Left boundary (i = 1) is zero, so no contribution is needed.
        matval(2:M, j) = U \ (L \ (matval(2:M, j+1) + aux));
    end

    % Return price by interpolating in the z grid at z = Sz and scaling back to S-space
    price = Min * interp1(vetz, matval(:, 1), Sz, 'linear', 'extrap');
    % Build the S grid: S = Min * z
    vetS = Min * vetz;
end
