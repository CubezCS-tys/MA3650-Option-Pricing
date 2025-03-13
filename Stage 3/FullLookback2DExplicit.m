function [price, S_grid, m_grid, V] = FullLookback2DExplicit(S0, r, q, sigma, T, Smax, dS, dm, dt)
% FullLookback2DExplicit prices a full floating-strike lookback call option
% using an explicit finite difference scheme on a 2D grid in (S, m).
%
% The full floating-strike lookback call has payoff:
%   V(T,S,m) = S - m,    for S >= m.
%
% Inputs:
%   S0    - initial asset price (and initial running minimum)
%   r     - risk-free interest rate
%   q     - dividend yield
%   sigma - volatility
%   T     - time to maturity
%   Smax  - maximum asset price on the grid
%   dS    - grid spacing in S
%   dm    - grid spacing in m
%   dt    - time step for the PDE solver
%
% Outputs:
%   price  - option price at t=0 for S = S0 and m = S0.
%   S_grid - vector of S grid points
%   m_grid - vector of m grid points
%   V      - solution matrix of option values at t=0 (rows: m, columns: S)

%% 1. Set up spatial grids
S_grid = 0:dS:Smax;        % S from 0 to Smax
m_grid = 0:dm:S0;          % m from 0 to S0 (running minimum cannot exceed S0 at inception)
Ns = length(S_grid);
Nm = length(m_grid);

% Initialize option value matrix V(t,S,m) at t = T (terminal condition)
V = nan(Nm, Ns);
for i = 1:Nm
    for j = 1:Ns
        if S_grid(j) >= m_grid(i)
            V(i,j) = S_grid(j) - m_grid(i);  % Terminal payoff
        end
    end
end

%% 2. Time-stepping: solve PDE backward from t = T to t = 0
Nt = round(T/dt);
dt = T / Nt;  % adjust dt accordingly

for n = 1:Nt
    V_new = V;  % copy current solution
    
    % Loop over interior grid points (in S and m)
    for i = 1:Nm
        for j = 2:Ns-1  % for S (avoid boundaries at S=0 and S=Smax)
            if S_grid(j) >= m_grid(i)
                % Second derivative in S (central difference)
                V_SS = (V(i, j+1) - 2*V(i, j) + V(i, j-1)) / dS^2;
                % First derivative in S (central difference)
                V_S = (V(i, j+1) - V(i, j-1)) / (2*dS);
                % Explicit Euler time stepping for the Black-Scholes operator:
                V_new(i,j) = V(i,j) - dt * (0.5 * sigma^2 * S_grid(j)^2 * V_SS + r * S_grid(j) * V_S - r * V(i,j));
            end
        end
    end
    
    % Boundary conditions:
    % (a) At S = Smax: assume V is given by extrapolation (copy the neighboring value)
    for i = 1:Nm
        V_new(i, end) = V_new(i, end-1);
    end
    
    % (b) Along the lower boundary S = m (for each m row): impose Neumann condition in m.
    % Approximate ∂V/∂m = 0 by copying the value from the next row.
    for i = 1:Nm-1
        % Find index j where S is closest to m_grid(i)
        [~, j_index] = min(abs(S_grid - m_grid(i)));
        if j_index >= 1 && j_index <= Ns
            V_new(i, j_index) = V_new(i+1, j_index);
        end
    end
    
    % (c) At m = 0 (first row): impose Neumann condition (copy from second row)
    V_new(1,:) = V_new(2,:);
    
    % Update solution for the next time step:
    V = V_new;
end

%% 3. Extract option price at t = 0 for S = S0 and m = S0
[~, j0] = min(abs(S_grid - S0));
[~, i0] = min(abs(m_grid - S0));
price = V(i0, j0);

end
