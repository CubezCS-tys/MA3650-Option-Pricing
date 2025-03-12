function [price, S_grid, m_grid, V] = FloatingStrikeLookbackFull(S0, r, T, sigma, Smax, dS, dm, dt)
% FullLookback2DExplicit prices a full floating-strike lookback call option
% using an explicit finite difference scheme on a 2D grid in (S, m).
%
% The full floating-strike lookback call has payoff:
%    V(T,S,m) = S - m,    for S >= m.
%
% The PDE in the region S>=m is:
%    V_t + 0.5*sigma^2*S^2*V_{SS} + r*S*V_S - r*V = 0.
%
% We define a grid:
%    S in [0, Smax] with spacing dS,
%    m in [0, S0] with spacing dm,
% and enforce V(t,S,m)=NaN (unused) when S < m.
%
% Boundary conditions:
%   - Terminal: V(T,S,m) = S - m,  for S>=m.
%   - At S = Smax: V(t,Smax,m) = Smax - m.
%   - Along S = m: impose a Neumann condition in m (approximate by
%     setting V(t,m,m) equal to the value just above in m).
%
% Inputs:
%   S0   - initial asset price (and initial minimum)
%   r    - risk-free interest rate
%   T    - time to maturity
%   sigma- volatility
%   Smax - maximum asset price to consider
%   dS   - asset price grid spacing
%   dm   - minimum grid spacing (for running minimum m)
%   dt   - time step (should be very small for stability)
%
% Outputs:
%   price - the option price at t = 0 for S=S0 and m=S0.
%   S_grid - vector of S grid points
%   m_grid - vector of m grid points
%   V      - matrix of option values at t=0 (rows correspond to m, cols to S)

%% Set up the grids
S_grid = 0:dS:Smax;  % S from 0 to Smax
m_grid = 0:dm:S0;    % m from 0 to S0 (initial minimum is S0)
Ns = length(S_grid);
Nm = length(m_grid);

% Initialize V at terminal time T
V = nan(Nm, Ns);  % We'll fill only valid points (S>= m)

for i = 1:Nm
    for j = 1:Ns
        if S_grid(j) >= m_grid(i)
            V(i,j) = S_grid(j) - m_grid(i);
        end
    end
end

% Number of time steps (backward in time)
Nt = round(T/dt);

%% Time stepping: explicit Euler method (backward in time)
% Loop backward from t = T to t = 0.
for n = 1:Nt
    V_new = V;  % Copy current solution
    % Loop over m-index and S-index for interior points.
    % We update only if S >= m.
    for i = 1:Nm
        for j = 2:Ns-1  % use central differences in S; j=1 and j=Ns are boundaries.
            if S_grid(j) >= m_grid(i)
                % Second derivative in S:
                V_SS = (V(i, j+1) - 2*V(i, j) + V(i, j-1)) / dS^2;
                % First derivative in S:
                V_S = (V(i, j+1) - V(i, j-1)) / (2*dS);
                % Explicit time step:
                V_t = - (0.5*sigma^2 * S_grid(j)^2 * V_SS + r*S_grid(j)*V_S - r*V(i,j));
                V_new(i,j) = V(i,j) + dt * V_t;
            end
        end
    end
    
    % Enforce boundary conditions:
    % At S = Smax: Dirichlet: V = Smax - m.
    for i = 1:Nm
        V_new(i, end) = Smax - m_grid(i);
    end
    
    % Along S = m (the lower boundary of the valid region):
    % For each m value, find the grid point in S closest to m.
    % Here we approximate the Neumann condition by setting V(i, j) equal to
    % the value at the next m level.
    for i = 1:Nm-1
        % Find the index j such that S_grid(j) is approximately equal to m_grid(i)
        % (since S_grid starts at 0, j=round(m/dS)+1 is a rough approximation)
        j = round(m_grid(i)/dS) + 1;
        if j < 1 || j > Ns, continue; end
        % Impose V_m ~ 0: we set V(i,j) = V(i+1,j)
        V_new(i,j) = V_new(i+1,j);
    end
    
    % Optionally, at m = 0 (first row) we can impose a Neumann condition:
    V_new(1,:) = V_new(2,:);
    
    % Update the solution for the next time step:
    V = V_new;
end

%% Interpolate the price at the initial condition: S = S0, m = S0.
% Find the indices where S is S0 and m is S0.
[~, j0] = min(abs(S_grid - S0));
[~, i0] = min(abs(m_grid - S0));
price = V(i0, j0);

end
