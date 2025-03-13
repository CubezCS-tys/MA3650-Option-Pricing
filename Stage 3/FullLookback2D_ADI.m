function [price, S_grid, m_grid, V] = FullLookback2D_ADI(S0, r, sigma, T, Smax, dS, dm, dt)
% FullLookback2D_ADI prices a full floating-strike lookback call option
% using an ADI scheme on a 2D grid in (S, m).
%
% The option value V(t,S,m) satisfies:
%   V_t + 0.5*sigma^2 S^2 V_SS + rS V_S - rV = 0,   for S>= m,
% with terminal condition at t=T:  V(T,S,m) = S - m.
% The boundary condition along S=m is approximated by a Neumann condition:
%   dV/dm = 0  when S = m.
%
% Inputs:
%   S0    - initial asset price (and initial running minimum)
%   r     - risk-free rate
%   sigma - volatility
%   T     - time to maturity
%   Smax  - maximum asset price on grid
%   dS    - grid spacing in S
%   dm    - grid spacing in m (m goes from 0 to S0)
%   dt    - time step for the ADI scheme
%
% Outputs:
%   price  - option price at t = 0 for S = S0 and m = S0.
%   S_grid - vector of S grid points
%   m_grid - vector of m grid points
%   V      - solution matrix at t = 0 (rows correspond to m, columns to S)

%% 1. Set up spatial grids
S_grid = 0:dS:Smax;            % S in [0, Smax]
m_grid = 0:dm:S0;              % m in [0, S0]
Ns = length(S_grid);
Nm = length(m_grid);

% Initialize option value matrix at terminal time T using the payoff V(T,S,m)=S-m.
V = nan(Nm, Ns);
for i = 1:Nm
    for j = 1:Ns
        if S_grid(j) >= m_grid(i)
            V(i,j) = S_grid(j) - m_grid(i);
        end
    end
end

%% 2. Set up coefficients for spatial derivatives in S direction
% For central differences:
%  V_S  ~ (V(j+1) - V(j-1))/(2*dS)
%  V_SS ~ (V(j+1) - 2*V(j) + V(j-1))/(dS^2)
%
% We use these coefficients in the implicit parts of the ADI scheme.
% (For simplicity, we assume the PDE operator is applied only in S.)
A = @(S) 0.5 * sigma^2 * S.^2;  % coefficient for V_SS
B = @(S) r * S;               % coefficient for V_S

%% 3. Time stepping using ADI (backward in time)
Nt = round(T/dt);
dt = T/Nt;  % adjust dt

for n = 1:Nt
    % --- ADI Step 1: Implicit in S, explicit in m ---
    V_half = V;  % initialize intermediate solution
    % For each fixed m (each row), solve implicit system in S.
    for i = 1:Nm
        % Build tridiagonal system for S at row i
        % For j = 2:Ns-1, the discretization:
        %   (I - theta*dt*L_S) V_half(i,j) = RHS, with theta = 0.5.
        % Here L_S V = A(S_grid(j)) * (V(j+1)-2V(j)+V(j-1))/(dS^2) + B(S_grid(j))*(V(j+1)-V(j-1))/(2*dS) - r*V(j)
        theta = 0.5;
        a = zeros(Ns-2,1); % lower diagonal
        b_diag = zeros(Ns-2,1); % main diagonal
        c = zeros(Ns-2,1); % upper diagonal
        RHS = zeros(Ns-2,1);
        
        for j = 2:Ns-1
            S_val = S_grid(j);
            A_val = A(S_val);
            B_val = B(S_val);
            
            a(j-1) = - theta*dt*( A_val/(dS^2) - B_val/(2*dS) );
            b_diag(j-1) = 1 + theta*dt*(2*A_val/(dS^2) + r);
            c(j-1) = - theta*dt*( A_val/(dS^2) + B_val/(2*dS) );
            
            % RHS from explicit part:
            RHS(j-1) = V(i,j) + (1-theta)*dt*( ...
                A_val*(V(i,j+1) - 2*V(i,j) + V(i,j-1))/(dS^2) + ...
                B_val*(V(i,j+1) - V(i,j-1))/(2*dS) - r*V(i,j) );
        end
        
        % Adjust for boundary conditions in S:
        % At j = 1 (S=0): V is known (often 0)
        % At j = Ns (S=Smax): use extrapolated value (we copy the neighbor)
        % Here we simply modify RHS:
        RHS(1) = RHS(1) - a(1)*V(i,1);
        RHS(end) = RHS(end) - c(end)*V(i,end);
        
        % Solve tridiagonal system for row i:
        V_row = thomas_algorithm(a, b_diag, c, RHS);
        % Update V_half(i,2:Ns-1)
        V_half(i,2:Ns-1) = V_row';
    end
    
    % --- ADI Step 2: Implicit in m, explicit in S ---
    V_new = V_half;
    % For each fixed S (each column), solve implicit system in m.
    for j = 1:Ns
        % We assume that in the m direction, the PDE has no second derivative
        % term (the PDE is independent of m in the interior) so the coupling in m
        % enters only through the boundary conditions.
        % Here we use a simple implicit step (or simply average with adjacent m-values)
        % to enforce a Neumann condition along m.
        
        % For interior m points i = 2:Nm-1, apply:
        % V_new(i,j) = 0.5*(V_half(i+1,j) + V_half(i-1,j));
        for i = 2:Nm-1
            V_new(i,j) = 0.5*(V_half(i+1,j) + V_half(i-1,j));
        end
        
        % For the boundaries in m:
        % At m = 0 (i=1) we set V_new(1,j) = V_new(2,j) (Neumann)
        V_new(1,j) = V_new(2,j);
        % At m = S0 (i=Nm), we simply copy the interior value:
        V_new(Nm,j) = V_new(Nm-1,j);
    end
    
    % Update solution for next time step:
    V = V_new;
end

%% 4. Extract Option Price at t=0 for S=S0 and m=S0
[~, j0] = min(abs(S_grid - S0));
[~, i0] = min(abs(m_grid - S0));
price = V(i0, j0);
end

%% Thomas Algorithm for tridiagonal system (helper function)
function x = thomas_algorithm(a, b, c, d)
% Solves a tridiagonal system Ax=d, where a, b, c are the sub-, diag-, and
% super-diagonals respectively.
    n = length(d);
    c_star = zeros(n,1);
    d_star = zeros(n,1);
    
    % Forward sweep:
    c_star(1) = c(1)/b(1);
    d_star(1) = d(1)/b(1);
    for i = 2:n
        denom = b(i) - a(i)*c_star(i-1);
        c_star(i) = c(i)/denom;
        d_star(i) = (d(i) - a(i)*d_star(i-1))/denom;
    end
    
    % Back substitution:
    x = zeros(n,1);
    x(n) = d_star(n);
    for i = n-1:-1:1
        x(i) = d_star(i) - c_star(i)*x(i+1);
    end
end
