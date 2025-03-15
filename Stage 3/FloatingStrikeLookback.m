function [uGrid, Xvec, tvec] = FloatingStrikeLookbackCall_CN(r, sigma, T, Xmax, NX, NT)
% FloatingStrikeLookbackCall_CN solves the 1D PDE for the dimension-reduced
% floating-strike lookback call option.
%
% The dimension reduction is achieved by setting X = S/m and u(t,X) = V(t,S,m)/m.
% The PDE in u(t,X) is:
%
%   u_t + 0.5*sigma^2*X^2*u_{XX} + r*X*u_X - r*u = 0,  for X > 1, 0 <= t < T,
%
% with terminal condition:
%   u(T,X) = X - 1,
%
% and boundary conditions:
%   u(t,1) = 0, and u(t,Xmax) = Xmax - exp(-r*(T-t)).
%
% INPUTS:
%   r     - risk-free rate
%   sigma - volatility
%   T     - time to maturity
%   Xmax  - maximum value of X (domain: X in [1, Xmax])
%   NX    - number of spatial steps (grid in X will have NX+1 points)
%   NT    - number of time steps (time grid has NT+1 levels)
%
% OUTPUTS:
%   uGrid - (NX+1 x NT+1) matrix, with uGrid(i,n) ~ u(t_n, X_i)
%   Xvec  - spatial grid vector in X (from 1 to Xmax)
%   tvec  - time grid vector (from 0 to T)

% Spatial grid
Xvec = linspace(1, Xmax, NX+1)';  % column vector from 1 to Xmax
dX = (Xmax - 1) / NX;

% Time grid
tvec = linspace(0, T, NT+1);  % from 0 to T
dt = T / NT;

% Preallocate solution matrix: each column is a time level.
uGrid = zeros(NX+1, NT+1);

% Terminal condition: u(T,X) = X - 1.
uGrid(:, end) = Xvec - 1;

% We solve for interior nodes i=2 to NX (there are NX-1 interior nodes)
N_int = NX - 1;

% Precompute coefficients at interior nodes (for i = 2,..., NX)
% Define arrays A, B, C of length (NX-1) corresponding to nodes i=2,...,NX.
A = zeros(N_int,1); % coefficient multiplying u_{i-1}
B = zeros(N_int,1); % coefficient multiplying u_i
C = zeros(N_int,1); % coefficient multiplying u_{i+1}
for i = 2:NX
    X_i = Xvec(i);
    % Finite difference approximations:
    % Second derivative: central difference ~ (u_{i+1} - 2*u_i + u_{i-1})/dX^2.
    % First derivative: central difference ~ (u_{i+1} - u_{i-1})/(2*dX).
    % Multiply these by the PDE coefficients.
    A(i-1) = 0.5 * sigma^2 * X_i^2 / dX^2 - r * X_i / (2*dX);
    B(i-1) = - sigma^2 * X_i^2 / dX^2 - r;
    C(i-1) = 0.5 * sigma^2 * X_i^2 / dX^2 + r * X_i / (2*dX);
end

% For Crank-Nicolson, we need to set up two matrices.
% The interior system has size N_int x N_int.
% We define:
%   lower (sub-diagonal): - (dt/2) * A(2:end)  [size: N_int-1]
%   main (diagonal):         1 - (dt/2) * B        [size: N_int]
%   upper (super-diagonal): - (dt/2) * C(1:end-1)   [size: N_int-1]
lower = - (dt/2) * A(2:end);      % from i = 3 to NX
main  = 1 - (dt/2) * B;           % for i = 2,...,NX
upper = - (dt/2) * C(1:end-1);    % from i = 2 to NX-1

LHS = diag(main) + diag(upper,1) + diag(lower,-1);

% Right-hand side operator (explicit part)
% lower_RHS: (dt/2) * A(2:end)
% main_RHS: 1 + (dt/2) * B
% upper_RHS: (dt/2) * C(1:end-1)
lower_RHS = (dt/2) * A(2:end);
main_RHS  = 1 + (dt/2) * B;
upper_RHS = (dt/2) * C(1:end-1);

RHS_mat = diag(main_RHS) + diag(lower_RHS, -1) + diag(upper_RHS, 1);

% LU factorization of LHS (remains constant for all time steps)
[L_mat, U_mat] = lu(LHS);

% Backward time-stepping: from time index NT to 1
for n = NT:-1:1
    % Boundary conditions at time t = tvec(n):
    % Left boundary: at X = 1, u(t,1)=0.
    uGrid(1, n) = 0;
    % Right boundary: at X = Xmax, use asymptotic behavior:
    uGrid(end, n) = Xmax - exp(-r*(T - tvec(n)));
    
    % Build the RHS vector for interior nodes i=2,...,NX:
    rhs = RHS_mat * uGrid(2:NX, n+1);
    
    % Adjust for the boundary contributions:
    % For the first interior node (i=2): contribution from left boundary u(1)
    rhs(1) = rhs(1) - lower(1)*uGrid(1, n);
    % For the last interior node (i=NX): contribution from right boundary u(NX+1)
    rhs(end) = rhs(end) - upper(end)*uGrid(end, n);
    
    % Solve the linear system for interior nodes:
    u_interior = U_mat \ (L_mat \ rhs);
    uGrid(2:NX, n) = u_interior;
end

end
