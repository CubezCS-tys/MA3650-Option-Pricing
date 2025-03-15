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

% We solve for interior nodes i=2 to NX (since i=1 and i=NX+1 are set by boundaries)
N_int = NX - 1;

% Precompute coefficients at interior nodes (they are time-independent here)
A = zeros(N_int,1); % coefficient for u_{i-1}
B = zeros(N_int,1); % coefficient for u_i
C = zeros(N_int,1); % coefficient for u_{i+1}
for i = 2:NX
    X_i = Xvec(i);
    % Second derivative term coefficient: 0.5*sigma^2*X_i^2/dX^2
    % First derivative term coefficient: r*X_i/(2*dX)
    A(i-1) = 0.5 * sigma^2 * X_i^2 / dX^2 - r * X_i / (2*dX);
    B(i-1) = - sigma^2 * X_i^2 / dX^2 - r;
    C(i-1) = 0.5 * sigma^2 * X_i^2 / dX^2 + r * X_i / (2*dX);
end

% For Crank-Nicolson, define the following matrices:
% The implicit (left-hand) operator is: I - dt/2 * L, and the explicit (right-hand)
% operator is: I + dt/2 * L.
% For interior node i (i=2,...,NX), L_i u = A(i-1)*u_{i-1} + B(i-1)*u_i + C(i-1)*u_{i+1}.
%
% So for each interior node we have:
%   (1 - dt/2*B(i-1)) u_i^n - dt/2*A(i-1)*u_{i-1}^n - dt/2*C(i-1)*u_{i+1}^n
% = (1 + dt/2*B(i-1)) u_i^{n+1} + dt/2*A(i-1)*u_{i-1}^{n+1} + dt/2*C(i-1)*u_{i+1}^{n+1}
%
% Assemble the tridiagonal matrix for the left-hand side (LHS):
lower = - (dt/2)*A;         % sub-diagonal (size: N_int-1)
main  = 1 - (dt/2)*B;         % main diagonal (size: N_int)
upper = - (dt/2)*C;         % super-diagonal (size: N_int-1)
LHS = diag(main) + diag(upper,1) + diag(lower,-1);

% For the right-hand side (RHS) operator, we need:
lower_RHS = (dt/2)*A;
main_RHS  = 1 + (dt/2)*B;
upper_RHS = (dt/2)*C;
RHS_mat = diag(main_RHS) + diag(lower_RHS, -1) + diag(upper_RHS, 1);

% LU factorization of LHS (remains constant for all time steps)
[L_mat, U_mat] = lu(LHS);

% Backward time-stepping
for n = NT:-1:1  % time index n goes from NT to 1 (backward)
    % Boundary conditions at time t = tvec(n):
    % Left boundary: at X=1, u(t,1)=0.
    uGrid(1, n) = 0;
    % Right boundary: at X=Xmax, use asymptotic behavior:
    uGrid(end, n) = Xmax - exp(-r*(T - tvec(n)));
    
    % Build the RHS vector for interior nodes i=2,...,NX:
    rhs = RHS_mat * uGrid(2:NX, n+1);
    
    % Adjust for the boundary values:
    % For the first interior node (i=2): contribution from u(1)
    rhs(1) = rhs(1) - lower(1)*uGrid(1, n);
    % For the last interior node (i=NX): contribution from u(NX+1)
    rhs(end) = rhs(end) - upper(end)*uGrid(end, n);
    
    % Solve for the interior nodes at time step n:
    u_interior = U_mat \ (L_mat \ rhs);
    uGrid(2:NX, n) = u_interior;
end

end
