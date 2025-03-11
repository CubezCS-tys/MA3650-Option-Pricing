function [price, zGrid, U] = FloatingStrikeCallCrankNicolson(S0, m, r, T, sigma, zMax, Nz, Nt)
%FLOATINGSTRIKECALLCRANKNICOLSON 
%  Price a floating-strike lookback call using Crank–Nicolson on the PDE
%
%    u_t + 0.5*sigma^2 * z^2 * u_{zz} + r*z*u_z - r*u = 0,
%
%  with terminal condition u(z,T) = max(z-1, 0) for 0 <= z <= zMax,
%  and boundary conditions:
%    u(0,t)   = 0,  (call payoff is zero if S=0)
%    u(zMax,t) ~ zMax - e^{-r (T - t)}  (as z-> large, payoff ~ z)
%
%  The ratio is z = S / m. Then the actual call price is C(S0,m,0) = m * u(S0/m, 0).
%
% INPUTS:
%   S0    - current underlying price
%   m     - 'minimum' value (the lookback reference), so z = S / m
%   r     - risk-free rate
%   T     - time to maturity
%   sigma - volatility
%   zMax  - maximum z on the truncated domain
%   Nz    - number of spatial steps
%   Nt    - number of time steps
%
% OUTPUTS:
%   price - the floating-strike call price at (S0,m,time=0)
%   zGrid - grid of z-values in [0, zMax]
%   U     - solution array U(i, n) = u(z_i, t_n),
%           where n=1 => t=T, n=Nt+1 => t=0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1) Set up grids
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dz    = zMax / Nz;
zGrid = linspace(0, zMax, Nz+1)';    % (Nz+1)-by-1 column vector: z=0..zMax

dt    = T / Nt;
tVec  = linspace(T, 0, Nt+1);       % from T down to 0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2) Initialize the solution array & set terminal condition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U = zeros(Nz+1, Nt+1);  
% At time T (index n=1),  u(z,T) = max(z - 1, 0)
U(:,1) = max(zGrid - 1, 0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3) Precompute PDE coefficients for the interior points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PDE: u_t + 0.5*sigma^2*z^2*u_{zz} + r*z*u_z - r*u = 0
% Discretize z in i=0..Nz.
% We'll build alpha_i, beta_i, gamma_i for i=1..Nz-1 as interior points.
alpha = zeros(Nz+1,1);
beta  = zeros(Nz+1,1);
gamma = zeros(Nz+1,1);

for i = 2:Nz
    zVal = zGrid(i);
    alpha(i) = 0.5 * sigma^2 * (zVal^2) / (dz^2) - 0.5 * r * zVal / dz;
    beta(i)  = - sigma^2 * (zVal^2) / (dz^2) - r;
    gamma(i) = 0.5 * sigma^2 * (zVal^2) / (dz^2) + 0.5 * r * zVal / dz;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4) Build the Crank-Nicolson matrices:
%%     (I + dt/2 * A) * U^{n+1} = (I - dt/2 * A) * U^{n} + BC adjustments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We'll build them as tri-diagonal. For i=2..Nz (1..(Nz-1) in "interior" sense).

% LHS = I + 0.5*dt*A  => subDiag, mainDiag, upDiag
subDiagL  = zeros(Nz+1,1);  
mainDiagL = ones(Nz+1,1);
upDiagL   = zeros(Nz+1,1);

% RHS = I - 0.5*dt*A
subDiagR  = zeros(Nz+1,1);  
mainDiagR = ones(Nz+1,1);
upDiagR   = zeros(Nz+1,1);

for i = 2:Nz
    mainDiagL(i) = 1 - 0.5*dt * beta(i);
    subDiagL(i)  =      -0.5*dt * alpha(i);
    upDiagL(i-1) =      -0.5*dt * gamma(i-1);  % be careful with indexing here

    mainDiagR(i) = 1 + 0.5*dt * beta(i);
    subDiagR(i)  =       0.5*dt * alpha(i);
    upDiagR(i-1) =       0.5*dt * gamma(i-1);
end

% But note that typically subDiag corresponds to i-1 row offset; let's
% rewrite more standardly:
%
% We'll do something simpler: set up the standard tri-di arrays for i=2..Nz-1
% then handle i=1, i=Nz with boundary conditions. Let's do it as a standard approach:

% Create vectors of length Nz-1 for interior only:
subL = zeros(Nz-1,1);
mainL= zeros(Nz-1,1);
upL  = zeros(Nz-1,1);
subR = zeros(Nz-1,1);
mainR= zeros(Nz-1,1);
upR  = zeros(Nz-1,1);

for i = 2 : Nz
    zVal = zGrid(i);
    if i<=Nz-1
        subL(i-1) = -0.5*dt * alpha(i);
        upL(i-1)  = -0.5*dt * gamma(i);
        subR(i-1) =  0.5*dt * alpha(i);
        upR(i-1)  =  0.5*dt * gamma(i);
    end
    mainL(i-1) = 1 - 0.5*dt * beta(i);
    mainR(i-1) = 1 + 0.5*dt * beta(i);
end

% Construct the two tridiagonal matrices for the interior (i=2..Nz):
LHS = spdiags([subL mainL upL], [-1 0 1], Nz-1, Nz-1);
RHS = spdiags([subR mainR upR], [-1 0 1], Nz-1, Nz-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5) Time stepping backwards: from n=1..Nt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n = 1:Nt
    t_n   = tVec(n);     % current time
    t_n1  = tVec(n+1);   % next time

    % "right-hand side" vector for interior:
    rhs = RHS * U(2:Nz, n);

    %----------------------
    % Apply boundary conditions
    %----------------------
    % i=1 => z=0:  u(0,t)=0 for a call
    U0 = 0;

    % i=Nz+1 => z=zMax: approximate large-z boundary
    %   For large z, a floating-strike call ~ z - e^{-r (T - t)}.
    %   So we set U(zMax,t_{n+1}) ~ zMax - e^{-r (T - t_{n+1})}.
    bc_zMax = zMax - exp(-r*(T - t_n1));

    % Incorporate these BC into the 'rhs':
    %    - subL(1) * U(1) because that would appear in row i=2 after expansion.
    % For row i=2 in interior system: subL(1)* U(1)
    rhs(1)   = rhs(1)   - subL(1)*U0;
    % For row i=Nz-1: upL(Nz-2)* U(Nz+1)
    rhs(end) = rhs(end) - upL(Nz-2)* bc_zMax;

    % Solve LHS * U(2:Nz, n+1) = rhs
    U(2:Nz, n+1) = LHS \ rhs;

    % Finally set the boundaries:
    U(1,   n+1) = U0;
    U(Nz+1,n+1) = bc_zMax;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 6) Extract the price at time 0 for S0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We marched from T down to 0, so U(:,end) is solution at t=0.
uAt0  = U(:, end);
z0    = S0/m;  % ratio
% Interpolate to get u(z0,0):
if z0 <= 0
    u0 = 0; 
elseif z0 >= zMax
    % if out of the grid, use the approx boundary:  z - e^{-rT}
    u0 = z0 - exp(-r*T); 
else
    u0 = interp1(zGrid, uAt0, z0, 'linear');
end

% Actual call price is  C(S0,m,0)= m * u0
price = m * u0;

end
