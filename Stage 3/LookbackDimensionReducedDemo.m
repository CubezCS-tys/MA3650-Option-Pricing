function LookbackDimensionReducedDemo()
% LOOKBACKDIMENSIONREDUCEDDEMO
%
% Demonstrates dimension reduction for a floating-strike lookback call
% option. We reduce the 2D PDE in (S, M) to a 1D PDE in (t, X) where X = S/M.
% Then we plot a 3D surface of u(t, X) = V(t,S,M)/M.

    % -----------------------------
    % 1) Model & Numerical Parameters
    % -----------------------------
    r     = 0.05;      % risk-free rate
    sigma = 0.2;       % volatility
    T     = 1.0;       % maturity (in years)
    Xmax  = 5.0;       % maximum X on the grid (X >= 1)
    NX    = 120;       % number of X steps
    NT    = 3000;      % number of time steps (explicit => large for stability)

    % -----------------------------
    % 2) Solve the PDE in 1D (t, X)
    % -----------------------------
    [uGrid, Xvec, tvec] = LookbackFloatingCallPDE_1D(r, sigma, T, Xmax, NX, NT);

    % -----------------------------
    % 3) 3D Plot of the solution
    %    We have uGrid(i,n) ~ u(tvec(n), Xvec(i)).
    % -----------------------------
    [TT, XX] = meshgrid(tvec, Xvec);  % Create a mesh for surf() dimension matching
    figure;
    surf(TT, XX, uGrid);        % dimension: [NX x (NT+1)] for each 
    shading interp;             % smoother shading
    xlabel('Time, t'); 
    ylabel('X = S / M');
    zlabel('u(t,X) = V(t,S,M)/M');
    title('Dimension-Reduced PDE Solution (Floating-Strike Lookback)');
    view(135, 30);              % adjust 3D viewing angle
    colorbar;                   % optional color scale

    % -----------------------------
    % 4) Example of extracting V(0,S0,M0)
    %    Suppose S0=100, M0=80 => X0=1.25
    % -----------------------------
    S0 = 100;
    M0 = 80;
    X0 = S0 / M0; 
    % find nearest X index
    [~, iX0] = min(abs(Xvec - X0));

    % PDE solution at t=0 is in column #1 of uGrid (time index 1 in a 1-based array)
    % but note that tvec goes from 0 to T => tvec(1) = 0, tvec(end)= T
    u0X = uGrid(iX0, 1);    % solution at (X0, t=0)

    V0 = M0 * u0X;  % dimension un-reduction: V(0,S0,M0)= M0 * u(0, X0)
    fprintf('Estimated option value at S0=%.2f, M0=%.2f: %.4f\n', S0, M0, V0);
end

% -------------------------------------------------------------------------
function [uGrid, Xvec, tvec] = LookbackFloatingCallPDE_1D(r, sigma, T, Xmax, NX, NT)
% LOOKBACKFLOATINGCALLPDE_1D
%
% Illustrative dimension-reduced PDE for a floating-strike lookback call.
% We define:
%   X = S / M  (>= 1),
%   u(t,X) = V(t,S,M) / M.
%
% Terminal condition at t=T : u(T,X)= (X - 1).
% Boundary conditions:
%   - at X=1 => u=0  (since S=M => V=0 => u=V/M=0)
%   - at X=Xmax => approximate u ~ X-1 (or a discount version). 
%
% PDE (example placeholder):
%   dU/dt = 0.5*sigma^2 * X^2 * d2U/dX^2 + r*X*dU/dX - r*U.
% In practice, the exact PDE for dimension-reduced lookback may differ slightly
% after deriving the SDE for X(t). Check references for the correct terms.
%
% Input:
%   r, sigma : market parameters
%   T        : maturity
%   Xmax     : maximum X
%   NX       : number of X steps
%   NT       : number of time steps
%
% Output:
%   uGrid : [NX x (NT+1)] array, where each column is u at a time level
%   Xvec  : grid in X (size NX)
%   tvec  : grid in time (size NT+1)

    % 1) X grid from 1 to Xmax
    Xvec = linspace(1, Xmax, NX);
    dX   = Xvec(2) - Xvec(1);

    % 2) Time grid from 0 to T
    tvec = linspace(0, T, NT+1);
    dt   = tvec(2) - tvec(1);   % uniform step

    % 3) Storage for the solution
    uGrid = zeros(NX, NT+1);

    % 4) Terminal condition at t=T => index = NT+1 in 1-based indexing
    %    so fill the last column
    for i = 1:NX
        uGrid(i, end) = Xvec(i) - 1;  % u(T,X)= X-1
    end

    % 5) Time-stepping backward: we want u(:,n) from u(:,n+1)
    for n = NT:-1:1
        % known solution at tvec(n+1)
        uOld = uGrid(:, n+1);
        % allocate new solution
        uNew = uOld; 

        for i = 2:NX-1
            Xval = Xvec(i);

            % PDE placeholder:
            % dU/dt = alpha * d2U/dX^2 + beta * dU/dX - r*U
            % with alpha= 0.5*sigma^2*X^2, beta= r*X.
            uip = uOld(i+1);
            uim = uOld(i-1);
            uic = uOld(i);

            d2udX2 = (uip - 2*uic + uim) / (dX^2);
            dudX   = (uip - uim) / (2*dX);

            alpha = 0.5 * sigma^2 * (Xval^2);
            beta  = r * Xval;

            dUdt  = alpha*d2udX2 + beta*dudX - r*uic;

            % Explicit Euler update
            uNew(i) = uic + dt * dUdt;
        end

        % Boundary conditions:
        uNew(1)   = 0;                   % at X=1 => u=0
        uNew(end) = Xvec(end) - 1;       % simple approximation at X=Xmax

        % Store
        uGrid(:, n) = uNew;
    end
end
