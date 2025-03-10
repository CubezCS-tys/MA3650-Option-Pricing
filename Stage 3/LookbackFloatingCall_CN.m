function LookbackFloatingCall_CN()
% LOOKBACKFLOATINGCALL_CN
% ----------------------------------------------------------
% This script derives and solves (via Crank–Nicolson) the PDE
% for a floating-strike lookback call option in the transformed
% variable z = S / m, where m is the running minimum of S.
%
% ----------------------------------------------------------
% 1) PROBLEM SETUP AND PDE DERIVATION
%
% For a floating-strike call, the payoff at T is:
%     Payoff = S(T) - m   (assuming S >= m)
%
% Define z = S / m  =>  S = m*z,  and let
%     V(S, m, t) = m * u(z, t).
%
% In the usual risk-neutral setting (with dividend yield q), S follows
%     dS = (r - q)*S dt + sigma*S dW,
% and for S > m, the minimum m is constant. One shows (via Itô’s lemma
% and standard arguments) that V(S,m,t) satisfies the Black–Scholes PDE
% for S>m. Substituting V = m*u(z,t), S = m*z, and m constant, we get
%
%     u_t + 0.5*sigma^2 * z^2 * u_{zz} + (r - q)*z*u_z - r*u = 0,
%
% on the domain z >= 1. The terminal condition at t = T comes from
% the payoff:
%
%     V(S,m,T) = S - m  =>  m*u(z,T) = m*(z - 1)  =>  u(z,T) = z - 1.
%
% Boundary conditions in z-space are typically:
%   - At z=1,  u(1,t) = 0,   (since S=m => payoff ~ 0)
%   - For large z,  u(z,t) ~ (z - 1).  In practice, we truncate at
%       z = zMax  and set  u(zMax,t) = zMax - 1.
%
% ----------------------------------------------------------
% 2) CRANK–NICOLSON FINITE-DIFFERENCE SCHEME
%
% We discretize the PDE
%
%   u_t = 0.5*sigma^2*z^2*u_{zz} + (r - q)*z*u_z - r*u
%
% backward in time from t=T to t=0 using a uniform grid in z
% (z_i, i=0..M) and uniform time steps (n=0..N).  We store the
% solution in a matrix U(i,n) = u(z_i, t_n).  The standard
% Crank–Nicolson method leads to tridiagonal systems at each time step.
%
% ----------------------------------------------------------
% 3) EXAMPLE USAGE
% We pick some parameters and solve for 0 <= t <= T, 1 <= z <= zMax,
% then plot the final solution (u(z,0)) vs. z and also convert to
% V(S,0) = m*u(z,0) vs. S = m*z.
%
% NOTE: This code is purely demonstrative. In production, you may
% refine the grid, boundary conditions, or the PDE operator.

    close all; clc;

    % -------------------------------
    % Model / Option Parameters
    % -------------------------------
    S0    = 100;      % Current underlying price
    m0    = 80;       % Running minimum so far
    r     = 0.1;     % Risk-free rate
    q     = 0.00;     % Dividend yield
    sigma = 0.4;     % Volatility
    T     = 1.0;      % Time to maturity

    % -------------------------------
    % Numerical Discretization
    % -------------------------------
    zMax  = 5;        % Maximum ratio = Smax/m0, e.g. we assume S can go ~5*m0
    M     = 200;      % Number of z steps
    N     = 200;      % Number of time steps
    dz    = (zMax - 1)/M;  % spacing in z
    dt    = T / N;         % spacing in time

    % Solve PDE with Crank-Nicolson
    [zGrid, Ufinal] = SolveLookbackCallCN(r, q, sigma, T, zMax, M, N);

    % -------------------------------
    % Interpolate to find solution at z0 = S0/m0
    % Then multiply by m0 => the actual option value at S0
    % -------------------------------
    z0   = S0 / m0;
    U0   = interp1(zGrid, Ufinal(:,1), z0, 'linear', 'extrap');
    V0   = m0 * U0;

    fprintf('Option value at S0=%.2f (m=%.2f) is %.4f\n', S0, m0, V0);

    % -------------------------------
    % Plot: 1) u(z,0) vs. z
    %       2) V(S,0) vs. S
    % -------------------------------
    figure(1);
    plot(zGrid, Ufinal(:,1), 'b-', 'LineWidth',1.5);
    xlabel('z = S/m');  ylabel('u(z,0)');  grid on;
    title('Floating-Strike Lookback Call in z-space (u vs. z)');

    % Convert to S = m0*z,  V = m0*u
    Sgrid = m0 * zGrid;
    Vgrid = m0 * Ufinal(:,1);

    figure(2);
    plot(Sgrid, Vgrid, 'r-', 'LineWidth',1.5);
    xlabel('Stock Price S');  ylabel('Option Value');
    title('Floating-Strike Lookback Call vs. S');  grid on;

end  % end of main function



function [zVals, U] = SolveLookbackCallCN(r, q, sigma, T, zMax, M, N)
% SolveLookbackCallCN:
%   Uses Crank–Nicolson to solve the PDE
%      u_t = 0.5*sigma^2 * z^2 * u_{zz}
%            + (r-q)*z * u_z
%            - r * u
%   for 1 <= z <= zMax, and  0 <= t <= T.
%
% Discretization:
%   - z in [1, zMax],  M steps => M+1 points
%   - t in [0, T],     N steps => N+1 time levels
%   - We'll store U(i,n) = u(z_i, t_n).
%   - We'll step backward in time from n=N (t=T) down to n=0 (t=0).
%
% Boundary/terminal conditions:
%   - Terminal: u(z,T) = z - 1,  for  z >= 1.
%   - z=1 boundary: u(1,t) = 0
%   - z=zMax boundary: u(zMax,t) = zMax - 1
%
% Returns:
%   zVals: (M+1)x1 vector of z grid points
%   U    : (M+1)x(N+1) matrix of the solution, so U(:,n+1) = u(z, t_n).

    % Setup grids
    dz = (zMax - 1)/M;
    zVals = linspace(1, zMax, M+1)';   % z_0=1, z_M=zMax

    dt = T / N;
    tVals = linspace(0, T, N+1);  % we will use tVals(N+1)=T

    % Initialize solution matrix
    U = zeros(M+1, N+1);

    % -- Terminal condition at t=T --
    %    u(z,T) = z - 1
    U(:, N+1) = zVals - 1;

    % -- Boundary conditions in z --
    % z=1 => U(1,n) = 0
    % z=zMax => U(M+1,n) = zMax - 1
    U(1,:) = 0;
    U(M+1,:) = zMax - 1;

    % Build coefficient arrays for the PDE operator
    % PDE:  u_t = A[u] = 0.5*sigma^2*z^2*u_{zz} + (r-q)*z*u_z - r*u
    %
    % We'll define discrete i=1..M+1.  The interior points are i=2..M.
    % For each i, define:
    %   alpha_i   = 0.5*sigma^2*z_i^2 / (dz^2)
    %   beta_i(+) = (r-q)*z_i / (2*dz)
    %   beta_i(-) = -(r-q)*z_i / (2*dz)
    %   plus the -r factor on the diagonal.
    %
    % Then in standard 2nd-order finite differences:
    %   u_{zz}(z_i) ~ [u_{i+1} - 2u_i + u_{i-1}] / dz^2
    %   u_z(z_i)   ~ [u_{i+1} - u_{i-1}] / (2 dz)

    iArr = (0:M)';   % so that z_i = 1 + i*dz, but we'll build operator for i=2..M-1

    % Precompute for each i:
    zMid = 1 + iArr*dz;    % same as zVals, but we'll skip boundary i=1 & i=M+1
    alpha = 0.5 * sigma^2 .* (zMid.^2) / (dz^2);
    % For the drift term:
    plusBeta  = ((r - q) .* zMid) / (2*dz);
    minusBeta = -plusBeta;

    % The reaction term -r on the diagonal
    % Summarize operator L:
    %   L[u_i] = alpha_i*(u_{i+1}-2u_i+u_{i-1})
    %          + plusBeta_i*(u_{i+1}) + minusBeta_i*(u_{i-1})
    %          - r*u_i
    %
    % We'll form the tri-di matrix for interior i=2..M.

    % For Crank–Nicolson, define:
    %   A = I - (dt/2)*L
    %   B = I + (dt/2)*L
    % so that:  A*u^{n} = B*u^{n+1}.

    % Build the interior dimension = M-1
    mainDiag = zeros(M-1,1);
    lowerDiag= zeros(M-2,1);
    upperDiag= zeros(M-2,1);

    for i = 2:M
        iLoc = i-1;  % local index in [1..M-1]
        a_i  = alpha(i);
        b_p  = plusBeta(i);
        b_m  = minusBeta(i);

        % Coeff for PDE: L[u_i]
        %   from alpha*(u_{i+1}-2u_i+u_{i-1})
        %     => +alpha for u_{i+1}, -2alpha for u_i, +alpha for u_{i-1}
        %   from plusBeta*(u_{i+1}) + minusBeta*(u_{i-1})
        %     => +b_p for u_{i+1}, +b_m for u_{i-1}
        %   from -r*u_i
        % => net for u_{i-1}: (alpha + b_m)
        % => net for u_i:    (-2 alpha - r)
        % => net for u_{i+1}:(alpha + b_p)

        c_im1 = a_i + b_m;       % lower
        c_i   = -2*a_i - r;      % main
        c_ip1 = a_i + b_p;       % upper

        mainDiag(iLoc) = c_i;
        if iLoc>1
            lowerDiag(iLoc-1) = c_im1;
        end
        if iLoc<(M-1)
            upperDiag(iLoc) = c_ip1;
        end
    end

    % Build L as a sparse tridiagonal matrix
    Lmat = spdiags([lowerDiag mainDiag upperDiag], [-1 0 1], M-1, M-1);

    % Build A, B for Crank–Nicolson
    Iint = speye(M-1);
    A = Iint - 0.5*dt * Lmat;
    B = Iint + 0.5*dt * Lmat;

    % We solve from n=N down to n=0
    for n = N:-1:1
        % The vector on the interior is U(2..M, n+1)
        uOld = U(2:M, n+1);

        % Right side = B*uOld + boundary adjustments
        rhs = B * uOld;

        % Boundary adjustments come from the known Dirichlet values at i=1 and i=M+1:
        %   i=1   => U(1,n)   = 0
        %   i=M+1 => U(M+1,n) = zMax - 1
        % Each influences the neighbors in the finite-diff stencils.
        %
        % We need to see how L includes i=1 and i=M+1 in the row i=2 and i=M, respectively.
        % Because the PDE for i=2 references i=1, etc.
        %
        % PDE row i=2 => depends on i=1 => c_im1*(U(1))
        % PDE row i=M => depends on i=M+1 => c_ip1*(U(M+1))
        %
        % For the CN approach: the boundary terms appear in both B*uOld and A*uNew,
        % so we add them to the RHS with the known boundary values from time n+1 and n.
        %
        % The coefficient c_im1 for row i=2 is: alpha(2)+minusBeta(2).
        % The coefficient c_ip1 for row i=M is: alpha(M)+plusBeta(M).
        % We'll handle them similarly.

        % Row i=2 => local index 1 => lowerDiag(1) => c_im1(2)
        c_i2 = alpha(2) + minusBeta(2);
        boundaryVal1 = U(1,n) + U(1,n+1);  % appear in B and A
        rhs(1) = rhs(1) - 0.5*dt*c_i2*boundaryVal1;

        % Row i=M => local index M-1 => upperDiag(M-2) => c_ip1(M)
        c_iM = alpha(M) + plusBeta(M);
        boundaryValM = U(M+1,n) + U(M+1,n+1); 
        rhs(end) = rhs(end) - 0.5*dt*c_iM*boundaryValM;

        % Solve A*uNew = rhs
        uNew = A \ rhs;

        % Store back into U(2..M,n)
        U(2:M, n) = uNew;
    end

end
