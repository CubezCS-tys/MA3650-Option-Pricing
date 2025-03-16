function [price, Xgrid, Uall] = LookbackCallFloating_CN(r, sigma, T, Xmax, Nx, Nt)
% LOOKBACKCALLFLOATING_CN Crank–Nicolson for a floating‐strike lookback call
%   dimension‐reduced PDE: u_t = 0.5*sigma^2 * X^2 * u_{XX} + r*X*u_X - r*u.
%   Domain: X in [1, Xmax],  t in [0, T].
%   Terminal payoff: u(T,X) = X - 1.
%   BCs: u(t,1)=0, and at Xmax we set u(t,Xmax)=Xmax-1 (Dirichlet).

    % ---- 1) Grid in X and t ----
    dx = (Xmax - 1) / Nx;
    Xgrid = linspace(1, Xmax, Nx+1)';  % Nx+1 points => i=1..(Nx+1)
    dt = T / Nt;

    % We'll store U(i,n): i=1..(Nx+1), n=1..(Nt+1).
    Uall = zeros(Nx+1, Nt+1);

    % ---- 2) Terminal condition at t=T ----
    % Uall(:,Nt+1) is u(T, X)
    for i = 1:(Nx+1)
        Uall(i, Nt+1) = max(Xgrid(i) - 1, 0); 
    end

    % ---- 3) Set up PDE coefficients on interior (i=2..Nx) ----
    iVec = (2:Nx)';
    Xi   = Xgrid(iVec);  % interior X-values

    alpha = 0.5 * sigma^2 .* (Xi.^2) / dx^2;  % for second derivative
    beta  = 0.5 * r .* Xi / dx;              % for first derivative

    % These form the standard 3‐point FD operator for: 0.5*sigma^2*X^2*u_{XX} + r*X*u_X - r*u
    a = alpha - beta;            % sub‐diag
    b = -2*alpha - r;            % main diag
    c = alpha + beta;            % super‐diag

    % ---- 4) Build Crank–Nicolson Matrices A and B ----
    %   A = I - 0.5*dt * L
    %   B = I + 0.5*dt * L
    mainA = 1 - 0.5*dt*b; 
    subA  = -0.5*dt * a(2:end);
    supA  = -0.5*dt * c(1:end-1);

    mainB = 1 + 0.5*dt*b;
    subB  =  0.5*dt * a(2:end);
    supB  =  0.5*dt * c(1:end-1);

    Amat = diag(mainA) + diag(subA, -1) + diag(supA, 1);
    Bmat = diag(mainB) + diag(subB, -1) + diag(supB, 1);

    % Factor A once (size: (Nx-1)x(Nx-1))
    [Lfactor, Ufactor] = lu(Amat);

    % ---- 5) Step backward in time: n=Nt->1 ----
    for n = Nt:-1:1
        Unp1 = Uall(:, n+1);     % known at t_{n+1}
        Un   = Uall(:, n);       % we will update interior (2..Nx)

        % Right‐hand side for interior
        rhs = Bmat * Unp1(2:Nx);

        % Boundary correction: 
        %   lower boundary i=1 => u=0
        %   top boundary i=Nx+1 => u=Xmax-1
        % At i=2, subA(1) references U(1), etc.
        % So we add 0.5*dt*a(1)*(Un(1)+Unp1(1)) to rhs(1) if needed:
        rhs(1) = rhs(1) + 0.5*dt * a(1) * (Un(1) + Unp1(1));  
        % But Un(1) is zero, so effectively that term is 0.

        % For i=Nx, superA(end) references U(Nx+1):
        BCtopNow  = (Xmax - 1);  % at time t_n
        BCtopNext = (Xmax - 1);  % at time t_{n+1}
        rhs(end) = rhs(end) + 0.5*dt * c(end) * (BCtopNow + BCtopNext);

        % Solve system A*Un_interior = rhs
        Un(2:Nx) = Ufactor \ (Lfactor \ rhs);

        % Reapply boundary conditions in space
        Un(1) = 0;               % X=1 => u=0
        Un(Nx+1) = Xmax - 1;     % X=Xmax => Dirichlet
        Uall(:, n) = Un;
    end

    % ---- 6) Final price if S(0)=m(0) => X(0)=1 => index i=1 ----
    price = Uall(1,1);
end


%% ----------------------------
%  Example of calling it:
%  (Put this in a separate script, e.g. "Test3.m".)
%
%  r      = 0.05;
%  sigma  = 0.2;
%  T      = 1.0;
%  Xmax   = 5.0;     % Must be > 1
%  Nx     = 200;
%  Nt     = 200;
%
%  [u0, Xgrid, U] = LookbackCallFloating_CN(r, sigma, T, Xmax, Nx, Nt);
%
%  % If S0 = m0 initially, the floating-strike lookback call value is:
%  S0        = 100;
%  callPrice = S0 * u0;
%  disp(['Floating-strike call = ', num2str(callPrice)]);
