function price = LookbackFloatingCall_CN(S0, M0, r, q, sigma, T, Nx, Xmax, Nt)
% LookbackFloatingCall_CN: Prices a floating-strike lookback call
%   via a one-dimensional Crank-Nicolson finite-difference method.
%
% Inputs:
%   S0    = current spot price
%   M0    = running minimum so far
%   r     = risk-free rate
%   q     = dividend yield
%   sigma = volatility
%   T     = time to maturity
%   Nx    = number of x-grid steps
%   Xmax  = maximum x-value in grid (x = S/M)
%   Nt    = number of time steps
%
% Output:
%   price = option value at time 0 (for S0, M0)

    % 1) Transform S, M -> x = S / M
    x0 = S0 / M0;  % The 'current' x we care about
    
    % 2) Set up the spatial grid in x
    %    We let x range from 1 to Xmax.  x=1 => S=M, the lowest ratio.
    dx = (Xmax - 1) / Nx;
    xGrid = linspace(1, Xmax, Nx+1)';  % column vector: x_0=1, x_Nx=Xmax
    
    % 3) Set up the time grid
    dt = T / Nt;
    tVals = linspace(0, T, Nt+1);  % tVals(Nt+1) = T
    
    % 4) Terminal condition: v(x,T) = x - 1
    %    We'll store v at each time step in a vector V.
    %    V(i) corresponds to v(x_i, t_n).
    vNow = max(xGrid - 1, 0);  % at t = T
    
    % 5) Coefficients for PDE: dv/dt = A(v)
    %    A(v) = 0.5*sigma^2*x^2*v'' + (r-q)*x*v' - r*v
    %    We'll build the finite-difference matrix for the interior nodes.
    
    % Precompute alpha_i, beta_i for each x_i
    % alpha_i => (r-q)*x_i
    % beta_i  => 0.5*sigma^2*x_i^2
    alpha = (r - q) * xGrid;
    beta  = 0.5 * sigma^2 * (xGrid.^2);
    
    % Build tri-di matrices for the operator L(v) in space.
    % We'll do i=2..Nx (interior), with boundary i=1 (x=1) and i=Nx+1 (x=Xmax).
    
    % Initialize arrays for sub-, main-, super-diagonal
    mainDiag = zeros(Nx+1,1);
    lowerDiag = zeros(Nx+1,1);
    upperDiag = zeros(Nx+1,1);
    
    % Finite difference stencils:
    % v''(x_i) ~ [v_{i+1} - 2 v_i + v_{i-1}] / dx^2
    % v'(x_i)  ~ [v_{i+1} - v_{i-1}] / (2 dx)
    %
    % => L(v_i) = beta_i/dx^2 * [v_{i+1} - 2v_i + v_{i-1}] 
    %           + alpha_i/(2 dx)* [v_{i+1} - v_{i-1}] 
    %           - r * v_i
    %
    % We'll fill the interior rows i=2..Nx.  For i=1 and i=Nx+1 we impose BCs.
    
    for i = 2:Nx  % interior points
        xi     = xGrid(i);
        alpha_i = alpha(i);
        beta_i  = beta(i);
        
        % Coeff of v_{i-1}, v_i, v_{i+1} in L(v_i)
        %    L(v_i) = c_{i-1}*v_{i-1} + c_i*v_i + c_{i+1}*v_{i+1}.
        
        c_im1 = beta_i/dx^2 - alpha_i/(2*dx);
        c_i   = -2*beta_i/dx^2 - r;
        c_ip1 = beta_i/dx^2 + alpha_i/(2*dx);
        
        lowerDiag(i) = c_im1;
        mainDiag(i)  = c_i;
        upperDiag(i) = c_ip1;
    end
    
    % Now incorporate boundary conditions:
    % We'll do Dirichlet at x=1 (i=1) and x=Xmax (i=Nx+1).
    % That means v(1,t) and v(Xmax,t) are set from BCs, not from PDE interior.
    %
    % For a floating-strike call:
    %  - At x=1, one might approximate v(1,t)=0 for all t < T 
    %    (since if S=M, the call is near 'intrinsic = 0', ignoring small time value).
    %  - At x=Xmax, we often approximate v(Xmax,t) ~ Xmax - 1 
    %    (large S vs. M => payoff ~ S - M => ratio ~ x => v ~ x - 1).
    %
    % A more refined boundary can be used, but we'll keep it simple here.
    
    % We'll handle these in the time-stepping by setting v(1) and v(Nx+1)
    % after each solve, rather than building them into the matrix rows.
    
    % 6) Construct the Crank-Nicolson matrices:
    % We want: v^{n+1} - v^n = (dt/2)*[ L(v^n) + L(v^{n+1}) ]
    % => (I - dt/2 * L) v^{n+1} = (I + dt/2 * L) v^n
    %
    % We'll build L as a tridiagonal matrix. Then define:
    % A = I - (dt/2)*L
    % B = I + (dt/2)*L
    %
    % Because the boundary is Dirichlet, we’ll fix the top/bottom rows after forming A,B.
    
    Ldiag_main = mainDiag;  % store for building L
    Ldiag_lower = lowerDiag;
    Ldiag_upper = upperDiag;
    
    % Build sparse L
    % interior rows i=2..Nx
    % row i => mainDiag(i), lowerDiag(i), upperDiag(i)
    % We'll keep row 1 and Nx+1 = 0 because they are Dirichlet boundaries.
    
    nnzmax = 3*(Nx+1); 
    rowIdx = [];
    colIdx = [];
    valL   = [];
    
    for i=2:Nx
        % main
        rowIdx = [rowIdx; i];
        colIdx = [colIdx; i];
        valL   = [valL; Ldiag_main(i)];
        % lower
        rowIdx = [rowIdx; i];
        colIdx = [colIdx; i-1];
        valL   = [valL; Ldiag_lower(i)];
        % upper
        rowIdx = [rowIdx; i];
        colIdx = [colIdx; i+1];
        valL   = [valL; Ldiag_upper(i)];
    end
    % Build L as sparse
    Lmat = sparse(rowIdx, colIdx, valL, Nx+1, Nx+1, nnzmax);
    
    % Identity
    Imat = speye(Nx+1);
    
    A = Imat - 0.5*dt * Lmat;
    B = Imat + 0.5*dt * Lmat;
    
    % 7) Backward time stepping: from n=Nt down to n=0
    vOld = vNow;  % at t = T
    for n = Nt:-1:1
        t_n   = tVals(n);   %#ok (not used explicitly, but can be helpful)
        t_n1  = tVals(n+1); %#ok
        
        % Right-hand side = B * vOld
        rhs = B * vOld;
        
        % Enforce boundary conditions in rhs:
        %   v(1) = 0, v(Nx+1) = (xGrid(end) - 1).
        %
        % The matrix multiplication B*vOld won't fix these automatically,
        % so we must adjust the rhs to incorporate Dirichlet conditions in A and B.
        
        % i=1 => v(1)=0
        rhs(1) = 0;  % since A(1,1)*v(1)= v(1), we want v(1) = 0 => modifies RHS
        % Also zero out row 1 of A to ensure v(1) is pinned
        A(1,1) = 1;   A(1,2:end) = 0;  
        
        % i=Nx+1 => v(Nx+1) = xGrid(end) - 1
        bcTop = xGrid(end) - 1;  
        rhs(Nx+1) = bcTop; 
        A(Nx+1,:) = 0; 
        A(Nx+1,Nx+1) = 1;
        
        % Solve A * vNew = rhs
        vNew = A \ rhs;
        
        % Impose boundary conditions explicitly (just to remove any roundoff):
        vNew(1)     = 0;
        vNew(Nx+1)  = bcTop;
        
        % Move on
        vOld = vNew;
    end
    
    % After finishing time steps, vOld = v(x, t=0)
    % We want the option value at x=x0 => we’ll interpolate if x0 not on the grid.
    vAtX0 = interp1(xGrid, vOld, x0, 'linear', 'extrap');
    
    % The actual option price is: V(S0, M0, 0) = M0 * v(x0, 0).
    price = M0 * vAtX0;

end
