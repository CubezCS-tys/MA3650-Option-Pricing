function [price, vetz, matval] = FloatingStrikeCallImplicit(S0, Min, r, T, sigma, Smax, dz, dt)
    % Ensure Smax is large enough for price movement
    Smax = max(Smax, Min * exp(4 * sigma * sqrt(T))); % Extend max stock price
    
    % Logarithmic transformation:
    % We'll treat z = ln(S/Min). Then S = Min * exp(z).
    zmin = log(1);            % This corresponds to S = Min
    zmax = log(Smax / Min);   % Upper boundary in log space
    
    % Ensure M is at least 3 for proper indexing
    M = max(round((zmax - zmin) / dz), 3);
    dz = (zmax - zmin) / M;
    
    % Ensure N is a minimum for stable time-stepping
    N = max(round(T / dt), 10);
    dt = T / N;
    
    % Debugging output
    disp(['[Implicit] Grid Size M = ', num2str(M), ', Time Steps N = ', num2str(N)]);
    
    % Logarithmic grid: z in [zmin, zmax]
    vetz = linspace(zmin, zmax, M+1)';
    % Real stock grid, if needed for reference:
    vetS = Min * exp(vetz); %#ok<NASGU>
    
    % Indices
    veti = 0:M; 
    vetj = 0:N;
    
    % Initialize solution matrix: size (M+1) x (N+1)
    matval = zeros(M+1, N+1);
    
    % Terminal condition at t = T (i.e., j = N+1):
    % For a floating-strike average, your code uses payoff = max(S/K - 1, 0).
    % The snippet simply had max(exp(z) - 1, 0).
    matval(:, N+1) = max(exp(vetz) - 1, 0); 
    
    % Boundary conditions in space for all time-steps:
    % 1) Lower boundary: S=Min => z=zmin => matval(1,:) = 0
    matval(1, :) = 0;
    % 2) Upper boundary: S ~ large => z=zmax => the code uses
    %       V = e^{zmax} - e^{-r (T - something)}
    % (It’s a proxy for "call payoff ~ S - PV(K) at large S".)
    matval(M+1, :) = exp(zmax) - exp(-r * (T - dt*vetj));
    
    % =========================
    % Fully Implicit Discretization
    % =========================
    %
    % PDE in log-space for V(t,z):
    %
    %   dV/dt + 0.5*sigma^2 * d^2V/dz^2
    %   + (r - 0.5*sigma^2)* dV/dz - r*V = 0.
    %
    % Fully implicit in time (Backward Euler):
    %   V^n - V^{n+1} = - dt * [ 0.5*sigma^2 * d^2V^n/dz^2 + ... - r V^n ]
    %
    % We typically form:
    %   (I + dt * L) V^n = V^{n+1},
    % where L is the spatial operator.
    %
    % The standard tridiagonal coefficients for uniform dz:
    %   alpha = 0.5*dt * [ sigma^2/dz^2  - (r - 0.5*sigma^2)/(2 dz) ]
    %   beta  = dt * ( sigma^2/dz^2 + r )
    %   gamma = 0.5*dt * [ sigma^2/dz^2 + (r - 0.5*sigma^2)/(2 dz) ]
    %
    % We'll store them as vectors so alpha(2), alpha(3:M), etc., are valid.
    
    alpha_val = 0.5*dt * (sigma^2/dz^2 - (r - 0.5*sigma^2)/(2*dz));
    beta_val  = dt * (sigma^2/dz^2 + r);
    gamma_val = 0.5*dt * (sigma^2/dz^2 + (r - 0.5*sigma^2)/(2*dz));
    
    alpha = alpha_val * ones(M+1, 1);
    beta  = beta_val  * ones(M+1, 1);
    gamma = gamma_val * ones(M+1, 1);
    
    % Interior grid points: i = 2..M
    % The matrix for (I + dt*L) has dimension (M-1) x (M-1).
    main_diag  = 1 + beta(2:M);
    sub_diag   = -alpha(3:M);
    super_diag = -gamma(2:M-1);
    
    % Build the tridiagonal matrix (call it A).
    A = diag(sub_diag, -1) + diag(main_diag) + diag(super_diag, 1);
    [Lfac, Ufac] = lu(A);  %# LU factor for efficiency
    
    % =============================
    % Backward time-stepping: j=N -> j=1
    % =============================
    for j = N:-1:1
        % Right-hand side = V^{j+1} at interior, plus boundary corrections
        rhs = matval(2:M, j+1);
        
        % Incorporate the boundary conditions at i=1 and i=M+1:
        %   - alpha(2)*matval(1,j)   enters the first row
        %   - gamma(M)*matval(M+1,j) enters the last row
        %
        % Because our matrix has sub/super diagonals labeled i-1 and i+1,
        % the boundary values feed into the first/last row of the system.
        rhs(1)   = rhs(1)   + alpha(2)  * matval(1,   j);
        rhs(end) = rhs(end) + gamma(M)  * matval(M+1, j);
        
        % Solve the tridiagonal system:
        matval(2:M, j) = Ufac \ (Lfac \ rhs);
    end
    
    % =================
    % Final: Interpolate at S0
    % =================
    % We want V(t=0,z) for z=ln(S0/Min)
    z0 = log(S0 / Min);
    price = Min * interp1(vetz, matval(:,1), z0, 'linear');
end
