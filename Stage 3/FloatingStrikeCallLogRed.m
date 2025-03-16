function [price, vetz, matval] = FloatingStrikeCallLogRed(S0, Min, r, T, sigma, Smax, dz, dt)

    % Ensure Smax is large enough for price movement
    Smax = max(Smax, Min * exp(4 * sigma * sqrt(T))); % Extend max stock price
    
    % Logarithmic transformation
    zmin = log(1);              % Lower boundary (corresponding to S=Min*exp(0) = Min)
    zmax = log(Smax / Min);     % Upper boundary in log space
    
    % Ensure M is at least 3 for proper indexing
    M = max(round((zmax - zmin) / dz), 3);
    dz = (zmax - zmin) / M;
    
    % Ensure N is reasonable for time-stepping
    N = max(round(T / dt), 10);
    dt = T / N;
    
    % Debugging output
    disp(['Grid Size M = ', num2str(M), ', Time Steps N = ', num2str(N)]);
    
    % Logarithmic grid
    vetz = linspace(zmin, zmax, M+1)';
    vetS = Min * exp(vetz);  %#ok<NASGU>  % Real-space stock grid if needed
    veti = 0:M;
    vetj = 0:N;
    
    % Initialize solution matrix
    matval = zeros(M+1, N+1);
    
    % Boundary conditions for the payoff and edges
    % At t = T, payoff = max(S - K, 0) but here floating strike => max(S/avg - 1,0).
    % The code as given uses max(e^z - 1, 0), so:
    matval(:, N+1) = max(exp(vetz) - 1, 0);  % Terminal payoff in transformed space
    
    % Lower boundary condition (zmin => S=Min)
    matval(1, :) = 0; 
    
    % Upper boundary condition (zmax => large S).  Typically you'd do something
    % like V ~ S - PV(K), but here is a version the code uses:
    matval(M+1, :) = exp(zmax) - exp(-r * (T - dt * vetj)); 
    
    % Coefficients for Crank-Nicolson in log space
    % — Make them vectors so alpha(3:M), etc. are valid.
    alpha_val = 0.25 * dt * (sigma^2 / dz^2 - (r - 0.5 * sigma^2) / (2 * dz));
    beta_val  = -0.5 * dt * (sigma^2 / dz^2 + r);
    gamma_val = 0.25 * dt * (sigma^2 / dz^2 + (r - 0.5 * sigma^2) / (2 * dz));
    
    alpha = alpha_val * ones(M+1, 1);
    beta  = beta_val  * ones(M+1, 1);
    gamma = gamma_val * ones(M+1, 1);
    
    % Build matrices M1 and M2 for the Crank-Nicolson scheme
    % These are (M-1) x (M-1) matrices for the interior nodes (i=2..M)
    main_diag_M1  = 1 - beta(2:M);
    sub_diag_M1   = -alpha(3:M);
    super_diag_M1 = -gamma(2:M-1);
    M1 = diag(sub_diag_M1, -1) + diag(main_diag_M1) + diag(super_diag_M1, 1);
    
    [L, U] = lu(M1);
    
    main_diag_M2  = 1 + beta(2:M);
    sub_diag_M2   = alpha(3:M);
    super_diag_M2 = gamma(2:M-1);
    M2 = diag(sub_diag_M2, -1) + diag(main_diag_M2) + diag(super_diag_M2, 1);
    
    % Solve backward in time
    aux = zeros(M-1, 1);
    for j = N:-1:1
        
        % We have M-1 interior nodes, so aux is length M-1
        % First entry picks up boundary at i=1, last picks up boundary at i=M+1
        aux(1)     = alpha(2)  * (matval(1,j)   + matval(1,j+1));
        aux(end)   = gamma(M)  * (matval(M+1,j) + matval(M+1,j+1));
        
        RHS = M2 * matval(2:M, j+1) + aux;
        matval(2:M, j) = U \ (L \ RHS);
    end
    
    % Return the price by interpolating along z at t=0
    price = Min * interp1(vetz, matval(:,1), log(S0 / Min), 'linear');
end
