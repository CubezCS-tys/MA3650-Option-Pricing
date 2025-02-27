function [price, gridS, matval] = CrankNicholsonPowerOptionTest( ...
    S0, K, r, T, sigma, p, Smax, dS, dt, optionType)
% CrankNicholsonPowerOption:
%   Prices a European power call/put using the Crank–Nicolson finite difference method.
%
%   Payoff:
%       - Power Call:  payoff = max(S^p - K, 0)
%       - Power Put:   payoff = max(K - S^p, 0)
%
%   INPUTS:
%       S0         : current underlying price
%       K          : strike
%       r          : risk-free rate
%       T          : time to maturity (in years)
%       sigma      : volatility
%       p          : exponent (S^p)
%       Smax       : maximum underlying price on the grid
%       dS         : spatial step size
%       dt         : time step size
%       optionType : 'call' or 'put'
%
%   OUTPUTS:
%       price      : the option value at t=0, S=S0
%       gridS      : vector of underlying prices, size (M+1) x 1
%       matval     : matrix of option values.  Each column is a time step,
%                    matval(:,1) is time=0, matval(:,end) is time=T

    % 1) Set up the grid sizes, adjusting if necessary
    M = round(Smax / dS);   % number of spatial steps
    dS = Smax / M;          % possibly updated
    N = round(T / dt);      % number of time steps
    dt = T / N;             % possibly updated
    
    % 2) Create storage for the option values
    %    matval(i,j) will be the option value at node (i, j)
    %      i -> index for S, j -> index for time
    matval = zeros(M+1, N+1);
    
    % Underlying price grid: S_i = i * dS, i=0..M
    gridS = linspace(0, Smax, M+1)';  
    
    % Time indices: j=0..N, each step is dt
    %   t_j = j * dt
    % We solve backwards in time from j=N (t=T) to j=0 (t=0)
    
    % 3) Final (terminal) condition at t = T
    switch lower(optionType)
        case 'call'
            % payoff = max(S^p - K, 0)
            matval(:, N+1) = max(gridS.^p - K, 0);
        case 'put'
            % payoff = max(K - S^p, 0)
            matval(:, N+1) = max(K - gridS.^p, 0);
        otherwise
            error('optionType must be either "call" or "put".');
    end
    
    % 4) Boundary conditions for 0 <= t <= T
    %    We'll define BC at S=0 (i=0) and at S=Smax (i=M).
    
    % For a power CALL: 
    %   - Lower boundary near S=0 => payoff ~ 0 if p > 0
    %   - Upper boundary near S=Smax => payoff ~ (Smax^p - K)*exp(-r*(T - t)) 
    %     if Smax^p >> K.  Or we can forcibly set it if we know Smax^p >> K
    %
    % For a power PUT:
    %   - Lower boundary near S=0 => payoff ~ K*exp(-r*(T - t))
    %   - Upper boundary near S=Smax => payoff ~ 0 if Smax^p >> K.
    
    for j = 0:N
        tau = j * dt;        % "forward" time, but we'll solve backwards
        t_j = T - tau;       % actual time from maturity
        
        switch lower(optionType)
            case 'call'
                % S=0 => payoff ~ max(0^p - K, 0) = 0, 
                %   at earlier times ~ 0 as well
                matval(1,  N+1-j) = 0;  
                
                % S=Smax => if Smax^p >> K, the option is deep in the money
                %   approximate => (Smax^p - K) * e^{-r * (T - t_j)}
                bigVal = (Smax^p - K) * exp(-r * (t_j));
                if bigVal < 0, bigVal = 0; end % ensure non-neg
                matval(M+1, N+1-j) = bigVal;
                
            case 'put'
                % S=0 => payoff ~ K at maturity, 
                %   discount => K * e^{-r * (T - t_j)} 
                matval(1, N+1-j) = K * exp(-r * (t_j));
                
                % S=Smax => if Smax^p >> K, payoff is ~0
                matval(M+1, N+1-j) = 0;
        end
    end
    
    % 5) Set up the CN coefficients alpha, beta, gamma.
    % We use i=0..M, so "interior" is i=1..M-1.  
    % Usually, i refers to S_i = i * dS.
    % alpha_i = 0.25 dt [ sigma^2 i^2 - r i ]
    % beta_i  = -0.5 dt [ sigma^2 i^2 + r ]
    % gamma_i = 0.25 dt [ sigma^2 i^2 + r i ]
    %
    % But we build them as vectors from i=0..M, then in the matrix we use i=2..M etc.
    
    i_vec  = (0:M)';
    alpha  = 0.25 * dt * ( sigma^2 * (i_vec.^2) - r * i_vec );
    beta   = -0.5 * dt * ( sigma^2 * (i_vec.^2) + r );
    gamma  = 0.25 * dt * ( sigma^2 * (i_vec.^2) + r * i_vec );
    
    % Build M1 and M2 (tridiagonal) for i=2..M
    % M1 * V^j = M2 * V^{j+1} + boundary terms
    M1 =  diag(1 - beta(2:M)    ) ...
        - diag(alpha(3:M),-1)  ...
        - diag(gamma(2:M-1), 1);
    
    M2 =  diag(1 + beta(2:M)    ) ...
        + diag(alpha(3:M),-1)   ...
        + diag(gamma(2:M-1), 1);
    
    % Factor M1 = L * U for fast solves
    [L,U] = lu(M1);
    
    % 6) Main time loop: solve backward in time
    for j = N:-1:1
        % "aux" is the known boundary contribution for the interior
        aux = zeros(M-1,1);
        
        % Add the known boundary from i=1 (S=0 side) -> alpha(2)* [V(1,j)+V(1,j+1)]
        aux(1)   = alpha(2)  * ( matval(1,j) + matval(1,j+1) );
        
        % Add the known boundary from i=M (S=Smax side) -> gamma(M)* [V(M+1,j)+V(M+1,j+1)]
        aux(end) = gamma(M)  * ( matval(M+1,j) + matval(M+1,j+1) );
        
        % Right-hand side = M2 * V^{j+1} + boundary
        rhs = M2 * matval(2:M, j+1) + aux;
        
        % Solve M1 * V^j = rhs
        matval(2:M, j) = U \ (L \ rhs);
    end
    
    % 7) The time-0 value is in matval(:,1).  
    % Interpolate to find the value at S0 (assuming S0 <= Smax).
    price = interp1(gridS, matval(:,1), S0, 'linear');
end
