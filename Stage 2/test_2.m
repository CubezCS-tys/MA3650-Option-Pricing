% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % Crank–Nicolson for a Standard European Call under Black–Scholes
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % clear; clc;
% % 
% % % 1) Parameters
% % Smax   = 2000;    % Maximum underlying price on the grid
% % K      = 150;    % Strike
% % r      = 0.05;   % Risk-free interest rate
% % sigma  = 0.2;    % Volatility
% % T      = 1.0;    % Maturity (in years)
% % N      = 2000;    % Number of S-steps
% % M      = 200;    % Number of time steps
% % 
% % dS = Smax / N;         % Spatial step size
% % dt = T / M;            % Time step size
% % S  = linspace(0, Smax, N+1).';  % Underlying price grid as a column vector
% % 
% % % We march backward in time: t_m = m*dt, m = 0..M, so that t_M = T.
% % 
% % % 2) Initialize the solution array
% % % V(i,m) = option value at S_i, time t_m
% % V = zeros(N+1, M+1);
% % 
% % % 3) Final (terminal) condition at t = T
% % % For a call: payoff = max(S - K, 0)
% % for i = 1:N+1
% %     V(i, M+1) = max(S(i) - K, 0);
% % end
% % 
% % % 4) Build coefficients for Crank–Nicolson
% % % The Black–Scholes PDE for a European call:
% % %   dV/dt + 0.5 * sigma^2 * S^2 * d^2V/dS^2 + r*S * dV/dS - r*V = 0
% % % Discretize S into i=0..N. We'll solve for i=1..N-1 (interior points).
% % 
% % % Arrays to hold the tridiagonal coefficients
% % a = zeros(N-1,1);
% % b = zeros(N-1,1);
% % c = zeros(N-1,1);
% % 
% % % For i=1..N-1, let S_i = i*dS.
% % % Typically in the standard CN scheme, you have:
% % %   alpha_i = 0.5 * dt * ( sigma^2 * S_i^2 / dS^2 - r * S_i / dS )
% % %   beta_i  = - dt * ( 0.5*sigma^2 * S_i^2 / dS^2 + 0.5*r )
% % %   gamma_i = 0.5 * dt * ( sigma^2 * S_i^2 / dS^2 + r * S_i / dS )
% % 
% % for i = 1:(N-1)
% %     Si       = i * dS;
% %     alpha_i  = 0.5 * dt * ( sigma^2 * Si^2 / dS^2 - r * Si / dS );
% %     beta_i   = - dt * ( 0.5 * sigma^2 * Si^2 / dS^2 + 0.5*r );
% %     gamma_i  = 0.5 * dt * ( sigma^2 * Si^2 / dS^2 + r * Si / dS );
% %     
% %     % For CN, we form two matrices Aplus and Aminus:
% %     %   (I + 0.5*A) and (I - 0.5*A), where A is the matrix from PDE.
% %     % We'll store the diagonals for them in a, b, c (for Aplus) and use +/- in the loop.
% %     
% %     % "Aplus" diag entries => +0.5*(alpha_i, beta_i, gamma_i)
% %     % "Aminus" diag entries => -0.5*(alpha_i, beta_i, gamma_i)
% %     
% %     % But it's usually simpler to store the "centered" version here and
% %     % then add or subtract in the iteration loop. For clarity, let's just
% %     % store the partial parts and create the final tri-di below:
% %     a(i) = alpha_i;   % sub-diagonal
% %     b(i) = beta_i;    % main diagonal
% %     c(i) = gamma_i;   % super-diagonal
% % end
% % 
% % % Build the tridiagonal matrices for CN:
% % Aplus  = zeros(N-1,N-1);
% % Aminus = zeros(N-1,N-1);
% % 
% % for i = 1:(N-1)
% %     % main diagonal
% %     Aplus(i,i)  = 1 - 0.5*b(i);
% %     Aminus(i,i) = 1 + 0.5*b(i);
% % end
% % 
% % for i = 1:(N-2)
% %     % sub-diagonal
% %     Aplus(i+1,i)  = -0.5*a(i+1);
% %     Aminus(i+1,i) =  0.5*a(i+1);
% %     
% %     % super-diagonal
% %     Aplus(i,i+1)  = -0.5*c(i);
% %     Aminus(i,i+1) =  0.5*c(i);
% % end
% % 
% % % 5) Boundary conditions
% % % For each time step, we need to set V(0,m) and V(N,m).
% % % - For a call, V(0,t) = 0  (cannot be worth more than zero if S=0).
% % % - As S->Smax, approximate V(Smax,t) ~ Smax - K*exp(-r*(T-t)) for large S.
% % 
% % % Impose them at each time step. 
% % % We'll do it inside the time loop or right after we solve the linear system.
% % 
% % V(1,:) = 0;  % at S=0, call is worthless
% % % At S=Smax => set a "linear in S" boundary with discount (rough approximation).
% % % More precisely: V(Smax,t) = Smax - K e^{-r(T-t)}
% % % We'll do something slightly simpler: keep it as Smax - K*exp(-r*(T - t_m)) 
% % % updated each time step.
% % 
% % % 6) Step back in time via Crank–Nicolson
% % % We already have V(:,M+1) = payoff at t=T.
% % % Now we move backward to m=M, M-1, ... , 1.
% % 
% % for m = M:-1:1
% %     % Current time = t_m = (m-1)*dt
% %     % For boundary condition at S=Smax, we set:
% %     t_m = (m-1)*dt;
% %     V(N+1,m) = Smax - K*exp(-r*(T - t_m));
% %     
% %     % Right-hand side for interior points:
% %     rhs = Aminus * V(2:N,m+1);
% %     
% %     % Incorporate boundary into rhs:
% %     %   sub-di affects V(1), super-di affects V(N+1).
% %     %   In "Aplus x = rhs", the first row depends on V(1),
% %     %   and the last row depends on V(N+1).
% %     % Let’s name them B0 and Bmax for boundary at S=0 and S=Smax:
% %     B0   = V(1,  m);     % lower boundary
% %     Bmax = V(N+1,m);     % upper boundary
% %     
% %     % sub-di contribution (row index i=1 => second node in domain)
% %     rhs(1)   = rhs(1)   + 0.5*a(1)* B0;
% %     % super-di contribution (row index i=N-1 => last interior node)
% %     rhs(end) = rhs(end) + 0.5*c(N-1)* Bmax;
% %     
% %     % Solve the linear system Aplus * V(2:N,m) = rhs
% %     V(2:N,m) = Aplus \ rhs;
% % end
% % 
% % % 7) The Option Value at t=0
% % % The array V(:,1) now holds the option value at time t=0 for S in [0,Smax].
% % % If you want the call price at S0 (e.g. S0 = K or another spot), pick the index.
% % 
% % S0    = K;                  % Example: price at S=K
% % i_S0  = round(S0 / dS) + 1; % +1 because indexing starts at 1 in MATLAB
% % callPrice = V(i_S0,1);
% % 
% % fprintf('----------------------------------\n');
% % fprintf('Crank–Nicolson Call Price at S=K: %f\n', callPrice);
% % fprintf('----------------------------------\n');
% % 
% % % 8) (Optional) Plot the value at t=0
% % figure;
% % plot(S, V(:,1), 'b-', 'LineWidth', 2);
% % xlabel('Underlying Price S');
% % ylabel('Option Value V(S,0)');
% % title('European Call: Crank–Nicolson, t=0');
% % grid on;

aux = zeros(M-1)