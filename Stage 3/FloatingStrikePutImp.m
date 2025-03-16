% % % function [price, vetS, matval] = FloatingStrikePutImp(S0,Max,r,T,sigma,Smax,dS,dt)
% % % %set up grid and adjust increments if necessary
% % % Sz = S0/Max;
% % % zmin = Smax/Max;
% % % M = round(zmin/dS);
% % % dz = zmin/M;
% % % N = round(T/dt);
% % % dt = T/N;
% % % matval = zeros(M+1,N+1);
% % % vetS = linspace(0,Smax,M+1)'; % Create a column vector (vetS) containing M+1 equally spaced asset prices from 0 to Smax.
% % % vetz = vetS./Max;
% % % veti = 0:M;
% % % vetj = 0:N;
% % % % set up boundary conditions
% % % matval(:,N+1) = max(1-vetz,0);
% % % matval(1,:) = exp(-r *(T-dt*vetj));
% % % % matval(M+1, :) = matval(M, :) / (1 - dz);
% % % matval(M+1, :) = 0;
% % % a = 0.5*(r*dt*veti-sigma^2*dt*(veti.^2));
% % % b = 1+sigma^2*dt*(veti.^2)+r*dt;
% % % c = -0.5*(r*dt*veti+sigma^2*dt*(veti.^2));
% % % coeff = diag(a(3:M),-1) + diag(b(2:M)) + diag(c(2:M-1),1);
% % % [L,U] = lu(coeff);
% % % % solve the sequence of linear systems
% % % aux = zeros(M-1,1);
% % % for j=N:-1:1
% % %    aux(1) = a(2) * matval(M+1,j);
% % %    matval(2:M,j) = U \ (L \ (matval(2:M,j+1) + aux));
% % % end
% % % 
% % % price = Max * interp1(vetz, matval(:,1), Sz);
% % function [price, vetS, matval] = FloatingStrikePutImp(S0, Max, r, T, sigma, Smax, dS, dt)
% % %FLOATINGSTRIKEPUTIMP
% % %   Prices a floating-strike put via a dimensionless approach: z = S / Max.
% % %
% % %   S0   = current asset price
% % %   Max  = reference maximum for the floating-strike
% % %   r    = risk-free rate
% % %   T    = maturity
% % %   sigma= volatility
% % %   Smax = spatial truncation for asset price
% % %   dS   = desired spacing in the ratio (or in S, adjusted to ratio)
% % %   dt   = time step
% % 
% %     %--------------------------
% %     % 1) Setup ratio grid z = S/Max
% %     %--------------------------
% %     zmax = Smax / Max;          % domain in ratio space
% %     M    = round(zmax / dS);    % number of steps in ratio
% %     dz   = zmax / M;            % adjusted ratio step
% %     N    = round(T / dt);       % number of time steps
% %     dt   = T / N;               % adjusted time step
% % 
% %     % Storage for PDE solution: size (M+1) x (N+1)
% %     matval = zeros(M+1, N+1);
% % 
% %     % Real asset prices & ratio:
% %     vetS = linspace(0, Smax, M+1)';  % from 0 to Smax
% %     vetz = vetS ./ Max;             % ratio z = S/Max
% % 
% %     % Indices for referencing
% %     veti = (0 : M)';
% %     vetj = (0 : N)';
% % 
% %     %--------------------------
% %     % 2) Terminal Condition (t = T)
% %     %    Floating-strike put => payoff = max(M - S, 0).
% %     %    Dimensionless: u(T,z) = max(1 - z, 0).
% %     %--------------------------
% %     matval(:, N+1) = max(1 - vetz, 0);
% % 
% %     %--------------------------
% %     % 3) Boundary Conditions in space
% %     %    a) z=0 (S=0): The put can be as high as Max,
% %     %       so dimensionless u(0) ~ e^{-r(T - t)}.
% %     %    b) z=zmax => S=Smax => large S => put worthless => u=0.
% %     %--------------------------
% %     % Lower boundary: i=1 => z=0
% %     matval(1, :) = exp(-r * (T - dt .* vetj));   % dimensionless
% % 
% %     % Upper boundary: i=M+1 => z=zmax => worthless
% %     matval(M+1, :) = 0;
% % 
% %     %--------------------------
% %     % 4) Coefficients for fully implicit PDE in ratio form
% %     %    PDE: u_t = 0.5*sigma^2 z^2 u_{zz} + (r z) u_z - r u
% %     %    with a standard 3‐point FD scheme => a, b, c.
% %     %--------------------------
% %     a = 0.5*( r*dt.*veti - sigma^2*dt.*(veti.^2) );
% %     b = 1 + (sigma^2*dt.*(veti.^2)) + r*dt;
% %     c = -0.5*( r*dt.*veti + sigma^2*dt.*(veti.^2) );
% % 
% %     % Build tridiagonal matrix for interior i=2..M
% %     coeff = diag(a(3:M), -1) + diag(b(2:M)) + diag(c(2:M-1), 1);
% % 
% %     % Factor for repeated solves
% %     [L,U] = lu(coeff);
% % 
% %     %--------------------------
% %     % 5) Backward time-stepping
% %     %--------------------------
% %     aux = zeros(M-1,1);
% %     for j = N : -1 : 1
% % 
% %         % Boundary from i=1 => z=0 => appears in the first row of the system
% %         % This row references 'a(2)*u(1,j)' (sub‐diagonal). Usually the sign is negative:
% %         aux(1) = -a(2) * matval(1,j);
% % 
% %         % Boundary from i=M+1 => z=zmax => appears in the last row => 'c(M)*u(M+1,j)'
% %         aux(M-1) = -c(M) * matval(M+1,j);
% % 
% %         % Solve tri‐diag system for interior (2..M)
% %         rhs = matval(2:M, j+1) + aux;
% %         matval(2:M, j) = U \ (L \ rhs);
% %     end
% % 
% %     %--------------------------
% %     % 6) Price at S0
% %     %    If S0=Max * z0 => z0 = S0/Max
% %     %    actual put value = Max * u(0,z0)
% %     %--------------------------
% %     z0 = S0 / Max;
% %     price = Max * interp1(vetz, matval(:,1), z0, 'linear');
% % end
% 
% 
% function [price, vetS, matval] = FloatingStrikePutImp(S0, Max, r, T, sigma, Smax, dS, dt)
% %FLOATINGSTRIKEPUTIMP
% %   Prices a floating-strike put via a dimensionless approach: z = S / Max.
% %
% %   S0    = Current asset price
% %   Max   = Running maximum reference for the floating-strike
% %   r     = Risk-free rate
% %   T     = Time to maturity
% %   sigma = Volatility
% %   Smax  = Maximum asset price for spatial truncation
% %   dS    = Desired spacing in ratio (or S) grid
% %   dt    = Time-step size
% %
% %   This code solves the PDE in the ratio z = S / Max:
% %   u_t = 0.5*sigma^2*z^2*u_{zz} + (r*z)*u_z - r*u,  for  0 <= z <= zmax,
% %   with payoff at t = T: u(T,z) = max(1 - z, 0).
% 
%     %--------------------------
%     % 1) Setup the ratio grid z = S/Max
%     %--------------------------
%     zmax = Smax / Max;          % domain in ratio space
%     M    = round(zmax / dS);    % number of ratio steps
%     dz   = zmax / M;            % adjusted ratio spacing
%     N    = round(T / dt);       % number of time steps
%     dt   = T / N;               % adjusted time step
% 
%     % Prepare the solution matrix: size (M+1)x(N+1)
%     matval = zeros(M+1, N+1);
% 
%     % Real stock-price grid and ratio grid
%     vetS = linspace(0, Smax, M+1)';  % from 0 to Smax
%     vetz = vetS ./ Max;             % ratio z = S/Max
% 
%     % Index vectors
%     veti = (0:M)';  
%     vetj = (0:N)';
% 
%     %--------------------------
%     % 2) Terminal Condition (t=T)
%     %    For a floating-strike put => payoff = max(Max - S, 0).
%     %    Dimensionless form => u(T,z) = max(1 - z, 0).
%     %--------------------------
%     matval(:, N+1) = max(1 - vetz, 0);
% 
%     %--------------------------
%     % 3) Boundary Conditions in space
%     %    - z=0  => S=0   => put can be up to Max => dimensionless: e^{-r(T-t)}
%     %    - z=zmax => S=Smax => worthless => 0
%     %--------------------------
%     % Lower boundary: i=1 => z=0
%     matval(1, :) = exp(-r * (T - dt .* vetj));  
% 
%     % Upper boundary: i=M+1 => z=zmax => worthless => 0
%     matval(M+1, :) = 0;
% 
%     %--------------------------
%     % 4) Coefficients for Fully Implicit Scheme
%     %    PDE: u_t = 0.5*sigma^2 z^2 u_{zz} + (r*z) u_z - r*u
%     %--------------------------
%     a = 0.5*( r*dt .* veti - sigma^2*dt .* (veti.^2) );
%     b = 1 + sigma^2*dt .* (veti.^2) + r*dt;
%     c = -0.5*( r*dt .* veti + sigma^2*dt .* (veti.^2) );
% 
%     % Build the tridiagonal matrix for interior points i=2..M
%     coeff = diag(a(3:M), -1) + diag(b(2:M)) + diag(c(2:M-1), 1);
% 
%     % Factorize once for repeated solves
%     [L, U] = lu(coeff);
% 
%     %--------------------------
%     % 5) Backward Time-Stepping
%     %    We march from j=N to j=1, using the known terminal condition at j=N+1
%     %--------------------------
%     aux = zeros(M-1, 1);
%     for j = N:-1:1
%         % Lower boundary: i=1 => the sub-diagonal involves a(2)*u(1)
%         aux(1)     = -a(2) * matval(1, j);
% 
%         % Upper boundary: i=M+1 => the super-diagonal involves c(M)*u(M+1)
%         aux(M-1)   = -c(M) * matval(M+1, j);
% 
%         % Right-hand side uses matval(2:M, j+1)
%         rhs = matval(2:M, j+1) + aux;
% 
%         % Solve the linear system for interior nodes 2..M
%         matval(2:M, j) = U \ (L \ rhs);
%     end
% 
%     %--------------------------
%     % 6) Interpolate Price at S0
%     %    If S0=Max*z0 => z0=S0/Max => dimensionless
%     %    Real put value => Max * u(0,z0)
%     %--------------------------
%     z0 = S0 / Max;
%     price = Max * interp1(vetz, matval(:,1), z0, 'linear');
% end
%


