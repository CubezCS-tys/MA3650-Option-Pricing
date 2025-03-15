

function [price, vetS, matval] = FloatingStrikeCallImp(S0,Min,r,T,sigma,Smax,dS,dt)
% set up grid and adjust increments if necessary
% Sz = S0/Min;
% zmax = Smax/Min;
% M = round(zmax/dS);
% dS = zmax/M;
% N = round(T/dt);
% dt = T/N;
% matval = zeros(M+1,N+1);
% vetS = linspace(0,Smax,M+1)'; % Create a column vector (vetS) containing M+1 equally spaced asset prices from 0 to Smax.
% vetz = vetS./Min;
% veti = 0:M;
% vetj = 0:N;
% % set up boundary conditions
% matval(:,N+1) = max(vetz-1, 0);
% matval(1,:) = 0;
% matval(M+1,:) = (zmax)-exp(-r *(T-dt*vetj));
% 
% % set up the coefficients matrix
% alpha = 0.25*dt*(sigma^2*(veti.^2) - r*veti);
% beta = -dt*0.5*(sigma^2*(veti.^2) + r);
% gamma = 0.25*dt*(sigma^2*(veti.^2) + r*veti);
% M1 = -diag(alpha(3:M) ,-1) + diag(1-beta(2:M)) - diag(gamma(2:M-1) ,1) ;
% [L,U] = lu(M1);
% M2 = diag(alpha(3:M) ,-1) + diag(1+beta(2:M)) + diag(gamma(2:M-1) ,1);
% %solve the sequence of linear systems
% aux = zeros(size(M2,2), 1);
% for j=N:-1:1
%     if length(aux)>1
%         aux(1) = alpha(2) * (matval(1,j)+matval(1,j+1));
%         aux(end) = gamma(end) * (matval(end,j)+matval(end,j+1));
%     else
%         aux = aux(1)+aux(end);
%     end
% 
%     matval(2:M,j) = U \ (L \ ((M2*(matval(2:M,j+1))+aux)));
% end
% 
% %return price, possibly by linear interpolation outside the grid
% 
% price = Min * interp1(vetz, matval(:,1), Sz);

% set up grid and adjust increments if necessary
M = round(Smax/dS);
dS = Smax/M;
N = round(T/dt);
dt = T/N;
matval = zeros(M+1,N+1);
vetS = linspace(0,Smax,M+1)';
veti = 0:M;
vetj = 0:N;
% set up boundary conditions
matval(:,N+1) = max(vetS-Min,0);
matval(1,:) = 0 ;
matval(M+1,:) = Smax-(Min)*exp(-r * dt*(N-vetj));
%set up the tridiagonal coefficients matrix
% alpha = 0.25*dt*(sigma^2*(veti.^2) - r*veti);
% beta = -dt*0.5*(sigma^2*(veti.^2) + r);
% gamma = 0.25*dt*(sigma^2*(veti.^2) + r*veti);
% M1 = -diag(alpha(3:M) ,-1) + diag(1-beta(2:M)) - diag(gamma(2:M-1) ,1) ;
% [L,U] = lu(M1);
% M2 = diag(alpha(3:M) ,-1) + diag(1+beta(2:M)) + diag(gamma(2:M-1) ,1);
% %solve the sequence of linear systems
% aux = zeros(size(M2,2), 1);
% for j=N:-1:1
%     if length(aux)>1
%         aux(1) = alpha(2) * (matval(1,j)+matval(1,j+1));
%         aux(end) = gamma(end) * (matval(end,j)+matval(end,j+1));
%     else
%         aux = aux(1)+aux(end);
%     end
% 
%     matval(2:M,j) = U \ (L \ ((M2*(matval(2:M,j+1))+aux)));
% end
% 
% %return price, possibly by linear interpolation outside the grid
% 
% price =  interp1(vetS, matval(:,1), S0);


a = 0.5*(r*dt*veti-sigma^2*dt*(veti.^2));
b = 1+sigma^2*dt*(veti.^2)+r*dt;
c = -0.5*(r*dt*veti+sigma^2*dt*(veti.^2));
coeff = diag(a(3:M),-1) + diag(b(2:M)) + diag(c(2:M-1),1);
[L,U] = lu(coeff);
% solve the sequence of linear systems
aux = zeros(M-1,1);
for j=N:-1:1
   aux(M-1) = - c(M) * matval(M+1,j);
   matval(2:M,j) = U \ (L \ (matval(2:M,j+1) + aux));
end
% return price, possibly by linear interpolation outside the grid
price = interp1(vetS, matval(:,1), S0);

% function [price, zGrid, matval] = FloatingStrikeCallImp(S0, m, r, T, sigma, Smax, dS, dt)
% % FLOATINGSTRIKECALLIMP  Crank–Nicolson for the floating-strike lookback call
% %   using the transformation  z = S / m ,  where m = min(S_t over [0,t]).
% %
% % INPUTS:
% %   S0    = current underlying price
% %   m     = the running minimum  (so z = S/m)
% %   r     = risk-free rate
% %   T     = time to maturity
% %   sigma = volatility
% %   Smax  = max S for the spatial grid
% %   dS    = mesh size in S
% %   dt    = time step in t
% %
% % OUTPUTS:
% %   price  = the computed option value at (S0,m)
% %   zGrid  = the spatial grid in z
% %   matval = the solution matrix,  size (M+1) x (N+1)
% 
%     %-------------------------
%     % 1) Setup the z-grid
%     %-------------------------
%     zmax = Smax / m;            % since z = S/m
%     M    = round(zmax / dS);    % number of z-steps
%     dz   = zmax / M;            % (possibly slightly different from dS)
% 
%     % time steps
%     N     = round(T / dt);
%     dt    = T / N;
% 
%     % Pre-allocate the solution matrix
%     matval = zeros(M+1, N+1);
% 
%     % Our z values (0..M) correspond to z in [0, zmax], but we shift to [1,zmax].
%     % We'll store index i=1 => z=1, i=M+1 => z=zmax.
%     zVec = linspace(0, zmax, M+1)';
% 
%     % However, the PDE domain is effectively z in [1, zmax].
%     % We'll just keep i=1 => z=0   as a dummy. We'll override boundary i=1 => z=1 in code.
%     % 
%     % For convenience, define iZ = 0..M. We'll interpret i=1 as z=1 in the boundary conditions.
% 
%     %-------------------------
%     % 2) Boundary & Terminal
%     %-------------------------
%     %  - Terminal condition at t=T:  u(z,T) = max(z-1,0)
%     %  - At z=1, for 0 <= t < T: u(1,t)=0
%     %  - At z=zmax:  u(zmax,t) ~ zmax - exp(-r*(T-t))  (large z approx)
% 
%     % Fill terminal condition:
%     for i = 1:M+1
%         zVal = zVec(i);
%         matval(i, N+1) = max(zVal - 1, 0);
%     end
% 
%     % Enforce z=1 => i1 ~ z=1. Find the index closest to z=1:
%     % We do that in the time loop, but set below for t=T as well.
% 
%     % Approx boundary at z=zmax for all times:
%     % We'll overwrite each time step in the loop, but initialize:
%     for j = 1:N+1
%         tVal = (j-1)*dt;
%         matval(M+1,j) = zmax - exp(-r*(T - tVal));
%     end
% 
%     %-------------------------
%     % 3) Crank–Nicolson Coeffs
%     %-------------------------
%     %  The PDE:  u_t = r u - (1/2) sigma^2 z^2 u_{zz} - r z u_z
%     %  or rearranged into CN form.  
%     %  Typically we define the standard "tri-di" scheme with:
%     %    alpha_i = 0.25*dt * [sigma^2 * i^2  - r i]
%     %    beta_i  = -0.5*dt * [sigma^2 * i^2 + r]
%     %    gamma_i = 0.25*dt * [sigma^2 * i^2 + r i]
%     %  Here, i ~ z-index. We'll do i=0..M, but the PDE is for 1 <= i <= M-1 interior.
% 
%     iIdx = (0:M)';   % integer indices
%     alpha = 0.25 * dt * ( sigma^2 * (iIdx.^2) - r * iIdx );
%     beta  = -0.5 * dt * ( sigma^2 * (iIdx.^2) + r );
%     gamma = 0.25 * dt * ( sigma^2 * (iIdx.^2) + r * iIdx );
% 
%     % Build the tri-di matrix for i=2..M
%     mainDiag = (1 - beta(2:M));
%     lowDiag  = -alpha(3:M);    % subdiagonal
%     upDiag   = -gamma(2:M-1);  % superdiagonal
% 
%     M1 = diag(mainDiag) + diag(lowDiag, -1) + diag(upDiag, 1);
%     [L,U] = lu(M1);
% 
%     mainDiag2 = (1 + beta(2:M));
%     lowDiag2  = alpha(3:M);
%     upDiag2   = gamma(2:M-1);
% 
%     M2 = diag(mainDiag2) + diag(lowDiag2, -1) + diag(upDiag2, 1);
% 
%     %-------------------------
%     % 4) Backward time loop
%     %-------------------------
%     for j = N:-1:1
%         % Known f(i) at time j+1 => matval(i,j+1), find matval(i,j)
%         % Enforce boundary at z=1 => matval(1,j)=0
%         matval(1,j) = 0;
% 
%         % Enforce boundary at z=zmax
%         matval(M+1,j) = zmax - exp(-r*(T - (j-1)*dt));
% 
%         % Build the RHS = M2* f(2..M, j+1) + boundary adjustments
%         rhs = M2 * matval(2:M,j+1);
% 
%         % Add inhomogeneous terms from boundaries, i=1 and i=M+1
%         % subdiagonal => depends on alpha(2), super => depends on gamma(M)
%         % 
%         % * i=2 row in M2 sees i=1 => alpha(2)? 
%         rhs(1)   = rhs(1)   + alpha(2)* matval(1,j)   + alpha(2)* matval(1,j+1);
%         % * i=M row in M2 sees i=M+1 => gamma(M+1)?
%         rhs(end) = rhs(end) + gamma(M)* matval(M+1,j) + gamma(M)* matval(M+1,j+1);
% 
%         % Solve for f(2..M, j)
%         matval(2:M,j) = U \ (L \ rhs);
%     end
% 
%     %-------------------------
%     % 5) Extract price at S0
%     %-------------------------
%     %  S0 => z0 = S0/m
%     z0 = S0 / m;
% 
%     % For safety, clamp z0 to [zVec(1), zVec(end)] if outside
%     if z0 < zVec(1),  z0=zVec(1);  end
%     if z0 > zVec(end),z0=zVec(end);end
% 
%     % Interpolate in z
%     valAt0 = interp1(zVec, matval(:,1), z0, 'linear');
% 
%     % That is the PDE solution u( z0, t=0 ).  For a floating-strike call,
%     % the actual option payoff is  m * u(z).  But typically if we've
%     % non-dimensionalized, check if you need to multiply by something.
%     % If u(t,z) is truly the “dimensionless” payoff, final option value is
%     %   price = m * valAt0.
%     price = m * valAt0;
% 
%     zGrid = zVec;  % return the z-grid
% end
