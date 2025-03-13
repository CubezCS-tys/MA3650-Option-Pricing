% function [price, vetS, matval] = FloatingStrikeCallImp(S0,Min,r,T,sigma,Smax,dS,dt)
% % set up grid and adjust increments if necessary
% Sz = S0/Min;
% z = Smax/Min;
% M = round((z)/dS);
% dS = (z)/M;
% N = round(T/dt);
% dt = T/N;
% matval = zeros(M+1,N+1);
% vetz = linspace(0,z,M+1)';
% veti = 0:M;
% vetj = 0:N;
% % set up boundary conditions
% matval(:,N+1) = vetz;
% matval(1,:) = 0;
% matval(M+1,:) = max(z-1, 0);
% % set up the tridiagonal coefficients matrix
% a = 0.5*(r*dt*veti-sigma^2*dt*(veti.^2));
% b = 1+sigma^2*dt*(veti.^2)+r*dt;
% c = -0.5*(r*dt*veti+sigma^2*dt*(veti.^2));
% coeff = diag(a(3:M),-1) + diag(b(2:M)) + diag(c(2:M-1),1);
% [L,U] = lu(coeff);
% % solve the sequence of linear systems
% aux = zeros(M-1,1);
% for j=N:-1:1
%    aux(M-1) = - c(M) * matval(M+1,j);
%    matval(2:M,j) = U \ (L \ (matval(2:M,j+1) + aux));
% end
% % return price, possibly by linear interpolation outside the grid
% price = Min*interp1(vetz, matval(:,1), Sz);
% matval;


function [price, vetz, matval] = FloatingStrikeCallImp(S0,Min,r,T,sigma,Smax,dS,dt)
% set up grid and adjust increments if necessary
Sz = S0/Min;
zmax = Smax/Min;
M = round(zmax/dS);
dS = zmax/M;
N = round(T/dt);
dt = T/N;
matval = zeros(M+1,N+1);
vetS = linspace(0,Smax,M+1)'; % Create a column vector (vetS) containing M+1 equally spaced asset prices from 0 to Smax.
vetz = vetS./Min;
veti = 0:M;
vetj = 0:N;
% set up boundary conditions
matval(:,N+1) = max(vetz-1,0);
matval(1,:) = 0;
matval(M+1,:) = (zmax)-exp(-r *(T-dt*vetj));
% set up the coefficients matrix
alpha = 0.25*dt*(sigma^2*(veti.^2) - r*veti);
beta = -dt*0.5*(sigma^2*(veti.^2) + r);
gamma = 0.25*dt*(sigma^2*(veti.^2) + r*veti);
M1 = -diag(alpha(3:M) ,-1) + diag(1-beta(2:M)) - diag(gamma(2:M-1) ,1) ;
[L,U] = lu(M1);
M2 = diag(alpha(3:M) ,-1) + diag(1+beta(2:M)) + diag(gamma(2:M-1) ,1);
%solve the sequence of linear systems
aux = zeros(size(M2,2), 1);
for j=N:-1:1
    if length(aux)>1
        aux(1) = alpha(2) * (matval(1,j)+matval(1,j+1));
        aux(end) = gamma(end) * (matval(end,j)+matval(end,j+1));
    else
        aux = aux(1)+aux(end);
    end

    matval(2:M,j) = U \ (L \ ((M2*(matval(2:M,j+1))+aux)));
end

%return price, possibly by linear interpolation outside the grid

price = Min * interp1(vetz, matval(:,1), Sz);

% function [price, vetz, matval] = FloatingStrikeCallImp(S0, Min, r, T, sigma, Smax, dS, dt)
% % Set up grid and adjust increments if necessary
% Sz = S0 / Min;
% z = Smax / Min;
% M = round(z / dS);
% dS = z / M;  % Adjust dS to fit the grid
% N = round(T / dt);
% dt = T / N;  % Adjust dt to fit the time steps
% matval = zeros(M+1, N+1);
% vetz = linspace(0, z, M+1)';  % Grid points in z-space
% veti = 0:M;  % Indices for z (0 to M)
% vetj = 0:N;
% 
% % Set up boundary conditions
% matval(:, N+1) = max(vetz - 1, 0);  % Payoff at maturity (z_T - 1)
% matval(1, :) = 0;  % Lower boundary (S=0)
% matval(M+1, :) = vetz(end) - exp(-r * dt * (N - vetj));  % Upper boundary (z - exp(-r(T-t)))
% 
% % Compute coefficients with correct dS scaling
% alpha = 0.25 * dt/dz * (sigma^2 * (veti.^2 * dS^2) - r * veti * dS);
% beta = -0.5 * dt/dz * (sigma^2 * (veti.^2 * dS^2) + r);
% gamma = 0.25 * dt/dz * (sigma^2 * (veti.^2 * dS^2) + r * veti * dS);
% 
% % Build tridiagonal matrices M1 (implicit) and M2 (explicit)
% % M1 is for (I - 0.5*dt*A), M2 for (I + 0.5*dt*A)
% M1 = -diag(alpha(3:M), -1) + diag(1 - beta(2:M)) - diag(gamma(2:M-1), 1);
% M2 = diag(alpha(3:M), -1) + diag(1 + beta(2:M)) + diag(gamma(2:M-1), 1);
% 
% % LU decomposition of M1 for efficient solving
% [L, U] = lu(M1);
% 
% % Solve backward in time
% aux = zeros(M-1, 1);  % Vector to store boundary contributions
% for j = N:-1:1
%     % Contributions from boundaries (lower and upper)
%     aux(1) = alpha(2) * (matval(1, j) + matval(1, j+1));
%     aux(end) = gamma(M) * (matval(M+1, j) + matval(M+1, j+1));
%     
%     % Solve M1 * matval(2:M, j) = M2 * matval(2:M, j+1) + aux
%     matval(2:M, j) = U \ (L \ (M2 * matval(2:M, j+1) + aux));
% end
% 
% % Interpolate to find the price at the initial scaled price Sz
% price = Min * interp1(vetz, matval(:, 1), Sz, 'linear', 0);
% end
