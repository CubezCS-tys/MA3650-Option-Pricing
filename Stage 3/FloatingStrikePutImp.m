function [price, vetz, matval] = FloatingStrikePutImp(S0,Max,r,T,sigma,Smax,dS,dt)
% set up grid and adjust increments if necessary
Sz = S0/Max;
zmin = Smax/Max;
M = round(zmin/dS);
dz = zmin/M;
N = round(T/dt);
dt = T/N;
matval = zeros(M+1,N+1);
vetS = linspace(0,Smax,M+1)'; % Create a column vector (vetS) containing M+1 equally spaced asset prices from 0 to Smax.
vetz = vetS./Max;
veti = 0:M;
vetj = 0:N;
% set up boundary conditions
matval(:,N+1) = max(1-vetz,0);
matval(1,:) = exp(-r *(T-dt*vetj));
matval(M+1, :) = matval(M, :) / (1 - dz);
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

price = Max * interp1(vetz, matval(:,1), Sz);