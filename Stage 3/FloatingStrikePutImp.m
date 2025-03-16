function [price, vetS, matval] = FloatingStrikePutImp(S0,Max,r,T,sigma,Smax,dS,dt)
%set up grid and adjust increments if necessary
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
% matval(M+1, :) = matval(M, :) / (1 - dz);
matval(M+1, :) = 0;
a = 0.5*(r*dt*veti-sigma^2*dt*(veti.^2));
b = 1+sigma^2*dt*(veti.^2)+r*dt;
c = -0.5*(r*dt*veti+sigma^2*dt*(veti.^2));
coeff = diag(a(3:M),-1) + diag(b(2:M)) + diag(c(2:M-1),1);
[L,U] = lu(coeff);
% solve the sequence of linear systems
aux = zeros(M-1,1);
for j=N:-1:1
   aux(1) = a(2) * matval(M+1,j);
   matval(2:M,j) = U \ (L \ (matval(2:M,j+1) + aux));
end

price = Max * interp1(vetz, matval(:,1), Sz);


