%function price = EuPutImpl(S0,K,r,T,sigma,Smax,dS,dt)
S0=50;
K=50;
r=0.05;
T=1;
sigma=0.2;
Smax=200;
dS=0.5;
dt=0.1;

%PUT
% % set up grid and adjust increments if necessary
% M = round(Smax/dS);
% dS = Smax/M;
% N = round(T/dt);
% dt = T/N;
% matval = zeros(M+1,N+1);
% vetS = linspace(0,Smax,M+1)';
% veti = 0:M;
% vetj = 0:N;
% % set up boundary conditions
% matval(:,N+1) = max(K-vetS,0);
% matval(1,:) = K*exp(-r*dt*(N-vetj));
% matval(M+1,:) = 0;
% % set up the tridiagonal coefficients matrix
% a = 0.5*(r*dt*veti-sigma^2*dt*(veti.^2));
% b = 1+sigma^2*dt*(veti.^2)+r*dt;
% c = -0.5*(r*dt*veti+sigma^2*dt*(veti.^2));
% coeff = diag(a(3:M),-1) + diag(b(2:M)) + diag(c(2:M-1),1);
% [L,U] = lu(coeff);
% % solve the sequence of linear systems
% aux = zeros(M-1,1);
% for j=N:-1:1
%    aux(1) = - a(2) * matval(1,j); % other term from BC is zero
%    matval(2:M,j) = U \ (L \ (matval(2:M,j+1) + aux));
% end
% % return price, possibly by linear interpolation outside the grid
% price = interp1(vetS, matval(:,1), S0)
% % plot(vetS,matval(:,1))



%CALL
M = round(Smax/dS);
dS = Smax/M;
N = round(T/dt);
dt = T/N;
matval = zeros(M+1,N+1);
vetS = linspace(0,Smax,M+1)';
veti = 0:M;
vetj = 0:N;
% set up boundary conditions
matval(:,N+1) = max(vetS.^p-K,0);
matval(1,:) = 0 ;
matval(M+1,:) = Smax-K*exp(-r * (T-dt));
% set up the tridiagonal coefficients matrix
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
price = interp1(vetS, matval(:,1), S0)
