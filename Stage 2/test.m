S0=85;
K=100;
r=0.05;
T=1;
sigma=0.2;
Smax=2*max(S0,K)*exp(r*T);
%dS=1;
%dt=0.1;
dS = 0.5;
dt = 0.01;
Sb = 80
M = round((Smax-Sb)/dS);
dS = (Smax-Sb)/M;
% M = round(Smax/dS);
% dS = Smax/M;
N = round(T/dt);
dt = T/N;
matval = zeros(M+1,N+1);
vetS = linspace(Sb,Smax,M+1)'
veti = vetS/dS;
vetj = 0:N;
% set up boundary conditions
matval(:,N+1) = max(K-vetS,0);
matval(1,:) = 0 ;
matval(M+1,:) = Smax-K*exp(-r * (T-dt));
%matval(M+1,:) = 0;


% M = round(Smax/dS);
% dS = Smax/M;
% N = round(T/dt);
% dt = T/N;
% matval = zeros(M+1,N+1);
% vetS = linspace(0,Smax,M+1)';
% veti = 0:M;
% vetj = 0:N;
% % set up boundary conditions
% matval(:,N+1) = max(vetS-K,0);
% matval(1,:) = 0 ;
% %matval(M+1,:) = Smax-K*exp(-r * (T-dt));
% matval(M+1,:) = (Smax-K)*exp(-r * dt * (N-vetj));
% % set up the coefficients matrix
% alpha = 0.25*dt*(sigma^2*(veti.^2) - r*veti);
% beta = -dt*0.5*(sigma^2*(veti.^2) + r);
% gamma = 0.25*dt*(sigma^2*(veti.^2) + r*veti);
% M1 = -diag(alpha(3:M) ,-1) + diag(1-beta(2:M)) - diag(gamma(2:M-1) ,1) ;
% [L,U] = lu(M1);
% M2 = diag(alpha(3:M) ,-1) + diag(1+beta(2:M)) + diag(gamma(2:M-1) ,1);
% %solve the sequence of linear systems
% lostval = zeros(size(M2,2), 1);
% for j=N:-1:1
%     if length(lostval)>1
%         lostval(1) = alpha(2) * (matval(1,j)+matval(1,j+1));
%         lostval(end) = gamma(end) * (matval(end,j)+matval(end,j+1));
%     else
%         lostval = lostval(1)+lostval(end);
%     end
% 
%     matval(2:M,j) = U \ (L \ (M2*(matval(2:M,j+1)+lostval)));
% end
% 
% %return price, possibly by linear interpolation outside the grid
% price = interp1(vetS, matval(:,1), S0)
% 
% 
% 
% 
% % M = round(Smax/dS);
% % dS = Smax/M;   
% % N = round(T/dt);
% % dt = T/N;
% % matval = zeros(M+1,N+1);
% % vetS = linspace(0,Smax,M+1)';
% % veti = 0:M;
% % vetj = 0:N;
% % % set up boundary conditions
% % matval(:,N+1) = max(vetS-K,0);
% % matval(1,:) = 0;
% % matval(M+1,:) = Smax - K*exp(-r*dt*(N-vetj));
% % % set up the tridiagonal coefficients matrix
% % a = 0.5*(r*dt*veti-sigma^2*dt*(veti.^2));
% % b = 1+sigma^2*dt*(veti.^2)+r*dt;
% % c = -0.5*(r*dt*veti+sigma^2*dt*(veti.^2));
% % coeff = diag(a(3:M),-1) + diag(b(2:M)) + diag(c(2:M-1),1);
% % [L,U] = lu(coeff);
% % % solve the sequence of linear systems
% % aux = zeros(M-1,1);
% % for j=N:-1:1
% %    aux(1) = - a(2) * matval(1,j); % other term from BC is zero
% %    matval(2:M,j) = U \ (L \ (matval(2:M,j+1) + aux));
% % end
% % 
% % % return price, possibly by linear interpolation outside the grid
% % price = interp1(vetS, matval(:,1), S0)
% 
% 
% 
% % 
% % % Parameters
% % S0 = 150;    % Current underlying price
% % K = 140;     % Strike on the power payoff
% % q = 1.4;       % Power exponent
% % r = 0.05;    % Risk-free rate
% % sigma = 0.8; % Volatility
% % T = 1;       % Time to maturity in years
% % 
% % % Price a power call
% % callPrice = powerOptionPrice(S0, K, q, r, sigma, T, 'call');
% % disp(['Power Call Price: ', num2str(callPrice)]);
% % 
% % % Price a power put
% % putPrice = powerOptionPrice(S0, K, q, r, sigma, T, 'put');
% % disp(['Power Put Price: ', num2str(putPrice)]);
% 
% % M = round(Smax/dS);
% % dS = Smax/M;
% % N = round(T/dt);
% % dt = T/N;
% % matval = zeros(M+1,N+1);
% % vetS = linspace(0,Smax,M+1)';
% % veti = 0:M;
% % vetj = 0:N;
% % % set up boundary conditions
% % matval(:,N+1) = max(vetS.^p-K,0);
% % matval(1,:) = 0 ;
% % matval(M+1,:) = 0;
% % % set up the tridiagonal coefficients matrix
% % a = 0.5*(r*dt*veti-sigma^2*dt*(veti.^2));
% % b = 1+sigma^2*dt*(veti.^2)+r*dt;
% % c = -0.5*(r*dt*veti+sigma^2*dt*(veti.^2));
% % coeff = diag(a(3:M),-1) + diag(b(2:M)) + diag(c(2:M-1),1);
% % [L,U] = lu(coeff);
% % % solve the sequence of linear systems
% % aux = zeros(M-1,1);
% % for j=N:-1:1
% %    aux(1) = - a(2) * matval(1,j); % other term from BC is zero
% %    matval(2:M,j) = U \ (L \ (matval(2:M,j+1) + aux));
% % end
% % % return price, possibly by linear interpolation outside the grid
% % price = interp1(vetS, matval(:,1), S0);
% % 
% % matt = transpose(matval);
% % mattf = flipud(transpose(matval));
