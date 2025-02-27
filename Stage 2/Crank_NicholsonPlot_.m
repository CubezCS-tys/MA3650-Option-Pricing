% Crank-Nicholson.m
clear;
clc;

S0=190;
K=100;
ri=0.05;
dS = 0.1;
dt = 0.005;
T = 0.1 ;
sigma=0.2;
Smax=2000;
r = dt / (dS)^2

M = round(Smax/dS);
dS = Smax/M;
N = round(T/dt);
dt = T/N;
sol = zeros(M+1,N+1);
vetS = linspace(0,Smax,M+1)';
veti = 0:M;
vetj = 0:N;
theta = 0.5


matval(:,N+1) = max(vetS-K,0);
matval(1,:) = 0;
matval(M+1,:) = Smax - K*exp(-ri*dt*(N-vetj));

F = -diag(ones(N-2,1),-1) + 2*diag(ones(N-1,1)) - diag(ones(N-2,1),1)
C= eye(N-1)+(r*theta/2)*F
D= eye(N-1)+(r*(theta-1))*F

for j=1:M
   sol(2:N,j+1) = C \ (D*sol(2:N,j));
end

% Interpolate the option price at a specific stock price (e.g., S0 = 100)

price = interp1(vetS, sol(:, 1), S0);

% Display the price
fprintf('Option price at S0 = %.2f: %.4f\n', S0, price);

% Plot the option price as a function of stock price
figure;
plot(vetS, sol(:, 1), 'b', 'LineWidth', 1.5);
xlabel('Stock Price (S)');
ylabel('Option Price (V)');
title('Put Option Price using Crank-Nicolson');
grid on;
      