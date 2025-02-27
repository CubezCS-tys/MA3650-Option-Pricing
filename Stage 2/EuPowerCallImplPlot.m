% S0=100;
% K=100;
% r=0.05;
% T=1;
% sigma=0.8;
% Smax=200;
% dS=0.5;
% dt=0.001;
% p = 2;

S0=150;
K=140;
r=0.05;
T=1;
sigma=0.2;
Smax=400;
%Smax = 2*max(S0,K)*exp(r*T);
dS=1;
dt=0.01;
p = 1.5;

% Run the Power Put Pricer
[price, vetS, matval] = EuPowerCallImpl(S0, K, r, T, sigma, Smax, dS, dt, p);

% Plot the option price surface
figure;
surf(linspace(0, T, size(matval,2)), vetS, matval);
title('Power Call Option Value Surface');
xlabel('Time to Maturity');
ylabel('Stock Price');
zlabel('Option Value');
shading interp;

% Plot the option value at t = 0
figure;
hold on;
plot(vetS, matval(:,1), 'LineWidth', 2);
title('Power Call Option Value at t = 0');
xlabel('Stock Price');
ylabel('Option Value');
grid on;

% Plot the option value at t = 0
plot(vetS, matval(:,end), 'LineWidth', 2);
xlabel('Stock Price');
ylabel('Option Value');
grid on;
legend('t = 0 (Today)', 't = T (Expiration)');
hold off;
% Display the computed option price
fprintf('The price of the power call option at S0 = %.2f is %.4f\n', S0, price);
