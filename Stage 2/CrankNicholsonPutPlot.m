

S0=50; % Initial Stock Price
K= 800; % Strike Price
r=0.05; % Risk-free interest rate
T=0.5; % Time to expiration (6 months)
sigma=0.2; % Volatility
Smax=1500; % Maxium Stock Price
dS=0.5; % Stock Step 
dt=0.001; % Time Step
p=1.5; % Exponxent



[price, vetS, matval] = CrankNicholsonPut(S0,K,r,T,sigma,Smax,dS,dt,p);

%Plot the option price surface
figure;
surf(linspace(0, T, size(matval,2)), vetS, matval);
title('Power Put Option Value Surface');
xlabel('Time to Maturity');
ylabel('Stock Price');
zlabel('Option Value');
shading interp;

% Plot the option value at t = 0
figure;
hold on;
plot(vetS, matval(:,1), 'LineWidth', 2);
title('Power Put CN Option Value');
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
%Display the computed option price
fprintf('The price of the power put option at S0 = %.2f is %.4f\n', S0, price);
