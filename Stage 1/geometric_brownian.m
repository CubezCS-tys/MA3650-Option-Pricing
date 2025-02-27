NSteps = 360;
T = 1;
S0 = 50;
mu = 0.01;
sigma = 0.3;
t = (0:1:NSteps)'/NSteps;
W = [0; cumsum(randn(NSteps,1))]/sqrt(NSteps);
t = t*T;
W = W*sqrt(T);
X = (mu-(sigma^2)/2)*t + 1 * W;
Y = S0*exp(X);
plot(t,Y);
title('Sample path of asset prices following the GBM') %with volatility', num2str(sigma), 'and drift', num2str(mu))
xlabel('Time')
ylabel('Asset Value')
