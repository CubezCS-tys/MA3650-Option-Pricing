% Parameters
S0 = 50;         % Initial asset price
mu = 0.05;       % Drift
sigma = 0.2;     % Volatility
T = 1;           % Time to expiration (1 year)
NSteps = 252;    % Number of time steps (daily steps)
NRepl = 1;       % Number of paths

% Generate a sample asset price path
S = AssetPaths(S0, mu, sigma, T, NSteps, NRepl);

% Floating Strike Lookback Call Option Payoff
S_min = min(S);     % Minimum price observed in the path
S_T = S(end);       % Price at expiration
call_payoff = max(0, S_T - S_min);

% Floating Strike Lookback Put Option Payoff
S_max = max(S);     % Maximum price observed in the path
put_payoff = max(0, S_max - S_T);

% Plot the asset price path
figure;
plot(1:length(S), S, 'LineWidth', 2);
hold on;
yline(S_min, '--', 'Color', 'b', 'LineWidth', 1.5); % Minimum line for Call
yline(S_max, '--', 'Color', 'r', 'LineWidth', 1.5); % Maximum line for Put
title('Asset Price Path with Min and Max Prices');
xlabel('Time Steps');
ylabel('Asset Price');
legend('Asset Path', 'Min Price (Call Strike)', 'Max Price (Put Strike)');
grid on;
hold off;

% Display Payoffs
fprintf('Floating Strike Lookback Call Payoff: %.2f\n', call_payoff);
fprintf('Floating Strike Lookback Put Payoff: %.2f\n', put_payoff);
