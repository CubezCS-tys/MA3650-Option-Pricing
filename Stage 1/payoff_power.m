% Parameters for the power option
K = 100;       % Strike price
p = 2;         % Power exponent
S = 0:0.1:20; % Range of underlying prices with finer increments

% Calculate the payoffs
power_call_payoff = max(S.^p - K, 0);
power_put_payoff = max(K - S.^p, 0);

% Critical underlying asset price where payoff starts to become positive
critical_S = K^(1/p);

% Plot Power Call Option Payoff
figure;
plot(S, power_call_payoff, 'LineWidth', 2);
hold on;
xline(critical_S, '--', 'Color', 'r', 'LineWidth', 1.5); % Vertical line at S = K^{1/p}
title(['Power Call Option Payoff (K = ' num2str(K) ', p = ' num2str(p) ')']);
xlabel('Underlying Asset Price (S)');
ylabel('Payoff');
legend('Payoff', ['S = K^{1/p} = ' num2str(critical_S)]);
grid on;
hold off;

% Plot Power Put Option Payoff
figure;
plot(S, power_put_payoff, 'LineWidth', 2);
hold on;
xline(critical_S, '--', 'Color', 'r', 'LineWidth', 1.5); % Vertical line at S = K^{1/p}
title(['Power Put Option Payoff (K = ' num2str(K) ', p = ' num2str(p) ')']);
xlabel('Underlying Asset Price (S)');
ylabel('Payoff');
legend('Payoff', ['S = K^{1/p} = ' num2str(critical_S)]);
grid on;
hold off;
