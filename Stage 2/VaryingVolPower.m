% Define the parameters common to all runs
S0 = 50;       % Initial Stock Price 
K = 800;        % Strike Price
r = 0.05;       % Risk-free interest rate
q = 0;          % Dividend yield (if any)
Smax = 1000;     % Maximum Stock Price
dS = 0.5;       % Stock Price step 
dt = 0.001;     % Time step 
p = 1.3;          % Exponent for the power option

% Define the sets of volatilities and expiration times to be tested
%vol_values = [0.15, 0.20, 0.25, 0.30, 0.35]; % 5 volatility values
vol_values = [0.1, 0.20, 0.5, 0.8];
T_values   = [0.5, 1.0];                      % 2 expiration times (e.g., 6 months and 1 year)

% Loop over each time to expiration
for j = 1:length(T_values)
    T = T_values(j);
    figure;
    hold on;
    
    % Loop over each volatility value
    for i = 1:length(vol_values)
        sigma = vol_values(i);
        % Run the Power Option Pricer
        [price, vetS, matval] = CrankNicholsonPut(S0, K, r, T, sigma, Smax, dS, dt, p);
        
        % Plot the option value at t = 0 (today)
        plot(vetS, matval(:,1), 'LineWidth', 2, ...
             'DisplayName', sprintf('\\sigma = %.2f', sigma));
         
        % Optionally, display the computed price for each case in the command window
        fprintf('T = %.1f, sigma = %.2f, Price: %.4f\n', T, sigma, price);
    end
    
    % Format the figure for the current time horizon
    title(sprintf('Power Call Option CN Value at t = 0 for T = %.1f', T));
    xlabel('Stock Price');
    ylabel('Option Value');
    legend('show');
    grid on;
    hold off;
end
