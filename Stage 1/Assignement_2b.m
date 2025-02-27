%

S0 = 50;                   % Initial asset price
mu = 0.05;                 % Drift
sigma_vals = [0.01, 0.05, 0.1, 0.2, 0.5]; % Volatilities
T_vals = [0.25, 0.5];      % Time horizons (3 and 6 months)
NSteps_3m = 90;            % Number of time steps for 3 months
NSteps_6m = 180;           % Number of time steps for 6 months
NRepl = 1;                 % Number of paths per time horizon (total 10 paths)
K = 75;
p =1.1;
% Initialize arrays to store results
volatility = [];
time_horizon = [];
power_call_payoffs = [];
power_put_payoffs = [];
lookback_call_payoffs = [];
lookback_put_payoffs = [];


% Loop through each time horizon
for t = 1:length(T_vals)
    T = T_vals(t);  % Current time horizon

    for s = 1:length(sigma_vals)
        sigma = sigma_vals(s);  % Current volatility

        % Set NSteps based on time horizon
        if T == 0.25
            NSteps = NSteps_3m;
        else
            NSteps = NSteps_6m;
        end

                % Generate asset price path
        SPaths = AssetPaths(S0, mu, sigma, T, NSteps, NRepl);
        S = SPaths(1, :);  % Extract the single path
        S_T = S(end); 
        
        % Calculate power option payoffs
        power_call_payoff = max(S_T.^p - K, 0);
        power_put_payoff = max(K - S_T.^p, 0);

        % Calculate floating strike lookback option payoffs
        S_min = min(S);    % Minimum price observed in the path
        call_payoff = max(0, S_T - S_min);  % Lookback call payoff
        S_max = max(S);    % Maximum price observed in the path
        put_payoff = max(0, S_max - S_T);  % Lookback put payoff

        % Store results in arrays
        volatility = [volatility; sigma];
        time_horizon = [time_horizon; T];
        power_call_payoffs = [power_call_payoffs; power_call_payoff];
        power_put_payoffs = [power_put_payoffs; power_put_payoff];
        lookback_call_payoffs = [lookback_call_payoffs; call_payoff];
        lookback_put_payoffs = [lookback_put_payoffs; put_payoff];

    end
end

% Create table of results
results_table = table(volatility, time_horizon, power_call_payoffs, power_put_payoffs, lookback_call_payoffs, lookback_put_payoffs, ...
                      'VariableNames', {'Volatility', 'Time Horizon', 'Power Call Payoff', 'Power Put Payoff', 'Lookback Call Payoff', 'Lookback Put Payoff'});

% Display the table in the Command Window
disp(results_table);

% Optionally, display the table in a figure
figure;
uitable('Data', results_table{:,:}, 'ColumnName', results_table.Properties.VariableNames, 'Position', [20 100 620 200]);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% % Parameters
% S0 = 50;                   % Initial asset price
% mu = 0.05;                 % Drift
% sigma_vals = [0.01, 0.05, 0.1, 0.2, 0.5]; % Volatilities
% T_vals = [0.25, 0.5];      % Time horizons (3 and 6 months)
% NSteps_3m = 90;            % Number of time steps for 3 months
% NSteps_6m = 180;           % Number of time steps for 6 months
% NRepl = 1;                 % Number of paths per time horizon (total 10 paths)
% K = 75;
% p = 1.1;
% 
% % Initialize arrays to store results
% volatility = [];
% time_horizon = [];
% power_call_payoffs = [];
% power_put_payoffs = [];
% lookback_call_payoffs = [];
% lookback_put_payoffs = [];
% 
% % Prepare figure for subplots
% figure;
% for s = 1:length(sigma_vals)
%     sigma = sigma_vals(s);  % Current volatility
% 
%     % Create subplot for each volatility
%     subplot(3, 2, s);  % Create a 3x2 grid of subplots
%     hold on;
% 
%     % Loop through each time horizon
%     for t = 1:length(T_vals)
%         T = T_vals(t);  % Current time horizon
% 
%         % Set NSteps based on time horizon
%         if T == 0.25
%             NSteps = NSteps_3m;
%         else
%             NSteps = NSteps_6m;
%         end
% 
%         % Generate asset price path
%         SPaths = AssetPaths(S0, mu, sigma, T, NSteps, NRepl);
%         S = SPaths(1, :);  % Extract the single path
%         S_T = S(end); 
%         
%         % Calculate power option payoffs
%         power_call_payoff = max(S_T.^p - K, 0);
%         power_put_payoff = max(K - S_T.^p, 0);
% 
%         % Calculate floating strike lookback option payoffs
%         S_min = min(S);    % Minimum price observed in the path
%         call_payoff = max(0, S_T - S_min);  % Lookback call payoff
%         S_max = max(S);    % Maximum price observed in the path
%         put_payoff = max(0, S_max - S_T);  % Lookback put payoff
% 
%         % Store results in arrays
%         volatility = [volatility; sigma];
%         time_horizon = [time_horizon; T];
%         power_call_payoffs = [power_call_payoffs; power_call_payoff];
%         power_put_payoffs = [power_put_payoffs; power_put_payoff];
%         lookback_call_payoffs = [lookback_call_payoffs; call_payoff];
%         lookback_put_payoffs = [lookback_put_payoffs; put_payoff];
% 
%         % Plot each path with a label for the time horizon
%         plot(1:NSteps+1, S, 'LineWidth', 1.5, 'DisplayName', ['T = ', num2str(T), ' years']);
%     end
% 
%     % Label the subplot
%     xlabel('Time Steps');
%     ylabel('Asset Price');
%     title(['Paths for \sigma = ', num2str(sigma)]);
%     legend show;
%     grid on;
%     hold off;
% end
% 
% % Adjust the layout to ensure subplots are clear
% sgtitle('Simulated Asset Price Paths with Different Volatilities and Time Horizons');
% 
% % Create table of results
% results_table = table(volatility, time_horizon, power_call_payoffs, power_put_payoffs, lookback_call_payoffs, lookback_put_payoffs, ...
%                       'VariableNames', {'Volatility', 'Time_Horizon', 'Power_Call_Payoff', 'Power_Put_Payoff', 'Lookback_Call_Payoff', 'Lookback_Put_Payoff'});
% 
% % Display the table in the Command Window
% disp(results_table);
% 
% % Optionally, display the table in a figure
% figure;
% uitable('Data', results_table{:,:}, 'ColumnName', results_table.Properties.VariableNames, 'Position', [20 100 640 200]);
