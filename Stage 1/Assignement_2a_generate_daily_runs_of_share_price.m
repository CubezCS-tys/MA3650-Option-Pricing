% S0 = 50;                   % Initial asset price
% mu = 0.05;                 % Drift
% sigma_vals = [0.01, 0.05, 0.1, 0.2, 0.5]; % Volatilities
% T_vals = [0.25, 0.5];      % Time horizons (3 and 6 months)
% NSteps_3m = 90;            % Number of time steps for 3 months
% NSteps_6m = 180;           % Number of time steps for 6 months
% NRepl = 1;                 % Number of paths per time horizon (total 10 paths)
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


S0 = 50;                   % Initial asset price
mu = 0.05;                 % Drift
sigma_vals = [0.01, 0.05, 0.1, 0.2, 0.5]; % Volatilities
T_vals = [0.25, 0.5];      % Time horizons (3 and 6 months)
NSteps_3m = 90;            % Number of time steps for 3 months
NSteps_6m = 180;           % Number of time steps for 6 months
NRepl = 1;                 % Number of paths per time horizon (total 10 paths)

% Prepare figure for subplots
figure;
% Loop through each time horizon
for t = 1:length(T_vals)
    T = T_vals(t);  % Current time horizon

    % Create subplot for each time horizon
    subplot(1, 2, t);  % Create a 1x2 grid of subplots
    hold on;

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

        % Plot each path with a label for the time horizon
        plot(1:NSteps+1, S, 'LineWidth', 1.5, 'DisplayName', ['\sigma = ', num2str(sigma)]);
    end

    % Label the subplot
    xlabel('Time Steps');
    ylabel('Asset Price');
    if T == 0.5
        title(['Paths for 6 months']);
    else
        title(['Paths for 3 months']);
    end
    legend show;
    grid on;
    hold off;
end

% Adjust the layout to ensure subplots are clear
sgtitle('Simulated Asset Price Paths with Different Volatilities and Time Horizons');


