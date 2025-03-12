function priceMC = MonteCarloPartialLookback(S0, fixedMin, r, sigma, T, numPaths, numSteps)
    % MonteCarloPartialLookback computes the partial lookback call price 
    % using Monte Carlo simulation. The payoff is max(S(T) - fixedMin, 0)
    % where fixedMin is provided.
    
    dt = T / numSteps;
    S = zeros(numSteps+1, numPaths);
    S(1,:) = S0;
    dW = sqrt(dt) * randn(numSteps, numPaths);
    
    for i = 2:numSteps+1
        S(i,:) = S(i-1,:) .* exp((r - 0.5*sigma^2)*dt + sigma*dW(i-1,:));
    end
    
    S_T = S(end,:);
    payoff = max(S_T - fixedMin, 0);
    priceMC = exp(-r*T) * mean(payoff);
end
