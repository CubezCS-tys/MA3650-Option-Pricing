%% Black-Scholes Power Call Option Function
function price = BlackScholesPowerCall(S, K, T, t, r, sigma, p)
    tau = T - t; % Time to maturity
    h1 = (log(S / K^(1/p)) + (r + (p - 0.5) * sigma^2) * tau) / (sigma * sqrt(tau));
    h2 = h1 - p * sigma * sqrt(tau);
    
    % Calculate the call price
    price = (S^p) * exp((p - 1) * (r + 0.5 * p * sigma^2) * tau) * normcdf(h1) - K * exp(-r * tau) * normcdf(h2);
end
