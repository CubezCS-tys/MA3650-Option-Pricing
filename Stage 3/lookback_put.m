% Lookback Put Option Pricing
function price = lookback_put(S, M, r, v, T)
    a1 = a_1(S, M, r, v, T);
    a2 = a_2(S, M, r, v, T);
    a3 = a_3(S, M, r, v, T);
    
    term1 = -S * norm_cdf(-a1);
    term2 = M * exp(-r * T) * norm_cdf(-a2);
    mult = (S * v^2) / (2 * r);
    term3 = norm_cdf(a1) - exp(-r * T) * (M / S)^((2 * r) / v^2) * norm_cdf(a3);
    
    price = term1 + term2 + mult * term3;
end
% Approximation to the cumulative normal distribution function
function cdf = norm_cdf(x)
    k = 1.0 / (1.0 + 0.2316419 * abs(x));
    k_sum = k * (0.319381530 + k * (-0.356563782 + k * (1.781477937 + k * (-1.821255978 + 1.330274429 * k))));
    cdf = 1 - (1 / sqrt(2 * pi)) * exp(-0.5 * x^2) * k_sum;
    
    if x < 0
        cdf = 1 - cdf;
    end
end

% Functions a_1, a_2, a_3
function result = a_1(S, H, r, v, T)
    result = (log(S / H) + (r + 0.5 * v^2) * T) / (v * sqrt(T));
end

function result = a_2(S, H, r, v, T)
    result = a_1(S, H, r, v, T) - v * sqrt(T);
end

function result = a_3(S, H, r, v, T)
    result = a_1(S, H, r, v, T) - (2 * r * sqrt(T) / v);
end


