function price = powerOptionPrice(S0, K, q, r, sigma, T, optionType)
% powerOptionPrice Prices a European power option under Black-Scholes
% without requiring the Statistics Toolbox.
%
% INPUTS:
%   S0        - Current underlying price
%   K         - Strike price (for the power payoff)
%   q         - Power exponent
%   r         - Risk-free interest rate (annualized)
%   sigma     - Volatility of underlying (annualized)
%   T         - Time to maturity in years
%   optionType- 'call' or 'put'
%
% OUTPUT:
%   price - The fair value of the power option at time 0.

    if nargin < 7
        optionType = 'call';
    end

    % Compute d2 and d1 as per formula
    d2 = (log(S0) - (1/q)*log(K) + (r - 0.5*sigma^2)*T) / (sigma*sqrt(T));
    d1 = d2 + q*sigma*sqrt(T);

    % Discount factor and intermediate term A
    discountFactor = exp(-r*T);
    A = S0^(q) * exp((q*(q-1)*sigma^2/2 + q*r)*T);

    switch lower(optionType)
        case 'call'
            price = discountFactor * ( A * myNormCdf(d1) - K * myNormCdf(d2) );
        case 'put'
            price = discountFactor * ( K * myNormCdf(-d2) - A * myNormCdf(-d1) );
        otherwise
            error('optionType must be either ''call'' or ''put''.');
    end
end

function p = myNormCdf(x)
    %MYNORMCDF Compute the standard normal CDF using the error function.
    p = 0.5*(1 + erf(x./sqrt(2)));
end
