function C0 = LookbackFloatingCallAnalytic(S0, M, r, T, sigma)
% LOOKBACKFLOATINGCALLANALYTIC 
%  Computes the closed-form price of a floating-strike lookback call 
%  in the Black–Scholes model (no dividends, continuous monitoring).
%
%  PAYOFF at T:   S(T) - min_{0..T} S(t)
%
%  SYNTAX:
%    C0 = LookbackFloatingCallAnalytic(S0, M, r, T, sigma)
%
%  INPUT:
%    S0    - Current underlying price
%    M     - Known running minimum so far (<= S0). 
%            If M == S0, this is a full lookback from inception.
%    r     - Risk-free interest rate (annualized)
%    T     - Time to maturity (in years)
%    sigma - Volatility (annualized)
%
%  OUTPUT:
%    C0    - The fair value of the floating-strike lookback call
%
%  NOTE:
%    Requires the Statistics Toolbox for normcdf (or define your own).

    % Check that M <= S0
    if M > S0
       warning('M > S0 is unusual for a floating-strike lookback call!');
    end

    % Precompute constants
    lambda = 2*r / sigma^2;

    % x1, x2 for the "up and in" part
    x1 = ( log(S0/M) + (r + 0.5*sigma^2)*T ) / (sigma * sqrt(T));
    x2 = x1 - sigma*sqrt(T);

    % y1, y2 for the reflection part
    y1 = ( log(S0/M) + (-r + 0.5*sigma^2)*T ) / (sigma * sqrt(T));
    y2 = y1 - sigma*sqrt(T);

    % Standard normal CDF (from Statistics Toolbox)
    Nx1 = normcdf(x1);
    Nx2 = normcdf(x2);
    Ny1 = normcdf(y1);
    Ny2 = normcdf(y2);

    % (S0/M)^(-lambda)
    ratioPow = (S0/M)^(-lambda);

    % Closed-form price
    C0 = S0 * ( Nx1 + ratioPow * Ny1 ) ...
         - M * exp(-r*T) * ( Nx2 + ratioPow * Ny2 );

end
