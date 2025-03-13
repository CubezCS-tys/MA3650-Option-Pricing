function priceMC = LookbackFloatingCallMC(S0, r, sigma, T, Npaths, Nsteps)
    % Monte Carlo simulation for a floating-strike lookback call
    dt = T/Nsteps;
    discount = exp(-r*T);

    payoffs = zeros(Npaths,1);

    for p = 1:Npaths
        S = S0;
        minS = S0;
        for n = 1:Nsteps
            % dS = r*S dt + sigma*S dW
            dW = sqrt(dt)*randn();
            S = S + r*S*dt + sigma*S*dW;
            if S < 0
                S = 0; % avoid negative, or you could allow a small cutoff
            end
            if S < minS
                minS = S;
            end
        end

        % payoff = S_T - min_{0..T} S_t
        payoff = S - minS;
        payoffs(p) = payoff;
    end

    priceMC = discount * mean(payoffs);
end
