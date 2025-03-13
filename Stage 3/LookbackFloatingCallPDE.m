function [Vgrid, Svec, Mvec] = LookbackFloatingCallPDE(r, sigma, T, Smax, NS, NM, NT)
    % LookbackFloatingCallPDE
    % Solve for a floating-strike lookback call option using a 2D PDE in (S,M).
    %
    % PDE: dV/dt + 1/2 sigma^2 S^2 d^2V/dS^2 + r S dV/dS - rV = 0 (for S > M)
    % Subject to:
    %   V(T, S, M) = S - M
    %   V(t, S=M, M) = 0
    %   and some far-field boundary condition at S = Smax
    %
    % r      = risk-free rate
    % sigma  = volatility
    % T      = time to maturity
    % Smax   = maximum S in the grid
    % NS     = number of S grid points
    % NM     = number of M grid points
    % NT     = number of time steps
    %
    % Returns:
    %   Vgrid - 2D array of option values at t=0 (Vgrid(i,j) ~ V(0, S_i, M_j))
    %   Svec, Mvec - The mesh grids for the S and M coordinates

    % 1) Setup spatial grids
    Svec = linspace(0, Smax, NS);  % asset price
    Mvec = linspace(0, Smax, NM);  % running minimum
    dS = Svec(2) - Svec(1);
    dM = Mvec(2) - Mvec(1);

    % 2) Setup time grid
    dt = T / NT;
    tvec = linspace(0, T, NT+1);
    % We'll step backward from T to 0
    % but store the solution primarily at t=0

    % 3) Initialize solution V(t, S, M). We store entire time for demonstration,
    %    but you can also store only the final slice if memory is an issue.
    V = zeros(NS, NM, NT+1);

    % 4) Terminal condition: at t = T, V(T,S,M) = S - M (but only for S >= M)
    for i = 1:NS
        for j = 1:NM
            if Svec(i) >= Mvec(j)
                V(i,j, end) = Svec(i) - Mvec(j);
            else
                V(i,j, end) = 0;  % infeasible region or no payoff
            end
        end
    end

    % 5) Main PDE loop, stepping backward in time:
    for n = NT:-1:1   % from NT down to 1
        t_now = tvec(n+1);
        t_prev = tvec(n);

        % We want to approximate V(:, :, n) from V(:, :, n+1).
        Vn_plus_1 = V(:, :, n+1);  % solution at next time (later in real time)

        % One can apply an ADI or Crank-Nicolson. For brevity, let's do a
        % simplistic explicit/implicit "mixed" step (not fully stable for large dt).
        % In production code, you'd likely use an ADI approach.

        %--- Copy old solution to start iteration
        Vn = Vn_plus_1;

        % We do a few relaxation sweeps to approximate the solution at time step n:
        for iter = 1:2
            % Sweep over the interior points (i=2..NS-1, j=2..NM-1) while skipping
            % S<M region.
            for i = 2:NS-1
                for j = 2:NM-1
                    if Svec(i) >= Mvec(j)
                        S_ij = Svec(i);
                        % Coefficients for PDE:
                        alpha = 0.5 * sigma^2 * S_ij^2;
                        beta  = r * S_ij;
                        % PDE is: dV/dt = alpha * d2V/dS^2 + beta * dV/dS - r V

                        % Finite differences in S-direction (second-order central in S):
                        V_Splus  = Vn(i+1,j);
                        V_Sminus = Vn(i-1,j);
                        V_S0     = Vn(i,j);
                        d2VdS2   = (V_Splus - 2*V_S0 + V_Sminus)/(dS^2);
                        dVdS     = (V_Splus - V_Sminus)/(2*dS);

                        % Explicit Euler update (very rough):
                        dV_dt = alpha * d2VdS2 + beta * dVdS - r * V_S0;
                        Vn(i,j) = V_S0 + dt * dV_dt;
                    else
                        % S < M region or out of domain
                        Vn(i,j) = 0;
                    end
                end
            end

            % Enforce boundary conditions:
            %  a) S=M => V=0
            for idx = 1:min(NS,NM)
                if abs(Svec(idx) - Mvec(idx)) < 1e-12
                    Vn(idx, idx) = 0;
                end
            end
            %  b) S=0 => V=0
            Vn(1, :, :) = 0;
            %  c) M=0 => running minimum = 0 => payoff could be large, but must handle carefully
            %     Here let's simply keep the PDE solution or impose that M=0 => payoff S - 0 = S if at maturity
            %     For t<T, no immediate boundary requirement except continuity. We'll just keep the PDE solution.
            %  d) S=Smax => far field approx => S>>M => V ~ S - M e^{-r(T-t)}
            Vn(end, :) = Smax - Mvec(:) * exp(-r*(T - t_prev));  % rough boundary guess

        end

        % Assign new solution
        V(:,:,n) = Vn;
    end

    % 6) The PDE solution at t=0 is V(:,:,1)
    Vgrid = V(:,:,1);
end
