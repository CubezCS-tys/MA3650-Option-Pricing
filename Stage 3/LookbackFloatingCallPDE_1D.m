function [uGrid, Xvec] = LookbackFloatingCallPDE_1D(r, sigma, T, Xmax, NX, NT)
    % LookbackFloatingCallPDE_1D: Illustrative 1D PDE for a floating-strike
    % lookback call under the dimension-reduction X = S/M.
    %
    % PDE:  du/dt + A(X)*du/dX + B(X)*d2u/dX^2 - r*u = 0 (example form)
    % with boundary conditions:
    %   u(T, X) = X - 1         (since V(T,S,M)/M = (S-M)/M = X - 1)
    %   u(t, X=1) = 0
    %   u(t, Xmax) ~ (Xmax - 1) or something approximate
    %
    % r, sigma = risk-free rate & volatility
    % T        = maturity
    % Xmax     = maximum X
    % NX       = number of X points
    % NT       = number of time steps
    %
    % Returns:
    %   uGrid   = 2D array [NX, NT+1], storing u(t_i, X_j)
    %   Xvec    = spatial grid in X

    % 1) Define X grid
    Xvec = linspace(1, Xmax, NX);
    dX   = Xvec(2) - Xvec(1);

    % 2) Define time grid
    dt = T / NT;
    tvec = linspace(0, T, NT+1);

    % 3) Allocate solution: uGrid(:, n) ~ u at time tvec(n), for X= Xvec
    uGrid = zeros(NX, NT+1);

    % 4) Terminal condition at t=T => n = NT+1 in 1-based indexing
    %    u(T, X) = X - 1
    for i=1:NX
        uGrid(i, end) = Xvec(i) - 1;
    end

    % 5) Main time-stepping loop, going backward from T -> 0
    for n = NT:-1:1
        % previous solution in time (actually the "known" future time in PDE sense)
        uOld = uGrid(:, n+1);

        % let's do an explicit Euler for demonstration
        % for each interior X (i=2..NX-1)
        uNew = uOld;  % copy to start

        for i = 2:NX-1
            Xval = Xvec(i);

            % PDE coefficients (example only) 
            % For dimension-reduced lookback, one might have something like:
            %   dU/dt = 1/2 sigma^2 * [some function of X] * d2U/dX^2
            %         + r* [some function of X]* dU/dX
            %         - r U
            % Insert your actual PDE terms here:
            
            % Dummy example: 
            %   dU/dt = 0.5 * sigma^2 * Xval^2 * d2U/dX^2
            %          + r * Xval * dU/dX
            %          - r * u
            % The actual PDE for dimension-reduced lookback is more involved,
            % but let's keep it short.
            
            % Finite difference approximations:
            uip = uOld(i+1);  % u at i+1
            uim = uOld(i-1);  % u at i-1
            uic = uOld(i);
            
            d2udX2 = (uip - 2*uic + uim) / (dX^2);
            dudX   = (uip - uim) / (2*dX);

            alpha = 0.5*sigma^2*(Xval^2);
            beta  = r*Xval;
            
            % PDE: dUdt = alpha*d2U/dX^2 + beta*dU/dX - r*u
            dUdt = alpha*d2udX2 + beta*dudX - r*uic;

            % Explicit Euler step:
            uNew(i) = uic + dt*dUdt;
        end

        % 6) Enforce boundary conditions at X=1 and X=Xmax
        %   - at X=1 => u=0
        uNew(1) = 0;

        %   - at X=Xmax => near expiry, u(T, X)= X-1 => approximate by
        %     u(t, Xmax)= (Xmax - 1)* e^{-r*(T - t)}, or something simpler
        %     For short time steps we might just do:
        uNew(end) = Xvec(end) - 1;  % or more refined boundary: (X-1)* e^{-r*(T - t_n)}

        % store the new time level
        uGrid(:, n) = uNew;
    end
end
