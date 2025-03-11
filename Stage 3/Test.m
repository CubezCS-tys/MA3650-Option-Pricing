% Parameters
r = 0.1;           % Risk-free rate
sigma = 0.4;        % Volatility
T = 1;              % Time to maturity
x_max = 5;          % Maximum stock price ratio (S/m)
N = 600;            % Number of space steps
M = 300;            % Number of time steps

% Grid setup
dx = (x_max - 1) / N;
dt = T / M;
x = 1:dx:x_max;     % Stock price ratio grid (S/m)
n = length(x);      % Number of grid points

% Terminal condition (H(x,T) = max(x - 1, 0))
H = max(x - 1, 0)';

% Precompute coefficients for the PDE
a = zeros(n, 1);
b = zeros(n, 1);
c = zeros(n, 1);

for i = 1:n
    xi = x(i);
    a(i) = ( (r - 0.5*sigma^2) * xi / (2*dx) ) + (0.5 * sigma^2 * xi^2 / dx^2);
    b(i) = (sigma^2 * xi^2 / dx^2) + r;
    c(i) = ( - (r - 0.5*sigma^2) * xi / (2*dx) ) - (0.5 * sigma^2 * xi^2 / dx^2);
end

% Time-stepping loop
for m = M:-1:1
    t = (m-1) * dt; % Current time
    
    % Apply Dirichlet boundary condition at x_max
    H(end) = x(end) - exp(-r * (T - t));
    
    % Construct LHS and RHS diagonals
    main_LHS = 1 - 0.5 * dt * b;
    upper_LHS = -0.5 * dt * c(1:end-1);
    lower_LHS = -0.5 * dt * a(2:end);
    
    main_RHS = 1 + 0.5 * dt * b;
    upper_RHS = 0.5 * dt * c(1:end-1);
    lower_RHS = 0.5 * dt * a(2:end);
    
    % Adjust for Neumann boundary condition at x=1 (i=1)
    main_LHS(1) = 1 - 0.5 * dt * (a(1) + b(1));
    upper_LHS(1) = -0.5 * dt * c(1);
    main_RHS(1) = 1 + 0.5 * dt * (a(1) + b(1));
    upper_RHS(1) = 0.5 * dt * c(1);
    
    % Construct RHS vector
    rhs = zeros(n, 1);
    for i = 1:n
        if i == 1
            rhs(i) = main_RHS(i) * H(i) + upper_RHS(i) * H(i+1);
        elseif i == n
            rhs(i) = H(i); % Dirichlet BC (already set)
        else
            rhs(i) = lower_RHS(i-1) * H(i-1) + main_RHS(i) * H(i) + upper_RHS(i) * H(i+1);
        end
    end
    
    % Solve the tridiagonal system for H(1:n-1)
    A = diag(main_LHS(1:n-1)) + diag(upper_LHS(1:n-2), 1) + diag(lower_LHS(1:n-2), -1);
    H_new = zeros(n, 1);
    H_new(1:n-1) = A \ rhs(1:n-1);
    H_new(n) = x(end) - exp(-r * (T - t)); % Enforce Dirichlet BC
    
    H = H_new;
end

% Plot the option value against stock price (S/m)
plot(x, H, 'LineWidth', 1.5)
xlabel('Stock Price (S)')
ylabel('Option Value')
title('Floating Strike Lookback Call Option Value')
grid on