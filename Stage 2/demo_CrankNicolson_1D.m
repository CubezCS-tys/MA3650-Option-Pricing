function demo_CrankNicolson_1D()
% DEMO_CRANKNICOLSON_1D
%
% A self-contained example of solving the 1D heat/diffusion equation:
%
%     u_t = alpha * u_{xx},   0 < x < 1,  0 < t <= T
%
% with non-zero boundary conditions:
%     u(0, t) = LeftBC(t)
%     u(1, t) = RightBC(t)
%
% using the Crank-Nicolson scheme in space (i=1..M+2) and time (n=0..N).
%
% We highlight how the boundary "aux" vector is built each time step
% and show an animation of the solution, marking boundary vs. interior.

    % -----------------------------
    % 1) PDE & Discretization Setup
    % -----------------------------
    alpha = 0.1;          % diffusion coefficient
    T     = 0.2;          % final time
    M     = 8;            % number of interior points
    dx    = 1/(M+1);      % spatial step (M+2 total points, i=0..M+1)
    xVals = linspace(0,1,M+2);  % x = 0, dx, 2dx, ..., 1

    N     = 40;           % number of time steps
    dt    = T/N;          % time step
    tVals = linspace(0,T,N+1);  % t = 0..T

    % For Crank-Nicolson, we define r = alpha * dt / (2 dx^2).
    r = alpha * dt / (2*dx^2);

    % -------------------------------
    % 2) Boundary Conditions & Init
    % -------------------------------
    % Let's define some time-dependent Dirichlet BC:
    %   u(0,t)   = 1 + sin(4*pi*t)         (Left boundary)
    %   u(1,t)   = 2                       (Right boundary)
    %
    % We'll define functions in local handles:
    LeftBC  = @(t) 1 + sin(4*pi*t);
    RightBC = @(t) 2;

    % Initial condition: let's do u(x,0) = sin(pi*x) (just an example)
    IC = @(x) sin(pi*x);

    % We'll store the full solution in a matrix "U".
    %   size(U) = (M+2) x (N+1)
    % Indices: U(i, n+1) means i=0..M+1 in space, n=0..N in time
    % We'll use 1-based indexing in MATLAB, so i=1..(M+2).
    U = zeros(M+2, N+1);

    % Fill initial condition at t=0
    for i = 1:(M+2)
        U(i,1) = IC(xVals(i));
    end

    % Impose the boundary condition at the initial time:
    U(1,1)      = LeftBC(0);  % i=1 => x=0
    U(M+2,1)    = RightBC(0); % i=M+2 => x=1

    % -------------------------------------
    % 3) Build Crank-Nicolson Matrices A,B
    % -------------------------------------
    % We solve (I + r*A)*U^{n+1} = (I - r*A)*U^n + boundary_terms
    % For standard second-difference, let's build A for interior i=2..(M+1).
    %   A is (M) x (M), indexing interior only.
    %   We'll call iInterior = 2:(M+1) in the full array.
    e = ones(M,1);
    A = spdiags([e -2*e e], -1:1, M, M);
    % A here represents the second difference. With factor r, we'd do r*A.
    % But let's keep it separate for clarity.

    I = speye(M);
    LHS = (I + r*A);   % factor on U^{n+1} (interior)
    RHS = (I - r*A);   % factor on U^{n}   (interior)

    % We'll do an LU decomposition for efficiency:
    [Lfac, Ufac] = lu(LHS);

    % -----------
    % 4) Solve
    % -----------
    % For n=0..N-1, step from t_n to t_{n+1}.
    for n = 1:N
        tNow   = tVals(n);
        tNext  = tVals(n+1);

        % Impose BC at current/next time:
        leftNow  = LeftBC(tNow);
        leftNext = LeftBC(tNext);
        rightNow  = RightBC(tNow);
        rightNext = RightBC(tNext);

        U(1,n)     = leftNow;     % for completeness
        U(M+2,n)   = rightNow;    %
        U(1,n+1)   = leftNext;
        U(M+2,n+1) = rightNext;

        % Extract the interior from the old time step: i=2..(M+1)
        Uold_interior = U(2:(M+1), n);

        % Build the RHS vector:
        %   RHS * Uold_interior
        b = RHS * Uold_interior;

        % 4a) Add "aux" from boundary values:
        %  Typically with Dirichlet BC: at i=2 => depends on left boundary,
        %  and at i=M+1 => depends on right boundary.
        %
        %  The standard CN derivation for non-zero BC leads to terms like:
        %       + r * [LeftBC(t_n) + LeftBC(t_{n+1})]
        % in the first interior node,
        %       + r * [RightBC(t_n) + RightBC(t_{n+1})]
        % in the last interior node.
        %
        % So let's do that:
        aux = zeros(M,1);  % same size as interior
        aux(1)   = r*(leftNow + leftNext);      % i=2
        aux(end) = r*(rightNow + rightNext);    % i=M+1

        b = b + aux;

        % 4b) Solve LHS * Unew_interior = b
        Unew_interior = Ufac \ (Lfac \ b);

        % 4c) Store in U
        U(2:(M+1), n+1) = Unew_interior;
    end

    % -----------------------------------------
    % 5) Visualization: Animate the Time Steps
    % -----------------------------------------
    figure('Name','Crank-Nicolson 1D Demo','Position',[100 100 900 500]);

    % We'll create two subplots:
    %   (A) The solution profile u(x) at each time step
    %   (B) A space-time grid with boundary highlights

    % (A) Plot the final solution for reference
    subplot(1,2,1);
    hLine = plot(xVals, U(:,1), 'b-o','LineWidth',1.5);
    hold on;
    hLeftBC  = plot(xVals(1), U(1,1),'rs','MarkerFaceColor','r');
    hRightBC = plot(xVals(end), U(end,1),'ks','MarkerFaceColor','k');
    ylim([min(U(:))-0.1, max(U(:))+0.1]);
    xlim([0,1]);
    grid on; box on;
    xlabel('x');
    ylabel('u(x,t)');
    title('Solution Profile Over Time');

    % We'll step through time and update hLine, hLeftBC, hRightBC.

    % (B) Build a 2D "grid" visualization in (i,n) or (x-index, time-index).
    % We'll just scatter the points for i=1..M+2, n=1..N+1,
    % then highlight boundary with color. For advanced arrow-labelling,
    % see the previous examples with annotation arrows.

    subplot(1,2,2);
    [iGrid, nGrid] = meshgrid(1:(M+2), 1:(N+1)); % i=1..(M+2), n=1..(N+1)
    iFlat = iGrid(:);
    nFlat = nGrid(:);
    scatter(iFlat, nFlat, 25,'filled','MarkerFaceColor','k','MarkerFaceAlpha',0.1);
    hold on;
    axis([0 M+3 0 N+2]);
    set(gca,'YDir','reverse');
    xlabel('Spatial index i');
    ylabel('Time step n');
    title('Space-Time Grid: boundary in red/blue');
    grid on; box on;

    % We can highlight boundary columns: i=1 and i=M+2
    scatter( ones(N+1,1),         (1:N+1)',  70,'ro','filled');
    scatter( (M+2)*ones(N+1,1),   (1:N+1)',  70,'bo','filled');

    % We'll also place text with the boundary values U(1,n) & U(M+2,n)
    % in an animation.

    % Create handles for boundary text:
    leftText  = text(2, 1, '', 'Color','r','FontWeight','bold');
    rightText = text(M+2-1, 1, '', 'Color','b','FontWeight','bold',...
                     'HorizontalAlignment','right');

    % 5a) Animate in time
    for n = 0:N
        % Update the solution profile:
        set(hLine, 'XData', xVals, 'YData', U(:,n+1));
        set(hLeftBC,  'XData', xVals(1), 'YData', U(1,n+1));
        set(hRightBC, 'XData', xVals(end), 'YData', U(end,n+1));
        subplot(1,2,1);
        title(sprintf('Solution at t=%.3f (n=%d)', tVals(n+1), n));

        % Update boundary text on the second subplot
        subplot(1,2,2);
        % Move the text near (i=1,n+1) and (i=M+2,n+1)
        set(leftText,  'Position',[2, n+1], ...
            'String', sprintf('%.2f', U(1,n+1)));
        set(rightText, 'Position',[M+2-1, n+1], ...
            'String', sprintf('%.2f', U(M+2,n+1)));

        drawnow;
        pause(0.2);
    end

end
