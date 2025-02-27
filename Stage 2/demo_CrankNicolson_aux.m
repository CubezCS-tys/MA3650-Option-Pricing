function demo_CrankNicolson_aux
% DEMO_VISUALIZE_AUX_GRID_WITH_VALUES
%
% Demonstrates a 2D (i,j) grid for a Crank-Nicolson scheme, highlighting
% which boundary points (in "space" i=1 or i=M+2) at times (j, j+1)
% feed into the 'aux' vector. Additionally displays the numeric values
% matval(1,j), matval(1,j+1), matval(M+2,j), matval(M+2,j+1).
%
% M+2 is total space points, N+1 is total time steps.
% We loop backward in time j = N..1, each time highlighting boundary.

    clear; clc; close all;

    %---------------------------
    % 1) Define grid dimensions
    %---------------------------
    M = 5;           % "interior" = i=2..M => total space pts = M+2 = 7
    N = 6;           % total time steps = N+1 => 7

    % Create a meshgrid of (i, j) for plotting
    [iGrid, jGrid] = meshgrid(1:(M+2), 1:(N+1));
    iFlat = iGrid(:);
    jFlat = jGrid(:);

    %-------------------------------------
    % 2) Example 'matval' array of values
    %-------------------------------------
    % In a real solver, matval(i,j) stores your PDE solution at space i, time j.
    % We'll just fill with random values for demonstration:
    matval = rand(M+2, N+1);

    %--------------------------------
    % 3) Basic figure: show all grid
    %--------------------------------
    figure('Name','Crank-Nicolson Boundary Points','Position',[100 100 800 600]);
    hold on;
    scatter(iFlat, jFlat, 40, 'filled', 'MarkerFaceColor','k','MarkerFaceAlpha',0.2);

    xlabel('Space index (i)');
    ylabel('Time index (j)');
    title('Grid of (i,j) with boundary aux contributions highlighted');
    axis([0 M+3 0 N+2]);  % a bit of padding
    set(gca,'YDir','reverse');  % so j=1 is near the top, if desired
    grid on; box on;

    %-------------------------------------------------
    % 4) Animate backward time steps j = N down to 1
    %-------------------------------------------------
    for j = N:-1:1
        % For each time step j, the boundary 'aux' uses:
        %   V(1, j) + V(1, j+1)       -> aux(1)
        %   V(M+2, j) + V(M+2, j+1)   -> aux(end)

        % Remove old highlight items before drawing new:
        delete(findall(gcf,'Tag','boundaryHighlight'));
        delete(findall(gcf,'Tag','boundaryArrow'));
        delete(findall(gcf,'Tag','boundaryValueText'));

        %---------------------------------------------
        % a) Highlight boundary points: i=1, i=M+2
        %---------------------------------------------
        % Left boundary: (i=1, j) and (i=1, j+1)
        leftX = [1, 1];
        leftY = [j, j+1];
        scatter(leftX, leftY, 100, 'ro','filled','Tag','boundaryHighlight');

        % Right boundary: (i=M+2, j) and (i=M+2, j+1)
        rightX = [M+2, M+2];
        rightY = [j, j+1];
        scatter(rightX, rightY, 100, 'bo','filled','Tag','boundaryHighlight');

        %----------------------------------------------
        % b) Label numeric values at these boundary pts
        %----------------------------------------------
        % left boundary (i=1)
        text(1 + 0.2, j,     sprintf('%.3f',matval(1,j)), ...
            'Color','r','Tag','boundaryValueText','FontWeight','bold');
        text(1 + 0.2, j+1,   sprintf('%.3f',matval(1,j+1)), ...
            'Color','r','Tag','boundaryValueText','FontWeight','bold');

        % right boundary (i=M+2)
        text((M+2) - 0.8, j,   sprintf('%.3f',matval(M+2,j)), ...
            'Color','b','Tag','boundaryValueText','FontWeight','bold',...
            'HorizontalAlignment','right');
        text((M+2) - 0.8, j+1, sprintf('%.3f',matval(M+2,j+1)), ...
            'Color','b','Tag','boundaryValueText','FontWeight','bold',...
            'HorizontalAlignment','right');

        %----------------------------------------------
        % c) Draw arrows to show where aux is added
        %----------------------------------------------
        % The first interior node is i=2 => gets aux(1)
        % The last interior node is i=M+1 => gets aux(end)
        yMid = j + 0.5;  % between j and j+1
        drawArrowInDataCoords([1.1, 2 - 0.1], [yMid, yMid], 'r', 0.8);
        drawArrowInDataCoords([ (M+2) - 0.1, (M+1) + 0.1 ], [yMid, yMid], 'b', 0.8);

        % Optional arrow label
        text(1.5, yMid - 0.1, 'aux(1)', ...
            'Color','r','Tag','boundaryArrow','FontWeight','bold','FontSize',10,...
            'HorizontalAlignment','center');
        text(M+1.5, yMid - 0.1, 'aux(end)', ...
            'Color','b','Tag','boundaryArrow','FontWeight','bold','FontSize',10,...
            'HorizontalAlignment','center');

        %-------------------------------------
        % d) Pause to visualize each time step
        %-------------------------------------
        drawnow;
        pause(0.8);
    end
end