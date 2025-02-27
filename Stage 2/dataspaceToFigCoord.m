function figPt = dataspaceToFigCoord(ax, xData, yData)
    % Get axes position in normalized figure coords
    axPos = get(ax,'Position');      % [x0 y0 width height]
    axLimits = axis(ax);            % [xMin xMax yMin yMax]

    % Normalized data range in x and y
    xNorm = (xData - axLimits(1)) / (axLimits(2) - axLimits(1));
    yNorm = (yData - axLimits(3)) / (axLimits(4) - axLimits(3));

    % Now map that to figure coords
    figX = axPos(1) + xNorm * axPos(3);
    figY = axPos(2) + yNorm * axPos(4);

    figPt = [figX, figY];
end