function drawArrowInDataCoords(x, y, colorSpec, headSize)
    % x, y are vectors [xStart, xEnd], [yStart, yEnd] in *data* units
    % colorSpec e.g. 'r', 'b', 'k'
    % headSize e.g. 0.05..1 (affects arrowhead size)
    %
    % We'll convert data coords -> normalized figure coords
    ax = gca;
    % [pixelX,pixelY] = ds2nfu(x,y) -> old approach
    % A simpler approach is ' annotation() ' with 'textarrow' 
    % but that needs normalized figure coords. We'll do a quick conversion:
    startPoint = dataspaceToFigCoord(ax, x(1), y(1));
    endPoint   = dataspaceToFigCoord(ax, x(2), y(2));

    annotation('arrow','Color',colorSpec,...
               'X',[startPoint(1), endPoint(1)], ...
               'Y',[startPoint(2), endPoint(2)], ...
               'LineWidth',1.5, ...
               'HeadLength',10*headSize,...
               'HeadWidth',7*headSize,...
               'Tag','boundaryArrow');
end