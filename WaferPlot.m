%#ok<*PROP>
classdef WaferPlot < handle & ... 
    matlab.mixin.CustomDisplay & ...
    matlab.mixin.Copyable & ...
    matlab.mixin.SetGet 
%WAFERPLOT - Class for plotting nice wafer surface maps.
%   More information will come soon.
%
%   (c) Pascal Muster, Germany, 2021

%% Properties
properties(Constant)
    Version string = '0.1a';
end

properties(Access=public)
    
    % Plot Objects
    Figure matlab.ui.Figure
    Tile matlab.graphics.layout.TiledChartLayout
    Axes matlab.graphics.axis.Axes
    ColorBar matlab.graphics.illustration.ColorBar
    Renderer = 'OpenGL';
    Handles % Drawed lines and handles
    
end

properties(SetObservable,Access=public)
    
    % Plot properties
    figSize (1,2) double {mustBePositive} = [16,15];
    fontSize (1,1) double {mustBePositive} = 12;
    waferSize (1,1) double {mustBePositive} = 150;
    waferUnit (1,:) char = 'mm'; 
    
    % Interpolation properties
    interpMethod (1,:) char ...
        {mustBeMember(interpMethod,{'rbf','griddata'})} = 'rbf';
    interpType (1,:) char  = 'linear';
    interpResolution (1,1) int16 ...
        {mustBeInRange(interpResolution,100,1000)} = 300;
    
    % Extrapolation properties
    extrapolate (1,1) logical = true;
    extrapolationAlpha (1,1) double ...
        {mustBeInRange(extrapolationAlpha,0,1)} = 0.65;
    
    plotType (1,:) char = 'surfc';
    
    contourLevels double = 10;
    contourLabels (1,1) logical = true
    
    xLabel (1,:) char = 'x';
    yLabel (1,:) char = 'y';
    cLabel (1,:) char = 'Value';
    
end

properties(Access=protected)
    DT % Delaunay Triangles
    CH % Convex Hull
    
    % Data Arrays
    xData (:,1) double {mustBeReal, mustBeFinite}
    yData (:,1) double {mustBeReal, mustBeFinite}
    zData (:,1) double {mustBeReal, mustBeFinite}
    
    % Grids
    X double
    Y double
    ZI double
    ZE double
    ZA double
end

%% Methods

methods(Access=public)
% Public Methods, Constructor
   
    function obj = WaferPlot(xData, yData, zData, varargin)
    %WAFERPLOT - Class for plotting nice wafer surface maps.
    %   More information will come soon.
    %
    %   (c) Pascal Muster, Germany, 2021
    
        if nargin == 0
            % Example data
            xData = linspace(-68,68,6);
            yData = linspace(-68,68,6);
            [xData, yData] = meshgrid(xData, yData);
            zData = 0.1*xData.^2-0.2*yData+0.1*randn(size(xData))+100;
            
            % Remove values outside of wafer radius
            zData(sqrt(xData.^2+yData.^2)>75) = nan;
            xData(isnan(zData)) = [];
            yData(isnan(zData)) = [];
            zData(isnan(zData)) = [];
            
            % Reshape to column vectors
            xData = reshape(xData,[],1);
            yData = reshape(yData,[],1);
            zData = reshape(zData,[],1);
        end
    
        % Save data arrays to object
        obj.xData = xData;
        obj.yData = yData;
        obj.zData = zData;
        
        % Delaunay Grid and Convex Hull
        obj.DT = delaunayTriangulation(xData, yData); % Delaunay Grid
        obj.CH = convexHull(obj.DT); % Convex Hull of Delaunay Grid
        
    end
    
end

methods(Access=public)
    
    function obj = interpolateGrid(obj)
    %INTERPOLATEGRID generates interpolation by griddata or rbf.
        
        [X, Y] = ndgrid(... % Create ndgrid in x and y about the full wafer size!
            linspace(-obj.waferSize/2, +obj.waferSize/2, obj.interpResolution)); 
        
        switch lower(obj.interpMethod)
            
            case 'griddata'
                
                % Scattered Interpolation by griddata
                ZI = griddata(obj.xData, obj.yData, obj.zData, X, Y, obj.interpType);
                ZE = nan(size(ZI));
                ZA = ZI;
        
            case 'rbf'
                
                % Scattered Interpolation by rbf
                rbf = rbfcreate([obj.xData, obj.yData]', obj.zData', ...
                    'RBFFunction', obj.interpType);
                ZA = rbfinterp([X(:), Y(:)]', rbf);
                
                % Extract interpolated and extrapolated values
                IDX = isnan(pointLocation(obj.DT, X(:), Y(:)));
                
                ZI = ZA; ZI(IDX) = nan;
                ZE = ZA; ZE(~IDX) = nan;
                
                % Reshape values to grid format
                ZI = reshape(ZI, size(X));
                ZE = reshape(ZE, size(X));
                ZA = reshape(ZA, size(X));
                ZE(sqrt(X.^2+Y.^2) > obj.waferSize/2) = nan;
                ZA(sqrt(X.^2+Y.^2) > obj.waferSize/2) = nan;
                
            otherwise
                error('Unknown interpolation method "%s".', obj.interpMethod);
        end
        
        if obj.extrapolate == false
            obj.ZE = nan(size(obj.ZE));
        end
        
        % Save X and Y ndgrids
        obj.X = X;
        obj.Y = Y;
        
        % Save interpolated (ZI) and extrapolated (ZE) values
        obj.ZI = ZI;
        obj.ZE = ZE;
        obj.ZA = ZA;
        
    end
end

methods(Access=public)
%Protected Methods

    function obj = createFigure(obj)
    %CREATEFIGURE generates basic WaferPlot figure and axes appearance.
    
        obj.Figure = figure( ...
            'Color', 'white', ...
            'Units', 'centimeters', ...
            'ToolBar', 'auto', ...
            'Renderer', obj.Renderer);
        
        obj.Figure.Position = [3,3,obj.figSize];
        
        % Figure Layout
        obj.Tile = tiledlayout(obj.Figure, 1, 1, 'Padding', 'compact');
        
        % Figure Axes Appearance
        obj.Axes = nexttile(obj.Tile);
        obj.Axes.FontSize = obj.fontSize;
        
        axtoolbar(obj.Axes, {'export','datacursor','brush'});
        
        obj.Axes.DataAspectRatio = [1,1,1];
        obj.Axes.XLabel.String = ['x [', obj.waferUnit, ']'];
        obj.Axes.YLabel.String = ['y [', obj.waferUnit, ']'];
        
        obj.Axes.XLim = round([-obj.waferSize/2, +obj.waferSize/2],-1);
        obj.Axes.YLim = round([-obj.waferSize/2, +obj.waferSize/2],-1);
        
        % No outline box
        obj.Axes.Box = 'off';
        obj.Axes.NextPlot = 'add';
        
        % Draw surface plot
        obj.drawSurface;
        
        % Draw grid
        obj.drawGrid(obj.Axes, obj.waferSize/2, ...
            'LineWidth', 0.8);
        
        % Draw wafer line circle
        obj.drawCircle(obj.Axes, obj.waferSize/2, '-', ...
            'LineWidth', 1.5, 'Color', 'black');
        
        obj.drawConvexHull(obj.Axes, obj.xData, obj.yData, '-', ...
            'LineWidth', 1.5, 'Color', 'black');
        
        % Draw measurement points
        obj.drawPoints(obj.Axes, obj.xData, obj.yData, 'ko');
        
        % Create despined axes
        offsetAxes(obj.Axes)
        
        obj.Axes.LineWidth = 1;
        obj.Axes.TickDir = 'out';
        obj.Axes.XColor = 'black';
        obj.Axes.YColor = 'black';
        
    end
    
    function obj = drawSurface(obj)
    %DRAWDELAUNAY draws a triangulated grid.
    
        switch obj.plotType
            
            case 'contourf'
                
                [~, h1] = contourf(obj.X, obj.Y, obj.ZI, obj.contourLevels, 'LineStyle', 'None');
                [~, h2] = contourf(obj.X, obj.Y, obj.ZE, h1.LevelList, 'LineStyle', 'None');
                [C, h3] = contour(obj.X, obj.Y, obj.ZI, h1.LevelList, 'Color', 'k');
                [~, h4] = contour(obj.X, obj.Y, obj.ZE, h1.LevelList, 'Color', 'k');
        
                % make extrapolation transparent
                drawnow; % update graphical objects
                h2Fills = h2.FacePrims;  % array of TriangleStrip objects
                [h2Fills.ColorType] = deal('truecoloralpha');  % default = 'truecolor'
                for i = 1 : numel(h2Fills)
                    h2Fills(i).ColorData(4) = 255*obj.extrapolationAlpha;   % default=255
                end
        
            case 'surfc'
                
                [~, h1] = contourf(obj.X, obj.Y, obj.ZI, 200, 'LineStyle', 'None');
                if ~all(isnan(obj.ZE), 'all')
                    [~, h2] = contourf(obj.X, obj.Y, obj.ZE, h1.LevelList, 'LineStyle', 'None');
                end
                [C, h3] = contour(obj.X, obj.Y, obj.ZI, obj.contourLevels, 'Color', 'k');
                if ~all(isnan(obj.ZE), 'all')
                    contour(obj.X, obj.Y, obj.ZE, h3.LevelList, 'Color', 'black');
                end
                
                if ~all(isnan(obj.ZE), 'all')
                    % make extrapolation transparent
                    drawnow; % update graphical objects
                    h2Fills = h2.FacePrims;  % array of TriangleStrip objects
                    [h2Fills.ColorType] = deal('truecoloralpha');  % default = 'truecolor'
                    for i = 1 : numel(h2Fills)
                        h2Fills(i).ColorData(4) = 255*obj.extrapolationAlpha;   % default=255
                    end
                end
                
            otherwise
                error('Unknown plot type %s', obj.plotType);
        end
                
        % contour labels
        clabel(C,h3)
        
        obj.ColorBar = colorbar(obj.Axes);
        obj.ColorBar.LineWidth = 1;
        obj.ColorBar.TickDirection = 'out';
        obj.ColorBar.Label.String = obj.cLabel;
        obj.ColorBar.FontSize = obj.fontSize;
        
        addlistener(obj,'contourLevels','PostSet',@obj.handlePropEvents);
        
    end
    
end

methods(Static=true, Access=public)
%Static Protected Methods

   function handlePropEvents(source,event)
   %HANDLEPROPEVENTS listens on SetObservable property changes.
      switch source.Name 
         case 'contourLevels'
            obj = event.AffectedObject;
            obj.drawSurface;
      end
   end

    function lineHandle = drawCircle(axes, radius, varargin)
    %DRAWCIRCLE draws a circle of specified radius. 
    
        phi = linspace(0,2*pi,2*365);
        phi(end+1) = phi(1);
        
        x = radius*cos(phi);
        y = radius*sin(phi);
        
        lineHandle = plot(axes, x, y, varargin{:});
        lineHandle.HandleVisibility = 'off';
        
        if nargout == 0
            clear lineHandle
        end
        
    end
    
    function lineHandle = drawDelaunayGrid(axes, x, y, varargin)
    %DRAWDELAUNAY draws a triangulated grid. 
    
        DT = delaunay(x, y);
        lineHandle = triplot(DT, x, y, 'Parent', axes, varargin{:});
        lineHandle.HandleVisibility = 'off';

    end

    function lineHandle = drawConvexHull(axes, x, y, varargin)
    %DRAWDELAUNAY draws a triangulated grid. 
        
        DT = delaunayTriangulation(x, y);
        CH = convexHull(DT);
        lineHandle = plot(axes, DT.Points(CH,1), DT.Points(CH,2), varargin{:});
        lineHandle.HandleVisibility = 'off';
        
    end
    
    function lineHandle = drawPoints(axes, x, y, varargin)
    %DRAWDELAUNAY draws a triangulated grid. 
        
        lineHandle = plot(axes, x, y, varargin{:});
        
    end
    
    function lineHandle = drawGrid(axes, radius, varargin)
    %DRAWGRID draws a major radial grid to the figure.
    
        xticks = axes.XTick(-radius < axes.XTick & axes.XTick < radius);
        yticks = axes.YTick(-radius < axes.YTick & axes.YTick < radius);
        
        if ~all(xticks==yticks)
            error('X- and Y-Ticks not equivalent. Aborting grid creation.');
        end
        
        % Draw circles
        phi = linspace(0,2*pi,2*365);
        phi(end+1) = phi(1);
        
        [color, varargin] = getPositional(varargin, 'Color', 0.75*[1,1,1]);
        
        for i = 1:numel(xticks)
            x = xticks(i)*cos(phi);
            y = xticks(i)*sin(phi);
            lineHandle = plot(axes, x, y, 'Color', color, varargin{:});
            lineHandle.HandleVisibility = 'off';
        end
        
        % Draw lines
        phi = 0:(30*pi/180):pi;
        phi(end) = [];
        
        t = linspace(-radius,radius);
        
        for i = 1:numel(phi)
            x = t*cos(phi(i));
            y = t*sin(phi(i));
            lineHandle = plot(axes, x, y, 'Color', color, varargin{:});
            lineHandle.HandleVisibility = 'off';
        end
        
    end
    
end
end

%% Additional Functions
function [value, remainingargs] = getPositional(args, name, default)
%GETPOSITIONAL extracts optional argument out of a cell array.

    if nargin < 3
        default = nan;
    end

    remainingargs = {};
    skipping = false;
    value = default;

    if iscell(args)
        
        % Extract cell value from array, if exists.
        for i = 1:length(args)
            if strcmpi(args{i}, name)
                value = args{i+1};
                skipping = true;
            elseif skipping
                skipping = false;
            else
                remainingargs{end+1} = args{i};
            end
        end

    elseif isstruct(args)
        
        % Extract fields and find field name.
        fields = fieldnames(args);
        idx = find(contains(lower(fields), lower(name)));

        if ~isempty(idx)
            value = args.(fields{idx});
        end

        % Don't change structure
        remainingargs = args;

    end

end
