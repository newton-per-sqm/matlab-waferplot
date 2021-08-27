%% Example Waferplot
% Summary of example objective

load Example

A = WaferPlot(xData,yData,zData);
A.contourLevels = [-100,40,100,200,300];
A.interpMethod = 'rbf';
A.interpolateGrid;
A.createFigure;



%% Section 1 Title
% Description of first code block
a = 1;

%% Section 2 Title
% Description of second code block
b = 2;
