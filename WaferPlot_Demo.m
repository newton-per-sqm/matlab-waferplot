%% Example Waferplot
% Summary of example objective

load Example

A = WaferPlot(xData,yData,zData);
A.contourLevels = [-100,-50,0,50,100,200,300,400];
A.interpMethod = 'rbf';
A.interpType = 'cubic';
A.interpolateGrid;
A.createFigure;


%% Section 1 Title
% Description of first code block
a = 1;

%% Section 2 Title
% Description of second code block
b = 2;
