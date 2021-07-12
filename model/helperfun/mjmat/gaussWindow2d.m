% Usage: g = gaussWindow2d(width,height,sdx,sdy,xCenter,yCenter)
% Creates a 2D gaussian. Adapted from mglMakeGaussian

function g = gaussWindow2d(width,height,sdx,sdy,xCenter,yCenter)

if ieNotDefined('xCenter'),xCenter = 0; end
if ieNotDefined('yCenter'),yCenter = 0; end

% get a grid of x and y coordinates
x = 0:width-1;
y = 0:height-1;
[xMesh,yMesh] = meshgrid(x,y);

% compute gaussian window
g = exp(-(((xMesh-xCenter).^2)/(2*(sdx^2))+((yMesh-yCenter).^2)/(2*(sdy^2))));
g(g(:)<0.01) = 0;
