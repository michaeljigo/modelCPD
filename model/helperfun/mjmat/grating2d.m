% Usage: g = grating2d(width,height,sf,angle,phase)
% Create a 2D grating. Adapted from mglMakeGrating

function g = grating2d(width,height,sf,angle,phase)

if ieNotDefined('angle'),angle = 0; end
if ieNotDefined('phase'),phase = 0; end

% make it so that angle of 0 is horizontal
angle = angle-90;

% calculate image paramters
phase = pi*phase/180;

% if height is nan, it means we should calculate a 1 dimensional grating
if isnan(height)
  % 1D grating (note we ignore orientation)
  x = round(-width/2):round(width/2);
  g = cosd(x*sf*2*pi+phase);
else
  % 2D grating
  % calculate orientation
  angle = pi*angle/180;
  a=cos(angle)*sf*2*pi;
  b=sin(angle)*sf*2*pi;

  % get a grid of x and y coordinates that has 
  % the correct number of pixels
  x = deg2rad(-round(width/2):round(width/2)-1);
  y = deg2rad(-round(height/2):round(height/2)-1);
  [xMesh,yMesh] = meshgrid(x,y);

  % compute grating
  g = cos(a*xMesh+b*yMesh+phase);
end
