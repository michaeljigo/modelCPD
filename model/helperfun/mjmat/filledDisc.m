% usage:    disc = filledDisc(width,radius,bgColor,discColor,<xCenter>,<yCenter>,
%                  <xDeg2pix>,<yDeg2pix>)
% by:       Michael Jigo
% date:     07/22/18
% purpose:  Create a 2D disc
%
% INPUT
% width        total width of image (width will also equal height)
% radius       radius of circle
% xCenter      horizontal center of circle (default=0)
% yCenter      vertical center of circle (default=0)
% bgColor      color outside radius of disc
% discColor    color within radius of disc

function disc = filledDisc(width,radius,bgColor,discColor,xCenter,yCenter,xDeg2pix,yDeg2pix)
%% Defaults
% defaults for xDeg2pix
if ~exist('xDeg2pix','var')
   if isempty(mglGetParam('xDeviceToPixels'))
      disp(sprintf('(makeGrating) mgl is not initialized'));
      return
   end
   xDeg2pix = mglGetParam('xDeviceToPixels');
end

% defaults for yDeg2pix
if ~exist('yDeg2pix','var')
   if isempty(mglGetParam('yDeviceToPixels'))
      disp(sprintf('(makeGrating) mgl is not initialized'));
      return
   end
   yDeg2pix = mglGetParam('yDeviceToPixels');
end

if ~exist('xCenter','var') || isempty(xCenter)
   xCenter = 0;
end

if ~exist('yCenter','var') || isempty(yCenter)
   yCenter = 0;
end

% get size in pixels
widthPixels = round(width*xDeg2pix);
heightPixels = round(width*yDeg2pix);
widthPixels = widthPixels+mod(widthPixels+1,2);
heightPixels = heightPixels+mod(heightPixels+1,2);

%% Create disk
x = -width/2:width/(widthPixels-1):width/2;
y = -width/2:width/(heightPixels-1):width/2;
[xMesh,yMesh] = meshgrid(x,y);
disc = double((xMesh-xCenter).^2+(yMesh-yCenter).^2<=radius.^2);
disc(disc==0) = bgColor;
disc(disc==1) = discColor;
