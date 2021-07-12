% usage:    window = raisedCosWindow_mgl(width,height,winWidth,winHeight,plateau,xDeg2pix,yDeg2pix)
% by:       Michael Jigo
% date:     05/25/18
% purpose:  Create a 2D raised cosine window using a similar format used by MGL.
%
% INPUT
% width        total width of image
% height       total height of image
% winWidth     breadth of cosine filter between minima
% winHeight    breadth of cosine filter between minima
% plateau      insert a plateau such that the window will extend through the width
%              and height of the images, with a gradual fade at the edges.

function window = raisedCosWindow_mgl(width,height,winWidth,winHeight,plateau,xDeg2pix,yDeg2pix,xCenter,yCenter)
% Currently only works with square windows
if width~=height
   error('Image must be a square; width and height must be equal.');
end

%% Defaults
if ~exist('plateau','var')
   plateau = 0;
end

% defaults for center position
if ~exist('xCenter','var')
   xCenter = 0;
end

if ~exist('yCenter','var')
   yCenter = 0;
end

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

% get size in pixels
widthPixels = round(width*xDeg2pix);
heightPixels = round(height*yDeg2pix);
widthPixels = widthPixels + mod(widthPixels+1,2);
heightPixels = heightPixels + mod(heightPixels+1,2);

%% Determine SF that would give window height and window width
winWidth = winWidth;
winHeight = winHeight;
winWidthPx = round(winWidth*xDeg2pix);
winHeightPx = round(winHeight*yDeg2pix);
winWidthPx = winWidthPx +  mod(winWidthPx+1,2);
winHeightPx = winHeightPx +  mod(winHeightPx+1,2);

% convert window widths to appropriate frequency
fWidth = 2*pi/(winWidth);
fHeight = 2*pi/(winHeight);

%% Create window
x = -width/2:width/(widthPixels-1):width/2;
y = -height/2:height/(heightPixels-1):height/2;
winX = (cos(fWidth*(x-xCenter))+1)/2;
winY = (cos(fHeight*(y-yCenter))+1)/2;

if plateau
   % make sure cosine function starts at 0 on left-most edge
   if winX(1)>0.5
      winX = (cos(fWidth*(x-xCenter)-pi)+1)/2;
      % make counter-phase cosine function
      counterX = (cos(fWidth*(x-xCenter))+1)/2;
   else
      counterX = (cos(fWidth*(x-xCenter)-pi)+1)/2;
   end

   if winY(1)>0.5
      winY = (cos(fHeight*(y-yCenter)-pi)+1)/2;
      % make counter-phase cosine function
      counterY = (cos(fHeight*(y-yCenter))+1)/2;
   else
      counterY = (cos(fHeight*(y-yCenter)-pi)+1)/2;
   end

   % sum with counter-phase cosine
   peakX = find(winX>=max(winX)-1e-6);
   peakX = min(peakX):max(peakX);
   % sum in-phase and counter-phase cosine to get plateau
   winX(peakX) = winX(peakX)+counterX(peakX);

   peakY = find(winY>=max(winY)-1e-6);
   peakY = min(peakY):max(peakY);
   % sum in-phase and counter-phase cosine to get plateau
   winY(peakY) = winY(peakY)+counterY(peakY);
else
   % display only 1 cycle of the cosinusoid
   winX(x<xCenter-winWidth/2 | x>xCenter+winWidth/2) = 0;
   winY(y<yCenter-winHeight/2 | y>yCenter+winHeight/2) = 0;
end
% take outer product to make 2D window
window = winY'*winX;
