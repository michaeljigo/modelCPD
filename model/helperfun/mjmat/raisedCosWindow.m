% usage:    window = raisedCosWindow(space,width)
% by:       Michael Jigo
% date:     05/25/18
% purpose:  Create a 2D raised cosine window
%
% INPUT
% space     scalar specifying the total width of image
% width     scalar specifying the width at 0 of the window
%           the 2D window will have a symmetric width and height
% plateau   sums raised cosine window with the inputted width until it covers space 

function window = raisedCosWindow(space,width,plateau)
if ~exist('plateau','var')
   plateau = 0;
end

%% Determine SF of cosinusoid based on width input
width = width;
f = 2*pi/(space*(width/space)); % width/space division converts from degrees to pixels

%% Create window
space = linspace(-space/2,space/2,space);
win = (cos(f*space)+1)/2;
if plateau
   % make sure cosine function starts at 0 on left-most edge
   if win(1)>0.5
      win = (cos(f*space-pi)+1)/2;
      % make counter-phase cosine function
      counterWin = (cos(f*space)+1)/2;
   else
      counterWin = (cos(f*space-pi)+1)/2;
   end
   % find the first and last peaks in the cosine function
   peak = find(win>=max(win)-1e-6);
   peak = min(peak):max(peak);
   % sum in-phase and counter-phase cosine to get plateau
   win(peak) = win(peak)+counterWin(peak);
else
   % only display 1 cycle
   win(space<-width/2 | space>width/2) = 0;
end
keyboard
% take outer produce to make 2D window
window = win'*win;
