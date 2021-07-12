% usage:    [cosine sine] = makeSinusoid2d(f,amplitude,phase,width,height,angle)
% by:       Michael Jigo
% date:     03/07/18
% purpose:  Create a 2D sinusoid.

function [c, s] = makeSinusoid2d(f,amplitude,phase,width,height,angle)

if ieNotDefined('phase')
    phase = 0;
end

if ieNotDefined('amplitude')
    amplitude = 1;
end

if ieNotDefined('width')
    width = 2^7;
end

if ieNotDefined('height')
    height = width;
end

if ieNotDefined('angle')
   angle = 90;
end

% convert f (cycles/time) to radians/sample
fs = width; % sampling rate, always power of 2
f = 2*pi*f/fs; % radians/sample

% set up time axis to give desired number of full cycles
t = 0:width-1;

% make it so that angle of 0 is horizontal
angle = angle-90;

% calculate orientation
angle = pi*angle/180;
a=cos(angle)*f;
b=sin(angle)*f;

% get a grid of x and y coordinates that has
% the correct number of pixels
[xMesh,yMesh] = meshgrid(t,0:height-1);

% compute grating
c = amplitude*cos(a*xMesh+b*yMesh+phase);
s = amplitude*sin(a*xMesh+b*yMesh+phase);
