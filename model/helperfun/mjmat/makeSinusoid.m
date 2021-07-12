% usage:    [cosine sine time] = makeSinusoid(f,amplitude,phase,width)
% by:       Michael Jigo
% date:     03/07/18
% purpose:  Create a 1D sinusoid.

function [c, s, t] = makeSinusoid(f,amplitude,phase,width)

if ~exist('phase','var')
    phase = 0;
end

if ~exist('amplitude','var')
    amplitude = 1;
end

if ~exist('width','var')
    width = 2^7;
end

% convert f (cycles/time) to radians/time (angular frequency)
fs = width; % sampling rate, always power of 2
f = 2*pi*f/fs; % radians/sample

% set up time axis to give desired number of full cycles
t = 0:width-1;
c = amplitude*cos(f*t+phase);
s = amplitude*sin(f*t+phase);
