function [ySave] = temporalFilter(x,n,tau,deltaT)
% y = temporalFilter(x,n,tau,deltaT)
%
% Cascade of exponential lowpass filters
%
% x: 3D input signal with time in the 3rd dimension
% y: output signal
% n: number of filters in the cascade
% tau: time constant
% 
% DJH 1/2019

ySave = zeros([size(x),n]);
y = zeros([size(x,1),size(x,2),n]);
for tt = 1:size(x,3)
  deltaY = (deltaT/tau) * (- y(:,:,1) + x(:,:,tt));
  y(:,:,1) = y(:,:,1) + deltaY;
  for nn = 2:n
    deltaY = (deltaT/tau) * (-y(:,:,nn) + y(:,:,nn-1));
    y(:,:,nn) = y(:,:,nn) + deltaY;
  end
  ySave(:,:,tt,:) = y;
end

return

%% Test/debug

clear all; 
%close all;

% Parameters
deltaT = 0.1;
n = [3,5];
tau = 12; % msec
sc = 1.7244; % so that response amplitude = 1 for 8Hz input

% Sampling
t = [-10:deltaT:1000-deltaT]';
[xIm,yIm,tIm] = meshgrid([1:3],[1:3],t);

% Input: sinusoidal flicker
% tf = 8/1000; % cycles/msec
% stimulus = sin(2*pi*tf*tIm);
% stimulus(:,:,t<0) = 0;

% Input: impulse
stimulus = (tIm == 0);

% Filtered output
tFiltResponses = sc * temporalFilter(stimulus,max(n),tau,deltaT);
y = tFiltResponses(:,:,:,n(1)) - tFiltResponses(:,:,:,n(2));
%std(1000/deltaT*y(2,2,:))
% 1/max(abs(y(2,2,[round(size(y,3)/2):end])))
% min(y(2,2,[round(size(y,3)/2):end]))
% max(y(2,2,[round(size(y,3)/2):end]))

% Plot it
figure(1); clf;
plot(t,squeeze(y(2,2,:)));
xlim([-10,200]);
%ylim([-1 1])
drawPublishAxis;

% Frequency response (this only makes sense for impulse input)
yCrop = squeeze(y(2,2,t>=0));
F = length(yCrop);
Fs = 1/(deltaT/1000);
frequencies = Fs/2 * linspace(0,1,F/2+1);
yFourierAmp = abs(fft(yCrop))/F;
yFourierAmp = yFourierAmp(1:F/2+1);
yFourierAmp = yFourierAmp / max(yFourierAmp);

% Plot it
figure(2); clf;
loglog(frequencies(1:100), yFourierAmp(1:100));
xlim([1,100]);
ylim([0.01,1]);
axis square;
drawPublishAxis;

% Preferred frequency
prefFreq = frequencies(find(yFourierAmp == 1))

