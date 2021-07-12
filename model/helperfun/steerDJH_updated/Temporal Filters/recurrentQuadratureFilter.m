function y = recurrentQuadratureFilter(x,lambda,omega,tau,n,deltaT)
% y = recurrentQuadratureFilter(x,lambda,omega,tau,n,deltaT)
%
% x: 3D input signal with time in the 3rd dimension (real valued)
% y: output signal (complex valued)
% omega: preferred temporal frequency
% lambda & tau: determine bandwidth
% n: number of filters in the cascade
%
% Computes:
%
%    tau dy/dt = -y + lambda (x + i imag(y)) + (1 - lambda) yhat
%    yhat = w y
%    w = 1 + i 2 pi tau omega
%
% Or:
%
%    tau dy_n/dt = -y_n + lambda y_(n-1) + (1 - lambda) yhat
%    yhat = w y_n
%
% Based on Eqs. 51 and 52 in:
% Heeger, D. J. and W. E. Mackey (2018). "ORGaNICs: A Theory of Working
% Memory in Brains and Machines." arXiv preprint arXiv:1803.06288.
% 
% DJH 1/2019

w = 1 + i*2*pi*tau*omega;

y = zeros([size(x),n]);
yTmp = zeros([size(x,1), size(x,2)]);
yHat = zeros([size(x,1), size(x,2)]);
deltaY = zeros([size(x,1), size(x,2)]);
T = size(x,3);

h = waitbar(0,'Computing...');
for tt = 1:T-1
  if (round(tt/100) == tt/100)
    waitbar(tt/T,h);
  end
  for nn = 1:n
    yHat(:,:) = w * y(:,:,tt,nn);
    yTmp(:,:) = y(:,:,tt,nn);
    if (nn==1)
      deltaY(:,:) = -yTmp + lambda * (x(:,:,tt) + i*imag(yTmp)) + (1-lambda) * yHat;
    else
      deltaY(:,:) = -yTmp + lambda * y(:,:,tt,nn-1) + (1-lambda) * yHat;
    end
    deltaY(:,:) = (deltaT./tau) * deltaY;
    y(:,:,tt+1,nn) = yTmp + deltaY;
  end
end
close(h);

return

%% Test/debug

clear all; 
%close all;

% Parameters
deltaT = 0.1; % msec
tau = 1; % msec
%lambda = 0.1;
lambda = 0.02;
omega = 8/1000; % cycles/msec
n = 3;

% Sampling
t = [-100:deltaT:1000-deltaT]';
[xIm,yIm,tIm] = meshgrid([1:3],[1:3],t);

% Input: sinusoidal flicker
% tf = 8/1000; % cycles/msec
% stimulus = sin(2*pi*tf*tIm);
% stimulus(:,:,t<0) = 0;
% yMax = 1.0;

% Input: impulse
stimulus = (tIm == 0) / deltaT;
yMax = lambda/tau;

% Filtered output
y = recurrentQuadratureFilter(stimulus,lambda,omega,tau,n,deltaT);
yCrop = squeeze(y(2,2,:,:));
yAmp = abs(yCrop);
yAmpIntegral = sum(yAmp) * deltaT

% Plot it
figure(1); clf;
subplot(n+1,1,1);
plot(t,squeeze(stimulus(2,2,:)));
ylim([0 1/deltaT]);
drawPublishAxis;
for nn = 1:n
  subplot(n+1,1,nn+1);
  plot(t,[real(yCrop(:,nn)) imag(yCrop(:,nn))]);
  hold on
  plot(t, yAmp(:,nn), 'k');
  hold off
  ylim([-yMax yMax]);
  title(['n = ',num2str(nn)]);
  drawPublishAxis;
end

% Frequency response (this only makes sense for impulse input)
yCrop2 = yCrop(t>=0,:);
F = size(yCrop2,1);
Fs = 1/(deltaT/1000);
frequencies = Fs/2 * linspace(0,1,F/2+1);
yFourierAmp = abs(fft(yCrop2))/F;
yFourierAmp = yFourierAmp(1:F/2+1,:);

% Plot it
figure(2); clf;
for nn = 1:n
  subplot(n,1,nn);
  yFourierAmpSc = yFourierAmp(:,nn) / max(yFourierAmp(:,nn));
  loglog(frequencies(2:101), yFourierAmpSc(2:101));
  xlim([1,100])
  ylim([0.01,1]);
  %axis square;
  % Preferred frequency
  prefFreq = frequencies(find(yFourierAmpSc == 1));
  title(['Preferred frequency = ',num2str(prefFreq)]);
  drawPublishAxis;
end
