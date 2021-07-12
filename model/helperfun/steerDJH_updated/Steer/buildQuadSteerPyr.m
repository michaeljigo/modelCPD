function pyr = buildQuadSteerPyr(im, pyrFRs)
% buildSteerBands: builds subbands multiscale of a multiscale
% image transform, given the frequency responses of the filters. 
%
% pyr = buildSteerBands(im, pyrFRs)
%    im: input image
%    pyrFRs = {hiFR, loFR, FRs}: see makeQqadSteerFRs
%
% Return values:
%    pyr = {highpass band, lowpass band, bandpass bands} 
%
% DJH 8/1996
% updated 1/2019

[M N] = size(im);
totalsize = M*N;

hiFR = pyrFRs{1};
loFR = pyrFRs{2};
FRs = pyrFRs{3};
numLevels = size(FRs,3);
numOrientations = size(FRs,4);

fourier = fftshift(fft2(im));

bands = zeros(size(FRs));
for level = 1:numLevels
  for orientation = 1:numOrientations;
    FR = FRs(:,:,level,orientation);
    band = real(ifft2(fftshift(fourier .* FR)));
    bandR = real(ifft2(fftshift(fourier .* real(FR))));
    bandI = real(ifft2(fftshift(fourier .* (i*imag(FR)))));
    bands(:,:,level,orientation) = bandR + i*bandI;
  end
end

hi = real(ifft2(fftshift(fourier.*hiFR)));
lo = real(ifft2(fftshift(fourier.*loFR)));

pyr = {hi, lo, bands};

