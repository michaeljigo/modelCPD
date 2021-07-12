function pyrFRs = makeQuadSteerFRs(dims,numLevels,numOrientations,bandwidth)
% Makes the frequency responses of the filters for a multiscale image
% transform.
%
% pyrFRs = makeQuadSteerFRs(dims,numLevels,numOrientations,bandwidth)
%    [M N]: image size
%    numLevels: number of levels/scales.
%    numOrientations: number of orientation subbands at each scale.
%    bandwidth: spatial frequency bandwidth in octaves
%
% Return values: pyrFRs = {hiFR, loFR, FRs};
%    hiFR: frequency response of the highpass subband.
%    loFR: frequency response of the lowpass subband.
%    FRs: cell array that contains frequency responses of the bandpass subbands.
%
% DJH 8/1996
% updated 1/2019

p = numOrientations-1;
Const = sqrt((2^(2*p)*(factorial(p))^2)/(factorial(2*p)*(p+1)));
[f1,f2] = freqspace(dims);
[wx,wy] = meshgrid(f1,f2);
r = sqrt(wx.^2 + wy.^2);
theta = atan2(wy,wx);

% bandpass bands
FRs = zeros([dims,numLevels,numOrientations]);
for level = 1:numLevels
  for orientation = 1:numOrientations;
    thetaOffset = (orientation-1) * pi/numOrientations; % ORIGINAL -- 01/28/21
    CtrFreq = pi/2^(level*bandwidth);
    band = i^p * Const * cos(theta-thetaOffset).^p .* logRaisedCos(r,CtrFreq,bandwidth);
    FRs(:,:,level,orientation) = abs(band) + band;
    %FRs(:,:,level,orientation) = abs(band);
    %FRs(:,:,level,orientation) = band;
  end
end

% hi band
hiFR = logRaisedCosHi(r,pi/(2^bandwidth),bandwidth);

% lo band
loFR = logRaisedCosLo(r,pi/(2^(bandwidth*numLevels)),bandwidth);

pyrFRs = {hiFR, loFR, FRs};
