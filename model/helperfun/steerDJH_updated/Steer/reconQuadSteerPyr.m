function result = reconQuadSteerPyr(pyr, pyrFRs, whichPart)
% reconSteerBands: reconstruct an image from the subband transform.
%
% result = reconSteerBands(pyr,freqResps)
%    pyr = {highpass band, lowpass band, bandpass bands}
%    pyrFRs = {hiFR, loFR, FRs};
%       hiFR: frequency response of the highpass subband.
%       loFR: frequency response of the lowpass subband.
%       FRs: cell array that contains frequency responses of the bandpass subbands
%    whichPart: 'real' or 'imaginary'

hiFR = pyrFRs{1};
loFR = pyrFRs{2};
FRs = pyrFRs{3};
M = size(FRs,1);
N = size(FRs,2);
numLevels = size(FRs,3);
numOrientations = size(FRs,4);

hi = pyr{1};
lo = pyr{2};
bands = pyr{3};

result = zeros(M,N);

for level = 1:numLevels
  for orientation = 1:numOrientations
    if strcmp(whichPart, 'real')
      freqResp = real(FRs(:,:,level,orientation));
      band = real(bands(:,:,level,orientation));
    else
      freqResp = i * imag(FRs(:,:,level,orientation));
      band = imag(bands(:,:,level,orientation));
    end
    freqBand = fftshift(fft2(band));
    result = result + real(ifft2(fftshift(freqBand.*conj(freqResp))));
  end
end
%hiFreqBand = fftshift(fft2(hi));
%result = result + real(ifft2(fftshift(hiFreqBand.*conj(hiFR))));
%loFreqBand = fftshift(fft2(lo));
%result = result + real(ifft2(fftshift(loFreqBand.*conj(loFR))));



