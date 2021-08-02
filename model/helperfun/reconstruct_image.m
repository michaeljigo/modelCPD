% Purpose:  Reconstruct an image after steerable pyramid decomposition.
%           Spatial frequencies and orientations can be selectively included in the reconstruction.

function [recon, nominalSF] = reconstruct_image(im,nori,freqBW,includeFreq,pxPerDeg)

% decompose
   nfreq = floor((log2(min(size(im)))-2)/freqBW);
   pyrFRs = makeQuadSteerFRs(size(im),nfreq,nori,freqBW); 
   pyr = buildQuadSteerPyr(im,pyrFRs);

% nominal spatial frequencies
   nfreq = floor((log2(min(size(im)))-2)/freqBW);
   freq_adjust    = mod(freqBW,1);
   if freq_adjust==0, freq_adjust = 1; end
   max_freq       = log2(pxPerDeg)-(1+freq_adjust); 
   channelFreq   = 2.^(max_freq-(0:freqBW:100));
   nominalSF      = round(channelFreq(1:nfreq)*10)./10;
  
% break function if includeFreq is empty
    if isempty(includeFreq)
        recon = nominalSF;
        return
    end

% reconstruct (based on reconQuadSteerPyr by David J. Heeger)
   recon = reconstruct(pyr,pyrFRs,'real',nominalSF,includeFreq);
   %recon = reconstruct(pyr,pyrFRs,'imaginary',nominalSF,includeFreq);

function result = reconstruct(pyr,pyrFRs,whichPart,nominalSF,includeFreq)
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

includeIdx = find(nominalSF>=min(includeFreq) & nominalSF<=max(includeFreq));

for level = includeIdx
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
