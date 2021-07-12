function result = normEnergies(pyr,sigma)
% normEnergies: computes contrast normalized local energy responses from a
% quadrature subband transform.
%
% result = normEnergies(pyr,sigma)
%    pyr: the image transform returned by buildQuadSteerPyr
%       pyr = {highpass band, lowpass band, bandpass bands} 
%    sigma: semi-saturation constant for the contrast normalization.
%
% result: [N,M,numLevels,numOrientations] array
%
% DJH, 8/96

if ~exist('sigma')
  sigma=0;
end

hi = pyr{1};
lo = pyr{2};
bands = pyr{3};
M = size(bands,1);
N = size(bands,2);
numLevels = size(bands,3);
numOrientations = size(bands,4);

result = zeros([size(bands,1) size(bands,2) size(bands,3)-1 size(bands,4)]);

energies = real(bands).^2 + imag(bands).^2;
normalizers = sum(energies,4);
for level = 1:numLevels-1
  normalizers(:,:,level) = blur(normalizers(:,:,level),level);
end

% Compute hipass normalizer (spatial blurring here because there
% is no quadrature pair).
hiNormalizer = blur(hi.^2);

% Normalize each bandpass energy band, ignoring top level
for level = 1:numLevels-1
  if (level > 1) & (level < numLevels)
    denominator = normalizers(:,:,level-1) + normalizers(:,:,level) + normalizers(:,:,level+1) + sigma^2;
  elseif (level == 1)
    denominator = hiNormalizer + normalizers(:,:,level) + normalizers(:,:,level+1) + sigma^2;
  end
  for orientation = 1:numOrientations
    result(:,:,level,orientation) = energies(:,:,level,orientation) ./ denominator;
  end
end

