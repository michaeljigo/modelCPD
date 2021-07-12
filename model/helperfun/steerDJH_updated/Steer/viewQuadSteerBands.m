function viewBands(pyr, wait, rng)
% viewBands: view the subbands of the transform
%
% viewBands( pyr, pind, wait)
%    pyr: the image transform
%    wait: time to wait between displaying successive subband images.
%    rng: [min max] for displayImage (default 'auto1')
%
% DJH, 8/96
% updated 1/2019

if ~exist('wait')
  wait=1;
end

if ~exist('rng')
  rng = 'auto1';
end

hi = pyr{1};
lo = pyr{2};
bands = pyr{3};
numLevels = size(bands,3);
numOrientations = size(bands,4);

displayImage(hi);
pause(wait);
for level = 1:numLevels
  for orientation = 1:numOrientations;
    band = bands(:,:,level,orientation);
    displayImage(band,rng);
    pause(wait);
  end
end
displayImage(lo);

