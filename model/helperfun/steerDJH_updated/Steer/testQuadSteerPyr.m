%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Compute frequency responses of "steerable pyramid filters"

clear all;
figure(1);
colormap gray(256);

% Make and display frequency responses
dims = [64 64];
numOrientations = 6;
bandwidth = 1/2;
numLevels = maxLevel(dims,bandwidth)
pyrFRs = makeQuadSteerFRs(dims,numLevels,numOrientations,bandwidth);
viewQuadSteerBands(pyrFRs,1/4,'auto1');

% Check Tiling:
tile = zeros(dims);
tile = tile + abs(pyrFRs{1}).^2;
tile = tile + abs(pyrFRs{2}).^2;
for level = 1:numLevels
  for orientation = 1:numOrientations
    tile = tile + 0.5 * (abs(pyrFRs{3}(:,:,level,orientation))).^2;
  end
end
clf
imagesc(tile,[.999 1.001])
max(max(tile))
min(min(tile))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Test code for buildQuadSteerPyr and reconSteerPyr (no subsampling)

clear all;
al = double(imread('einstein.pgm'));
im = al(32:159,64:191);
displayImage(im);
colormap gray(256)

% Build and display pyramid
numOrientations = 6;
bandwidth = 1;
dims = size(im);
numLevels = maxLevel(dims,bandwidth);
pyrFRs = makeQuadSteerFRs(dims,numLevels,numOrientations,bandwidth);
pyr = buildQuadSteerPyr(im,pyrFRs);
viewQuadSteerBands(pyr,1/4,'auto1');

% View energies
for level = 1:numLevels
  for orientation = 1:numOrientations;
    band = pyr{3}(:,:,level,orientation);
    displayImage(abs(band).^2,'auto1');
    pause(1/4);
  end
end

% Reconstruct from real part
reconR = reconQuadSteerPyr(pyr, pyrFRs, 'real');
mse(double(im),reconR)
% Reconstruct from imaginary part
reconI = reconQuadSteerPyr(pyr, pyrFRs, 'imaginary');
mse(double(im),reconI)
displayImage(reconR + i*reconI);

% Check means (all but lowpass band should be zero)
mean2(im)
mean2(pyr{1})
mean2(pyr{2})
for level = 1:numLevels
  for orientation = 1:numOrientations
    mean2(pyr{3}(:,:,level,orientation))
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Normalized energy pyramid

clear all;
al = double(imread('einstein.pgm'));
im = al(32:159,64:191);

% N = 128;
% k = 16;
% c = 0.1;
% [xIm, yIm] = meshgrid([1:N],[1:N]);
% im = 2 * c * sin(2*pi*k/N*xIm) + 1;
% displayImage(im);

% Divide by mean and subtract 1 so that it is a "contrast"
% image.
im = im/mean2(im) - 1;

% Compute frequency responses
numOrientations = 6;
bandwidth = 1;
dims = size(im);
numLevels = maxLevel(dims,bandwidth);
pyrFRs = makeQuadSteerFRs(dims,numLevels,numOrientations,bandwidth);

% Build pyramid
pyr = buildQuadSteerPyr(im,pyrFRs);

% Normalize it
sigma = 0.1;
nEnergies = 1/4 * normEnergies(pyr,sigma);
%max(nEnergies(:))
mean(mean(sum(sum(nEnergies,4),3)))

% Display it
for level = 1:numLevels-1
  for orientation = 1:numOrientations;
    band = nEnergies(:,:,level,orientation);
    %['level: ', num2str(level), ' orientation: ', num2str(orientation), ' mean: ', num2str(mean2(band)), ' max: ', num2str(max2(band))]
    displayImage(band,[0 1/4]);
    pause(1/4);
  end
end

