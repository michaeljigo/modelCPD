% Purpose:  Take input matrix of pixel intensitities and decompose it into energy at individual SFs and orientation subbands using steerable pyramid.
%           Each channel is the filter response centered at a given spatial frequency and orientation.
%           
%           This function calls Eero Simoncelli's steerable pyramid (Simoncelli, Freeman, Adelson & Heeger, 1992, Shiftable multiscale transforms) 
%           that has been updated by David J. Heeger such that spatial frequency subbands are not subsampled.
%
% Input:    im             2D matrix 
%           px_per_deg     spatial sampling rate of filters
%
% Output:   energy         4D matrix of energy computed at several spatial frequency and orientation channels
%           channel        information about the filter parameters
%
% By:       Michael Jigo

function [energy channel] = decompose_image(im,px_per_deg,varargin)

%% Set pyramid parameters
in = {'n_freq' 'n_ori' 'bw'};
val = {1e3 6 1};
channel = parseOptionalInputs(in,val,varargin);


%% Load image
if ischar(im)
   % input is path to file; load file
   im = pgmRead(im);
else
   % input is a matrix with corresponding pixel intensities
end


%% Compute maximum # of bands that can be made given image size
max_nfreq = floor((log2(min(size(im)))-2)/channel.bw);
channel.n_freq = max_nfreq;


%% Compute center frequency and orientation of each "channel"
% orientations
channel.ori = rad2deg(pi*(0:channel.n_ori-1)./channel.n_ori);

% frequencies
freq_adjust = mod(channel.bw,1);
if freq_adjust==0, freq_adjust = 1; end
max_freq = log2(px_per_deg)-(1+freq_adjust); % maximum possible frequency, given nyquist limit: 2 = 1 (from initial subsample) + 1 (from nyquist limit of subsampled image)
channel.freq = 2.^(max_freq-(0:channel.bw:100));
channel.freq = channel.freq(1:channel.n_freq);


%% Decompose and compute energy
% decompose
pyrFRs = makeQuadSteerFRs(size(im),channel.n_freq,channel.n_ori,channel.bw); % BANDWIDTHS ARE FWHM (full width at half maximum)
pyr = buildQuadSteerPyr(im,pyrFRs);

% compute energy
pyr = cellfun(@(x) (real(x).^2+imag(x).^2),pyr,'UniformOutput',0); % energy

% create 4D energy matrix
energy = zeros([channel.n_ori channel.n_freq size(im)]);
for f = 1:channel.n_freq
   for o = 1:channel.n_ori
      energy(o,f,:,:) = squeeze(pyr{3}(:,:,f,o));
   end
end
