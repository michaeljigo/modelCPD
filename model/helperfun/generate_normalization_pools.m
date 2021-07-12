% Purpose:  Specify the kernels that will define the extent of normalization pools in spatial (x,y), orientation, and spatial frequency dimensions.
%
%           - Spatial kernel will be a cosine function whose width will be a scaled copy of the width of the filter's frequency.
%           - Orientation and frequency kernels will be defined by logical vectors indicating which orientation/frequency channels will be included in their respective normalization pools.
%
% Inputs:   channel        steerable pyramid subband information
%           supdrive       structure of suppressive drive parameters
%           im_params      image parameters
%           ecc_im         eccentricity indices
%
% Output:   supdrive       suppressive drive structure updated with normalization pools
%
% By:       Michael Jigo

function supdrive = generate_normalization_pools(channel,im_params,supdrive,ecc_im)

%% Initialize fields
supdrive.space_kernels = [];
supdrive.ori_pool = [];
supdrive.freq_pool = [];

%% Space (i.e., surround suppression)
% For the spatial dimension, suppression will be computed as the 
% convolution of raised cosine kernels whose widths vary with SF

% determine number of degrees in a cycle for each filter frequency
deg_per_cycle = 1./channel.freq;

% kernel widths will be a scaled copy of the cycle widths
space_pool_widths = deg_per_cycle*supdrive.sup_space;

% initialize matrix
if numel(im_params.im_size)==1
   % input is a square
   im_params.im_size = repmat(im_params.im_size,1,2);
end

% generate raised-cosine kernels for each frequency
im_size_pix = min(im_params.px_per_deg*im_params.im_size);
im_halfwidth = min(im_params.im_size)/2;
space = linspace(-im_halfwidth,im_halfwidth,im_size_pix);

% pre-allocate matrix
for f = 1:numel(channel.freq)
   kern = make_cosine_fun(space,0,space_pool_widths(f),1);
   supdrive.space_kernels(f,:,:) = kern'*kern;
end


%% Orientation
% For orientaiton, suppression will be specified by logical vectors representing which orientation
% channels will be included in the suppression pool

% pre-allocate matrix
supdrive.ori_pool = zeros(numel(channel.ori),numel(channel.ori));
for o = 1:numel(channel.ori)
   % orientation is a circular dimension, so wrap the pool around at boundaries
   bounds = [-supdrive.sup_ori supdrive.sup_ori]/2+channel.ori(o);
   bounds = mod(bounds,180);
   if range(bounds)<1e-3
      bounds = [0 180];
   end
   % define orientation normalization pools
   if diff(bounds)<0
      supdrive.ori_pool(o,:) = channel.ori>=bounds(1) | channel.ori<=bounds(2);
   else
      supdrive.ori_pool(o,:) = channel.ori>=bounds(1) & channel.ori<=bounds(2);
   end
end


%% Frequency
% For spatial frequency, suppression will be specified by logical vectors representing which SF bands
% will be included in the suppression pool

% pre-allocate matrix
supdrive.freq_pool = zeros(numel(channel.freq),numel(channel.freq));

% fill matrix with pool indices
for f = 1:numel(channel.freq)
   bounds = 2.^([-supdrive.sup_freq supdrive.sup_freq]/2+log2(channel.freq(f)));
   supdrive.freq_pool(f,:) = channel.freq>=bounds(1) & channel.freq<=bounds(2);
end
