% Purpose:  Generate image defining eccentricity at each pixel in image. Image is always centered on the nominal eccentricity.
%
% Input:    channel     steerable pyramid channel information
%           im_params   parameters of input image
%           sim_params  simulation parameters (i.e., eccentricities to be tested)
%
% Output:   ecc_im      eccentricity images
%           eccmap      extra information regarding components of eccentricity images
%
% By:       Michael Jigo

function [ecc_im eccmap] = generate_eccentricity_image(channel,im_params,sim_params)

%% Compute image size
if numel(im_params.im_size)==1
   % input is a square
   im_params.im_size = repmat(im_params.im_size,1,2);
end
im_size_pix = round(im_params.px_per_deg*im_params.im_size);
im_halfwidth = im_params.im_size./2;


%% Generate eccentricity matrix for each center eccentricity
% pre-allocate eccentricity matrix
ecc_im = zeros([numel(sim_params.ecc) 1 1 im_size_pix(1) im_size_pix(2)]);

% generate matrix for each eccentricity
width = linspace(-im_halfwidth(2),im_halfwidth(2),im_size_pix(2));
height = linspace(-im_halfwidth(1),im_halfwidth(1),im_size_pix(1));
[x,y] = meshgrid(width,height); %y = flipud(y);
% convert to polar coordinates
[ph,r] = cart2pol(x,y);
for e = 1:numel(sim_params.ecc)
   %ecc_im(e,1,1,:,:) = r+sim_params.ecc(e);
   ecc_im(e,1,1,:,:) = sqrt((x+sim_params.ecc(e)).^2+y.^2);
   allr(:,:,e) = r+sim_params.ecc(e);
end
% duplicate eccentricity image in each frequency
ecc_im = repmat(ecc_im,[1 numel(channel.ori) numel(channel.freq) 1 1]);

% create structure for mapping pixels to degrees
eccmap.x = x;
eccmap.y = y;
eccmap.r = allr;
eccmap.ph = ph;
