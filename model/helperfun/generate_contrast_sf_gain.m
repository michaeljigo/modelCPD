% Purpose:  Generate contrast and SF gain.
%
%           - Eccentricity changes at individual pixels based on distance from center of image. 
%           - Images will be centered at a given eccentricity, defined by the simulation parameters.
%           - SF gain are defined by a log-parabola whose parameters (amplitude, bandwidth, and peak frequency) will vary with eccentricity based on inputted parameters. 
%           - Contrast gain is defined by an exponential function across eccentricity
%
% Inputs:   channel        steerable pyramid subband information
%           stimdrive      structure of stimulus drive and contrast gain parameters
%           im_params      image parameters
%           ecc_im         eccentricity indices
%
% Outputs:  sf_gain        SF gain across eccentricity
%           contrast_gain  contrast gain across eccentricity
%
% By:       Michael Jigo

function [sf_gain contrast_gain] = generate_contrast_sf_gain(channel,stimdrive,im_params,ecc_im)

%% Generate weights at each pixel of eccentricity matrix
% evaluate SF gain at each eccentricity spanned by image (log parabola); maximum of function always 1
sf_gain = evaluate_log_parabola('amp_max',0,'amp_slope',0,'freq_max',stimdrive.freq_max,'freq_slope',stimdrive.freq_slope,'bw_max',stimdrive.bw_max,'bw_slope',stimdrive.bw_slope,...
   'ecc',ecc_im,'channel_freq',channel.freq,'channel_ori',channel.ori);

% evaluate contrast gain at each eccentricity spanned by image
[~, contrast_gain] = evaluate_log_parabola('amp_max',stimdrive.cg_max,'amp_slope',stimdrive.cg_slope,'ecc',ecc_im,'channel_freq',channel.freq,'channel_ori',channel.ori);
contrast_gain = 1./contrast_gain; % compute inverse such that small values correspond to high contrast sensitivity

% turn infs into nans
sf_gain(sf_gain==inf) = nan;
