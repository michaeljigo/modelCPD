% Purpose:        Generate model output using imAmodel (IMage-computable Attention model). 
%                 Set model parameters (or use defaults) and evaluate an image at a range of eccentricities (inputted or default values).
% By:             Michael Jigo
% Last updated:   04/22/21
%
% Input:          im          2D image
%                 varargin    optional inputs, see list below 
%
% Output:         popresp     population response model output (4D matrix: [ecc ori sf x y])
%                 params      model evaluation parameters

function [popresp params] = imAmodel(im,varargin)
addpath(genpath('./helperFun'));

%%%% Set default parameters 
in = {'px_per_deg' ...        % pixels per degree 
      'im_size' ...           % input image size (in degrees of visual angle)
      'ecc' ...               % eccentricity(ies) to simulate (scalar or vector)
      'cg_max' ...            % contrast gain maximum
      'cg_slope' ...          % contrast gain slope
      'freq_max' ...          % intrinsic spatial frequency maximum
      'freq_slope' ...        % intrinsic spatial frequency slope
      'bw_max' ...            % intrinsic bandwidth maximum
      'bw_slope' ...          % intrinsic bandwidth slope along eccentricity
      'attn_freq_max' ...     % attention spatial frequency maximum
      'attn_freq_slope' ...   % attention spatial frequency slope
      'attn_bw' ...           % spatial frequency bandwidth of attention
      'attn_amp_max' ...      % amplitude amplitude maximum
      'attn_amp_slope' ...    % attention amplitude slope
      'attn_spread' ...       % spatial spread of attention
      'sf_profile' ...        % 'narrow', 'broad' or 'space_only'
      'spatial_profile' ...   % function defining spatial spread of attention 'center' or 'center_surround'
      'use_attn' ...          % evaluate attention parameters
      'verbose' ...           % display output 
      'channel_freq_bw' ...   % bandwidth (FWHM in octaves) of pyramid spatial frequency channels
      'channel_n_ori' ...     % # of orientation channels in pyramid (FWHM=360/channel_n_ori)
      'stimdrive' ...         % structure of stimulus drive and contrast gain parameters
      'supdrive' ...          % structure of suppressive drive (suppressive pools only) parameters
      'attn' ...              % structure of attention parameters
      'preprocess_image' ...  % 0=do nothing; 1=change range of image, then window 
      'energy' ...            % allows energy to be input, so that pyramid decomposition is skipped
      'decompose_only' ...    % only decompose input image then return function
      'attended_ecc' ...      % location at which attention is centered (nan=follow the stimulus ecc)
      'model_variant'};       % model variant that will be fit: 'main_model', 'minus_sf', 'minus_ori', 'minus_space', 'minus_context', 'minus_sum'

val = {32 ...                 % px_per_deg
       [] ...                 % im_size
       0:10 ...               % ecc
       1.9  ...               % cg_max
       -0.05 ...              % cg_slope
       2.1 ...                % freq_max
       -0.7 ...               % freq_slope
       1.6 ...                % bw_max
       0 ...                  % bw_slope
       2.3 ...                % attn_freq_max
       -0.2 ...               % attn_freq_slope
       2.9 ...                % attn_bandwidth
       4.3 ...                % attn_amp_max
       0 ...                  % attn_amp_slope
       4 ...                  % attn_spread
       'narrow' ...           % sf_profile
       'center' ...           % spatial_profile
       0 ...                  % use_attn
       0 ...                  % verbose
       1 ...                  % channel_freq_bw
       6 ...                  % channel_n_ori
       [] ...                 % stimdrive
       [] ...                 % supdrive
       [] ...                 % attn
       0 ...                  % preprocess_image
       [] ...                 % energy
       0 ...                  % decompose_only
       nan ...                % attended_ecc
       'main_model'};         % model_variant

%% Parse varargin inputs
p = parseOptionalInputs(in,val,varargin);
if isempty(p.im_size)
   p.im_size = size(im)./p.px_per_deg;
end



%%%% Preprocessing and steerable pyramid decomposition
%% Preprocess image
if p.preprocess_image
   % change into full-contrast image
   im = double(im);
   im = change_range(im,0,1,min(im(:)),max(im(:)));

   % window image with a cosine function to eliminate edge effects
   spacex = size(im,1); spacey = size(im,2);
   spacex = linspace(-spacex/2,spacex/2,spacex);
   spacey = linspace(-spacey/2,spacey/2,spacey);
   [winx winax] = make_cosine_fun(spacex,0,max(spacex)*(3/2)*0.7,1); % 1D across x
   [winy winay] = make_cosine_fun(spacey,0,max(spacey)*(3/2)*0.7,1); % 1D across y
   winx = winx+winax;
   winy = winy+winay;
   win = winx'*winy; % make 2D

   % do windowing
   im = im.*win; 
end


%% Decompose input image and compute energy
if isempty(p.energy)
   % decompose input image
   [energy.im, energy.channel] = decompose_image(im,p.px_per_deg,'bw',p.channel_freq_bw,'n_ori',p.channel_n_ori);
   if p.verbose; fprintf('Image decomposed.\n'); end
else
   % a decomposed image was input, so...
   % check that channel parameters match desired values
   ori_check = isequal(p.energy.channel.n_ori,p.channel_n_ori);
   freq_check = isequal(p.energy.channel.bw,p.channel_freq_bw);
   if ori_check==1 && freq_check==1
      energy = p.energy;
   else
      error('Steerable pyramid parameters do not match.');
   end
end


%% Return function if preprocessing and decomposing is all that's desired
if p.decompose_only
   popresp = [];
   params.energy = energy;
   return
end



%%%% Prepare model parameters
%% Initialize parameters
[stimdrive, supdrive, attn] = init_parameters('cg_max',p.cg_max,'cg_slope',p.cg_slope,'freq_max',p.freq_max,'freq_slope',p.freq_slope,'bw_max',p.bw_max,'bw_slope',p.bw_slope,'attn_freq_max',p.attn_freq_max,'attn_freq_slope',p.attn_freq_slope,'attn_bw',p.attn_bw,'attn_amp_max',p.attn_amp_max,'attn_amp_slope',p.attn_amp_slope,'attn_spread',p.attn_spread,'spatial_profile',p.spatial_profile); 

% If parameter structures were input, use those instead
if ~isempty(p.stimdrive)
   stimdrive = p.stimdrive;
end
if ~isempty(p.supdrive)
   supdrive = p.supdrive;
end
if ~isempty(p.attn)
   attn = p.attn;
end


%% Create eccentricity image(s)
[ecc_im eccmap] = generate_eccentricity_image(energy.channel,p,p);


%% Compute contrast gain & SF gain
% here, contrast gain values are transferred to the suppressive drive structure; will improve in future versions
[stimdrive.sf_gain supdrive.contrast_gain] = generate_contrast_sf_gain(energy.channel,stimdrive,p,ecc_im);


%% Create normalization pools
supdrive = generate_normalization_pools(energy.channel,p,supdrive,ecc_im);


%% Generate attention field matrices
if p.use_attn
   attn = generate_attention_modulation(attn,energy.channel,p,ecc_im,eccmap);
end
if p.verbose; fprintf('Parameters initialized.\n'); end



%%%% Evaluate model
t = tic;
if p.verbose; fprintf('Evaluating eccentricity....\n'); end
popresp = nan([numel(p.ecc) size(energy.im)]);
for e = 1:numel(p.ecc)
   % contrast gain for current eccentricity
   ecc_exc = stimdrive;
   ecc_exc.sf_gain = squeeze(ecc_exc.sf_gain(e,:,:,:,:));

   % sf gain for current eccentricity
   ecc_inh = supdrive;
   ecc_inh.contrast_gain = squeeze(ecc_inh.contrast_gain(e,:,:,:,:));

   % atteniton modulation for current eccentricity
   if p.use_attn
      attn_modulation = squeeze(attn.modulation(e,:,:,:,:));
   else
      attn_modulation = 1;
   end

   % compute normalized image
   popresp(e,:,:,:,:) = compute_model_response(energy.im,ecc_exc,ecc_inh,attn_modulation,p.model_variant);
   if p.verbose; fprintf('%i/%i\n'); end
end
t = toc(t);
if p.verbose; fprintf('Done. Took %.1f seconds\n',t); end


%% Package output
params.stimdrive  = stimdrive;
params.supdrive   = supdrive;
params.attn       = attn;
params.channel    = energy.channel;
params.energy     = energy;
params.eccmap     = eccmap;
params.options    = p;


