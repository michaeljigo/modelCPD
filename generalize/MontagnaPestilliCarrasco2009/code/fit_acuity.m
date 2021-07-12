% Purpose:  Fit image-computable Attention model (imAmodel) to acuity thresholds from Montagna, Pestilli & Carrasco, 2009.
%
% By:       Michael Jigo
%           05.25.21

function fit_acuity(varargin)
addpath(genpath('~/apps'));
addpath(genpath('../../../../modelCPD_v4'));

%% Set default parameters
in = {'sf_profile' ...     % 'narrow' or 'broad' or 'space_only'
   'model_variant' ...     % 'main_model', 'minus_sf', 'minus_ori', 'minus_space', 'minus_context', 'minus_sum'
   'attn_type' ...         % 'exo' or 'endo'
   'px_per_deg' ...        % pixels per degree of images
   'im_width' ...          % side length of image in degrees
   'display_fit'};         % 1=display fits; 0=don't display fits

val = {'narrow' ...        % sf_profile
   'main_model' ...        % model_variant
   'exo' ...               % attn_type
   32 ...                  % px_per_deg
   3 ...                   % im_width
   1};                     % display_fit

p = parseOptionalInputs(in,val,varargin); 


%% Load experimental data
   load('../data/behavior/acuity.mat');

   % extract the data for the desired attention type
   attnidx = find(ismember(params.attntype,p.attn_type));
   data = data(attnidx);


%% Make and decompose Landolt squares
   p.gap_size = linspace(0,30,10); % in arc minutes
   for g = 1:numel(p.gap_size)
      % create target and no-target images
      landolt.targ(:,:,g) = make_landolt_square('gap_size',p.gap_size(g),'im_size',p.im_width,'px_per_deg',p.px_per_deg);
      landolt.notarg(:,:,g) = imrotate(landolt.targ(:,:,g),180); % no-target is flipped variant of target

      % decompose images
      [~,energy.targ(g)] = imAmodel(landolt.targ(:,:,g),'decompose_only',1,'px_per_deg',p.px_per_deg);
      [~,energy.notarg(g)] = imAmodel(landolt.notarg(:,:,g),'decompose_only',1,'px_per_deg',p.px_per_deg);
   end


%% Initialize parameters
   % set up parameters for pruning
   p.exp_list = {'jc20'};
   p.param_list = {'cg_max' 'cg_slope' 'freq_max' 'freq_slope' 'bw_max' 'attn_freq_max' 'attn_freq_slope' 'attn_bw' 'attn_amp_max' 'attn_spread'};
   p.full_index = ones(1,numel(p.param_list));
   [stimdrive, supdrive, attn, model_params, p.model_bnd, plaus_bnd, p.fixed_index] = prune_unprune_parameters(p);
   
   % set new initial parameters to speed up optimization
   switch p.attn_type
      case 'endo'
         stimdrive.cg_max = 1.5;
         stimdrive.cg_slope = -0.0912;
         stimdrive.freq_max = 3;
         stimdrive.freq_slope = -0.9;
         stimdrive.bw_max = 1.1;

         attn.attn_freq_max = 2;
         attn.attn_freq_slope = -0.3;
         attn.attn_bw = 3.5;
         attn.attn_amp_max = 15;
         attn.attn_spread = 0.6;


         model_params(1) = 1.5;
         model_params(2) = -0.0912;
         model_params(3) = 3;
         model_params(4) = -0.9;
         model_params(5) = 1.1;
         model_params(6) = 2;
         model_params(7) = -0.3;
         model_params(8) = 3.5;
         model_params(9) = 15;
         model_params(10) = 0.6;
      case 'exo'
         stimdrive.cg_max = 2;
         stimdrive.cg_slope = -0.07;
         stimdrive.freq_max = 1.8118;
         stimdrive.freq_slope = -0.7777;
         stimdrive.bw_max = 1.4001;

         attn.attn_freq_max = 2.8;
         attn.attn_freq_slope = -0.03;
         attn.attn_bw = 1.5;
         attn.attn_amp_max = 10;
         attn.attn_spread = 0.6003;

         model_params(1) = 2;
         model_params(2) = -0.07;
         model_params(3) = 1.8118;
         model_params(4) = -0.7777;
         model_params(5) = 1.4001;
         model_params(6) = 2.8;
         model_params(7) = -0.04;
         model_params(8) = 1.5;
         model_params(9) = 8;
         model_params(10) = 0.6003;
   end

   % scale model parameters to range 0-1
   model_params = (model_params-p.model_bnd(:,1))./(p.model_bnd(:,2)-p.model_bnd(:,1)); % random starting points
   p.unscaled_model_bnd = p.model_bnd;

   % scale plausible bounds, for use with BADS
   for col = 1:2
      plaus_bnd(:,col) = (plaus_bnd(:,col)-p.model_bnd(:,1))./(p.model_bnd(:,2)-p.model_bnd(:,1));
   end
   
   % change hard bounds to reflect scaling
   p.model_bnd = repmat([0 1],[size(p.model_bnd,1),1]); 
   
   % change nans in plaus bnd to [0 1]
   for r = 1:size(plaus_bnd,1)
      if any(isnan(plaus_bnd(r,:)))
         plaus_bnd(r,:) = [0 1];
      end
   end
   

%% Optimize
   % BADS fitting options
   options = bads('defaults');
   options.TolMesh = 1e-4;
   options.TolFun  = 1e-2;
   options.MaxIter = 4;
      
   
   % run optimization
   objective = @(params)fit_acuity_objective(stimdrive,supdrive,attn,data,landolt,energy,p,params);
   %[fit_params fit_err] = bads(objective,model_params',p.model_bnd(:,1)',p.model_bnd(:,2)',plaus_bnd(:,1)',plaus_bnd(:,2)',[],options); % BADS


%% Generate best-fitting values from best fits
fit_params = model_params';
[fit_err, model] = fit_acuity_objective(stimdrive,supdrive,attn,data,landolt,energy,p,fit_params);

% rescale parameters
fit_params = fit_params'.*(p.unscaled_model_bnd(:,2)-p.unscaled_model_bnd(:,1))+p.unscaled_model_bnd(:,1);

% put parameters back in structure
[stimdrive,supdrive,attn] = prune_unprune_parameters(p,fit_params',stimdrive,supdrive,attn);


%% Create output variable
out.sse = fit_err;
out.stimdrive = stimdrive;
out.supdrive = supdrive;
out.attn = attn;
out.model = model;
out.p = p;
out.data = data;
out.fitparams = fit_params;

%% save
savedir = '../data/fitted_parameters/';
if ~exist(savedir,'dir')
   mkdir(savedir);
end
keyboard
filename = sprintf('%s_restricted_within_texture.mat',p.attn_type);
save([savedir,filename],'out');


%% display and save fit
if p.display_fit
   %display_acuity_fits('sf_profile',p.sf_profile,'model_variant',p.model_variant,'attn_type',p.attn_type);
end
