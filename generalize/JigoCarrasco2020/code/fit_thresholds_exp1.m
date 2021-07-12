% Purpose:  Fit image-computable Attention model (imAmodel) to contrast thresholds from Experiment 1 in Jigo & Carrasco, 2020.
%
% By:       Michael Jigo
%           05.25.21

function fit_thresholds_exp1(varargin)
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
   5 ...                   % im_width
   1};                     % display_fit

p = parseOptionalInputs(in,val,varargin); 


%% Load experimental data
   load('../data/behavior/exp1.mat');

   % extract the data for the desired attention type
   attnidx = find(ismember(params.attntype,p.attn_type));
   data = data(attnidx);

   % use all SFs but the highest
   freq_idx = 1:5;
   data.freq = data.freq(freq_idx);
   data.crf.thresh = data.crf.thresh(:,:,freq_idx,:);


%% Make and decompose gratings that match stimuli in the experiment
   p.contrasts = logspace(-3,0,7);
   for f = 1:numel(data.freq)
      for c = 1:numel(p.contrasts)
         % create target and no-target images
         tmp45 = make_sinusoid2D('im_width',p.im_width,'window_width',data.grating_diameter,'freq',data.freq(f),'amp',p.contrasts(c),'ori',45,'px_per_deg',p.px_per_deg);
         tmp135 = make_sinusoid2D('im_width',p.im_width,'window_width',data.grating_diameter,'freq',data.freq(f),'amp',p.contrasts(c),'ori',135,'px_per_deg',p.px_per_deg);
         grating.targ(:,:,f,c) = tmp45;
         grating.notarg(:,:,f,c) = tmp135;

         % decompose images
         [~,energy.targ(f,c)] = imAmodel(tmp45,'decompose_only',1,'px_per_deg',p.px_per_deg);
         [~,energy.notarg(f,c)] = imAmodel(tmp135,'decompose_only',1,'px_per_deg',p.px_per_deg);
      end
   end


%% Initialize parameters
   % set up parameters for pruning
   p.exp_list = {'jc20'};
   p.param_list = {'cg_max' 'cg_slope' 'freq_max' 'freq_slope' 'bw_max' 'attn_freq_max' 'attn_freq_slope' 'attn_bw' 'attn_amp_max' 'attn_spread'};
   p.full_index = ones(1,numel(p.param_list));
   [stimdrive, supdrive, attn, model_params, p.model_bnd, plaus_bnd, p.fixed_index] = prune_unprune_parameters(p);

   % set new initial parameters to speed up optimization
   switch p.attn_type
      case 'exo'
         %stimdrive.cg_max = 2.5008;
         %stimdrive.cg_slope = -0.0337;
         %stimdrive.freq_max = 1.0981;
         %stimdrive.freq_slope = -0.2029;
         %stimdrive.bw_max = 1.3965;

         %attn.attn_freq_max = 1.7030;
         %attn.attn_freq_slope = -0.1;
         %attn.attn_bw = 2.2881;
         %attn.attn_amp_max = 1.4571;
         %attn.attn_spread = 5.3725;

         %model_params(1) = 2.5008;
         %model_params(2) = -0.0337;
         %model_params(3) = 1.0981;
         %model_params(4) = -0.2029;
         %model_params(5) = 1.3965;
         %model_params(6) = 1.7030;
         %model_params(7) = -0.0312;
         %model_params(8) = 2.2881;
         %model_params(9) = 1.4571;
         %model_params(10) = 5.3725;
      case 'endo'
         %stimdrive.cg_max = 2.4764;
         %stimdrive.cg_slope = -0.0358;
         %stimdrive.freq_max = 1.0867;
         %stimdrive.freq_slope = -0.1925;
         %stimdrive.bw_max = 1.4663;

         %attn.attn_freq_max = 2.8050;
         %attn.attn_freq_slope = -0.2;
         %attn.attn_bw = 1.6798;
         %attn.attn_amp_max = 1.1763;
         %attn.attn_spread = 5.2111;

         %model_params(1) = 2.4764;
         %model_params(2) = -0.0358;
         %model_params(3) = 1.0867;
         %model_params(4) = -0.1925;
         %model_params(5) = 1.4663;
         %model_params(6) = 2.805;
         %model_params(7) = -0.2;
         %model_params(8) = 1.6798;
         %model_params(9) = 1.1763;
         %model_params(10) = 5.2111;
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
   options.MaxIter = 7;
      
   
   % run optimization
   objective = @(params)fit_thresholds_objective_exp1(stimdrive,supdrive,attn,data,grating,energy,p,params);
   [fit_params fit_err] = bads(objective,model_params',p.model_bnd(:,1)',p.model_bnd(:,2)',plaus_bnd(:,1)',plaus_bnd(:,2)',[],options); % BADS


%% Generate best-fitting values from best fits
[~, modelcs] = fit_thresholds_objective_exp1(stimdrive,supdrive,attn,data,grating,energy,p,fit_params);

% rescale parameters
fit_params = fit_params'.*(p.unscaled_model_bnd(:,2)-p.unscaled_model_bnd(:,1))+p.unscaled_model_bnd(:,1);

% put parameters back in structure
[stimdrive,supdrive,attn] = prune_unprune_parameters(p,fit_params',stimdrive,supdrive,attn);


%% Create output variable
out.sse = fit_err;
out.stimdrive = stimdrive;
out.supdrive = supdrive;
out.attn = attn;
out.modelcs = modelcs;
out.p = p;
out.data = data;
out.fitparams = fit_params;

%% save
savedir = '../data/fitted_parameters/';
if ~exist(savedir,'dir')
   mkdir(savedir);
end
filename = sprintf('exp1_%s_restricted_within_texture.mat',p.attn_type);
save([savedir,filename],'out');


%% display and save fit
if p.display_fit
   %display_exp1_fits('sf_profile',p.sf_profile,'model_variant',p.model_variant,'attn_type',p.attn_type);
end
