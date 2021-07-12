% Purpose:  Fit image-computable Attention model (imAmodel) to Neutral contrast thresholds and cueing effects from Experiment 2 in Jigo & Carrasco, 2020.
%
% By:       Michael Jigo
%           05.26.21

function fit_thresholds_exp2(varargin)
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
   8 ...                   % im_width
   1};                     % display_fit

p = parseOptionalInputs(in,val,varargin); 


%% Load experimental data
   load('../data/behavior/exp2.mat');

   % extract the data for the desired attention type
   attnidx = find(ismember(params.attntype,p.attn_type));

   % average the neutral conditions across experimental conditions
      % performance
      tmp_neut(:,:,:,1) = squeeze(data(1).performance(:,1,:,:));
      tmp_neut(:,:,:,2) = squeeze(data(2).performance(:,1,:,:));
      data(attnidx).performance(:,1,:,:) = mean(tmp_neut,4);

      % thresholds
      clear tmp_neut
      tmp_neut(:,:,:,1) = log(data(1).thresh);
      tmp_neut(:,:,:,2) = log(data(2).thresh);
      tmp_neut          = exp(mean(tmp_neut,4));
      data(attnidx).thresh = tmp_neut;
      
   % update data structure
   data = data(attnidx);


   % compute group-average thresholds
   avg_thresh = squeeze(exp(mean(log(data.thresh),1)));
      

   % use SFs, except highest
   freq_idx = 1:7;
   data.freq = data.freq(freq_idx);
   data.performance = data.performance(:,:,freq_idx,:);
   data.attn_effect = data.attn_effect(:,freq_idx,:);
   data.thresh = data.thresh(:,freq_idx,:);



%% Make and decompose gratings with distractors
   for f = 1:numel(data.freq)
      % Create targets and no-targets with distractors
         p.im_center = p.im_width/2; % define image center
         for e = 1:numel(data.ecc)
            % non-windowed gratings set to the group-average contrast threshold at each eccentricity
            [grat45(:,:,e),~,imap] = make_sinusoid2D('im_width',p.im_width,'window_width',0,'freq',data.freq(f),'amp',avg_thresh(f,e),'ori',45,'px_per_deg',p.px_per_deg);
            grat135(:,:,e) = make_sinusoid2D('im_width',p.im_width,'window_width',0,'freq',data.freq(f),'amp',avg_thresh(f,e),'ori',135,'px_per_deg',p.px_per_deg);
         
            % create windows at each eccentricity
            whorz = make_cosine_fun((imap.x(1,:)+p.im_center),data.ecc(e),data.grating_diameter,1);
            wvert = make_cosine_fun(imap.y(:,1),0,data.grating_diameter,1);
            window(:,:,e) = wvert*whorz;
         end

         % create non-target (135 deg tilted gratings at each eccentricity)
         grating.notarg(:,:,f) = (grat135(:,:,1).*window(:,:,1))+(grat135(:,:,2).*window(:,:,2));

         % create target at each eccentricity
         % 2 deg
         grating.targ(:,:,f,1) = (grat45(:,:,1).*window(:,:,1))+(grat135(:,:,2).*window(:,:,2));
         % 6 deg
         grating.targ(:,:,f,2) = (grat45(:,:,2).*window(:,:,2))+(grat135(:,:,1).*window(:,:,1));


      % Decompose images
      [~,energy.notarg(f)] = imAmodel(grating.notarg(:,:,f),'decompose_only',1,'px_per_deg',p.px_per_deg,'channel_freq_bw',0.5);
      for e = 1:numel(data.ecc)
         [~,energy.targ(f,e)] = imAmodel(grating.targ(:,:,f,e),'decompose_only',1,'px_per_deg',p.px_per_deg,'channel_freq_bw',0.5);
      end
   end


%% Initialize parameters
   % set up parameters for pruning
   p.exp_list = {'jc20'};
   p.param_list = {'cg_max' 'cg_slope' 'freq_max' 'freq_slope' 'bw_max' 'bw_slope' 'attn_freq_max' 'attn_freq_slope' 'attn_bw' 'attn_amp_max' 'attn_amp_slope' 'attn_spread'}; 
   p.full_index = ones(1,numel(p.param_list));
   [stimdrive, supdrive, attn, model_params, p.model_bnd, plaus_bnd, p.fixed_index] = prune_unprune_parameters(p);
   
   % set new initial parameters to speed up optimization
   %stimdrive.cg_max = 2.5356;
   %stimdrive.cg_slope = -0.0406;
   %stimdrive.freq_max = 1.1978;
   %stimdrive.freq_slope = -0.1789;
   %stimdrive.bw_max = 1.3132;

   %model_params(1) = 2.5356;
   %model_params(2) = -0.0406;
   %model_params(3) = 1.1978;
   %model_params(4) = -0.1789;
   %model_params(5) = 1.3132;

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
   options.MaxIter = 30;
      
   
   % run optimization
   objective = @(params)fit_thresholds_objective_exp2(stimdrive,supdrive,attn,data,grating,energy,p,params);
   [fit_params fit_err] = bads(objective,model_params',p.model_bnd(:,1)',p.model_bnd(:,2)',plaus_bnd(:,1)',plaus_bnd(:,2)',[],options); % BADS


%% Generate best-fitting values from best fits
[~, model] = fit_thresholds_objective_exp2(stimdrive,supdrive,attn,data,grating,energy,p,fit_params);

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
filename = sprintf('exp2_%s.mat',p.attn_type);
save([savedir,filename],'out');


%% display and save fit
if p.display_fit
   %display_exp2_fits('sf_profile',p.sf_profile,'model_variant',p.model_variant,'attn_type',p.attn_type);
end

