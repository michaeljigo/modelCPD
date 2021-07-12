% Purpose:  This function will allow one to fit the model simultaneously to several experiments.
%           Fits will be configured in a manner identical to that described in the manuscript.
% 
% By:       Michael Jigo

function out = fit_exps(varargin)
addpath(genpath('../helperfun'));
addpath(genpath('../texture'));
addpath(genpath('~/apps'));

%% Set default parameters
in = {'sf_profile' ...   % 'narrow' or 'broad' or 'space_only'
   'model_variant' ...     % 'main_model', 'minus_sf', 'minus_ori', 'minus_space', 'minus_context', 'minus_sum'
   'attn_type' ...         % 'neutral' or 'involuntary' or 'voluntary'
   'display_fit'};         % 1=display fits; 0=don't display fits

val = {'narrow' ...        % sf_profile
   'main_model' ...        % model_variant
   'involuntary' ...       % attn_type
   1};                     % display_fit

p = parseOptionalInputs(in,val,varargin); 


%% Save path
savedir = sprintf('../data/fitted_parameters/%s/',p.model_variant);
if ~exist(savedir,'dir')
   mkdir(savedir);
end
filename = sprintf('%s_%s.mat',p.attn_type,p.sf_profile);


%% Choose expeirments to fit 
switch p.attn_type
   case 'neutral'
      exp_list = {'yc98_exp1' 'yc98_exp2' 'tc02' 'clh06' 'yc08' 'ymc08_exp2' 'ymc08_exp1' 'ymc08_exp3' 'ymc08_exp4' 'bc17'};
   case 'involuntary'
      exp_list = {'yc98_exp1' 'yc98_exp2' 'tc02' 'clh06' 'yc08' 'ymc08_exp2'};
   case 'voluntary'
      exp_list = {'ymc08_exp3' 'ymc08_exp4' 'ymc08_exp1' 'bc17'};
end
p.exp_list = exp_list;




%% Choose model parameters to be optimized
switch p.attn_type
   case 'neutral'
      param_list = {'cg_max' 'cg_slope' 'freq_max' 'freq_slope' 'bw_max'}; 
   case {'involuntary' 'voluntary'}
      if strcmp(p.sf_profile,'space_only')
         param_list = {'cg_max' 'cg_slope' 'freq_max' 'freq_slope' 'bw_max' 'attn_amp_max' 'attn_spread'};
      else
         param_list = {'cg_max' 'cg_slope' 'freq_max' 'freq_slope' 'bw_max' 'attn_freq_max' 'attn_freq_slope' 'attn_bw' 'attn_amp_max'};
      end
end
p.param_list = param_list;




%% Create indexing matrix
% this matrix determines which parameters are shared among experiments, and which are fixed
if numel(p.exp_list)>1
   for e = 1:numel(p.exp_list)
      param_index = ones(1,numel(param_list));
      switch p.exp_list{e}
         case {'yc98_exp1' 'yc98_exp2' 'tc02' 'clh06' 'yc08' 'ymc08_exp2'} % involuntary attention experiments
            % cg_max is free
            free_index = ismember(param_list,{'cg_max'});
            param_index(free_index) = e;

            % set spread to be free among experiments for space_only
            if strcmp(p.sf_profile,'space_only')
               free_index = ismember(param_list,{'attn_spread'});
               param_index(free_index) = e;
            end

            anchor = find(ismember(p.exp_list,'yc98_exp1')); % anchor attention parameters to be fixed to a single experiment
            if strcmp(p.attn_type,'neutral')
               fix_params = {'attn_freq_max' 'attn_freq_slope'};
               for a = 1:numel(fix_params)
                  pindex = ismember(param_list,fix_params{a});
                  param_index(pindex) = anchor;
               end
            end

         case {'ymc08_exp3' 'ymc08_exp4' 'ymc08_exp1' 'bc17'} % voluntary attention experiments
            % cg_max is free
            free_index = ismember(param_list,{'cg_max'});
            param_index(free_index) = e;

            % set spread to be free among experiments
            if strcmp(p.sf_profile,'space_only')
               free_index = ismember(param_list,{'attn_spread'});
               param_index(free_index) = e;
            end

            if strcmp(p.attn_type,'neutral')
               anchor = find(ismember(p.exp_list,'ymc08_exp3')); % anchor attention parameters to be fixed to a single experiment
               fix_params = {'freq_max'};
               for a = 1:numel(fix_params)
                  pindex = ismember(param_list,fix_params{a});
                  param_index(pindex) = anchor;
               end
            end
      end
      full_index(e,:) = param_index;
   end
else
   full_index = ones(1,numel(param_list)); % points to the experiment a given parameter will take its value from
end
% add neccessary values to the "p" structure
p.param_list = param_list;
p.full_index = full_index;




%% Initialize model parameters
% create random starting points
[stimdrive, supdrive, attn, model_params, p.model_bnd, plaus_bnd, p.fixed_index] = prune_unprune_parameters(p);

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




%% Load datasetes and images
for e = 1:numel(p.exp_list)
   tmp = load(sprintf('../data/behavior/%s.mat',p.exp_list{e}));
   data(e) = tmp.data;

   if strcmp(p.attn_type,'neutral')
      data(e).dprime = data(e).dprime(1,:);
      data(e).sem = data(e).sem(1,:);
   end
end



%% Decompose images to speed optimization
for e = 1:numel(p.exp_list)
   [~,energy.targ(e)] = imAmodel(data(e).stim.targ,'decompose_only',1);
   [~,energy.notarg(e)] = imAmodel(data(e).stim.notarg,'decompose_only',1);
end



%% Optimize
% BADS fitting options
options = bads('defaults');
options.MaxIter = 10;


% run optimization
objective = @(params)fit_exps_objective(stimdrive,supdrive,attn,data,energy,p,params);
[fit_params fit_err] = bads(objective,model_params',p.model_bnd(:,1)',p.model_bnd(:,2)',plaus_bnd(:,1)',plaus_bnd(:,2)',[],options); % BADS



%% Generate best-fitting values from best fits
[~, model_resp] = fit_exps_objective(stimdrive,supdrive,attn,data,energy,p,fit_params);

% rescale parameters
fit_params = fit_params'.*(p.unscaled_model_bnd(:,2)-p.unscaled_model_bnd(:,1))+p.unscaled_model_bnd(:,1);

% put parameters back in structure
[stimdrive,supdrive,attn] = prune_unprune_parameters(p,fit_params',stimdrive,supdrive,attn);


%% Create output variable
out.sse = fit_err;
out.stimdrive = stimdrive;
out.supdrive = supdrive;
out.attn = attn;
out.model_resp = model_resp;
out.p = p;
out.data = data;
out.fitparams = fit_params;

%% save
save([savedir,filename],'out');


%% display and save fit
if p.display_fit
   display_fits('sf_profile',p.sf_profile,'model_variant',p.model_variant,'attn_type',p.attn_type);
end
