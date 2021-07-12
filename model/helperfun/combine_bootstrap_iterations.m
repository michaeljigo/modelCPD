% Purpose:  Combine files across bootstrap iterations for ease of analysis.
%           Additionally, compute AIC and R2 measures for each iteration, as well as their bootstrapped confidence intervals
%
% By:       Michael Jigo
% Edited:   07.06.21
%
% Input:    inputs are all varargin variables and must be entered as a pair (e.g., 'sf_profile','match')
%
%           sf_profile        'match' (exo fit with narrow SF or endo fit with broad SF)  or 'mismatch' (exo and endo fit with the opposite SF profile)
%
%           model_variant     'main_model'      best-fitting model
%                             'minus_ori'       without cross-orientation suppression
%                             'minus_sf'        without cross-frequency suppression
%                             'minus_space'     without surround suppression
%                             'minus_context'   without all contextual modulation
%                             'minus_sum'       without spatial summation
%
%           attn_type         'exo', 'endo' or 'neutral'
% 
% Output:   out               structure containing concatenated data
%              fit_err        matrix of error of model fits
%              fit_params     matrix of all free parameters fit to the data
%              model_resp     structure containing all model-derived performance
%              exp_data       behavioral data that were fit
%              free_params    structure containing free parameters, formatted into easy-to-view structure
%              aic            AIC scores
%              bic            BIC scores
%              r2             R2 score
%              sim_params     parameters that controlled teh simulation
%              targfreq       target frequency for each texture stimulus

function out = combine_bootstrap_iterations(varargin)

%% Set default parameters
   in = {'sf_profile' ...        % 'match'=exo+narrow,endo+broad; 'mismatch'=exo+broad,endo+narrow
         'model_variant' ...     % 'main_model'; 'local_suppression'=narrow normalization pool; 'no_pooling'=no spatial summation 
         'attn_type'} ;          % 'exo'/'endo'/'neutral'
   val = {'match' ...
         'main_model' ...
         'neutral'};
   dp = parseOptionalInputs(in,val,varargin); 


%% Load bootstrap files
   % generate filename for saved bootstrap files
   if strcmp(dp.attn_type,'neutral')
      filename = sprintf('%s_%s',dp.attn_type,dp.model_variant);
   else
      filename = sprintf('%s_%s_%s',dp.attn_type,dp.sf_profile,dp.model_variant);
   end


   % get files from each bootsrap iteration
   datadir  = sprintf('../data/bootstrap_samples/%s',filename);
   files    = dir(sprintf('%s/*.mat',datadir));


   % loop through and concatenate bootstrap samples
   for f = 1:numel(files)
      load(sprintf('%s/%s',datadir,files(f).name));
      numiter = str2num(files(f).name(1:end-4)); % current iteration that's loaded

      % update parameter names
      newnames = {'attn_freq_max' 'attn_freq_slope' 'attn_bw' 'attn_amp_max' 'attn_spread'};
      oldnames = {'freq_ecc_intercept' 'freq_ecc_slope' 'bandwidth' 'amp_ecc_intercept' 'spread'};
      for n = 1:numel(newnames)
         idx = ismember(p.param_list,oldnames{n});
         if any(idx)
            p.param_list{idx} = newnames{n};
         end
      end
   
      if f==1
         % initialize matrices for each experiment
         all_fit_params = nan(numel(fit_params),p.bootstrap);
         all_fit_err    = nan(1,p.bootstrap);
         for e = 1:numel(p.exp_list)
            % these structures will hold model-derived d-prime and r2 for each iteration
            all_model_resp(e).dprime   = nan([size(exp_data(e).dprime) p.bootstrap]);
            r2.(p.exp_list{e})         = nan(1,p.bootstrap);
            for pa = 1:numel(p.param_list)
               % this structure stores the parameter values for each iteration
               free_params(e).(p.param_list{pa}) = nan(1,p.bootstrap);
            end
         end
      end


      % format fitted parameters into structures
      for e = 1:numel(p.exp_list)
         % initialize parameter structure for each experiment
            if f==1
               [stimdrive(e) supdrive(e) attn(e)] = init_parameters;
            end
         
         % unprune best-fitting parameters
            [stimdrive_fit supdrive_fit attn_fit] = prune_unprune_parameters(p,fit_params',stimdrive,supdrive,attn);
   
         % store parameters that were optimized
            for pa = 1:numel(p.param_list)
               switch p.param_list{pa}
                  case {'cg_max' 'cg_slope' 'freq_max' 'freq_slope' 'bw_max'}
                     free_params(e).(p.param_list{pa})(numiter) = stimdrive_fit(e).(p.param_list{pa});
                  case {'attn_freq_max' 'attn_freq_slope' 'attn_bw' 'attn_amp_max' 'attn_spread'}
                     free_params(e).(p.param_list{pa})(numiter) = attn_fit(e).(p.param_list{pa});
               end
            end
         
         % store model-derived performance
         all_model_resp(e).dprime(:,:,numiter)  = model_resp(e).dprime;

         % store # of observations in each experiment
         nobs(e)                                = numel(exp_data(e).dprime); 
   
         % SSE
         model_sse(numiter,e)                   = nansum((exp_data(e).dprime(:)-model_resp(e).dprime(:)).^2);
   
         % R2
         rss                                    = nansum((exp_data(e).dprime(:)-model_resp(e).dprime(:)).^2);
         tss                                    = nansum((exp_data(e).dprime(:)-mean(exp_data(e).dprime(:))).^2);
         r2.(p.exp_list{e})(numiter)            = 1-(rss./tss);

      end
   end

   % compute AIC and BIC scores
      % sum SSE across experiments
      nobs        = sum(nobs);          % # of observations that were fit
      nfree       = size(fit_params,1); % # of free parameters
      total_sse   = nansum(model_sse,2);
      tmp_aic     = nobs*log(total_sse./nobs)+(2*nfree);
      tmp_bic     = nobs*log(total_sse)-nobs*log(nobs)+(nfree)*log(nobs);

      % compute bootstrapped confidence interval for each model comparison metric
      [aic.ci aic.boot] = get_bootstrap_ci(tmp_aic,0.025,1e3);
      [bic.ci bic.boot] = get_bootstrap_ci(tmp_bic,0.025,1e3);


   % load the distribution of target frequencies
      if ismember(dp.attn_type,{'exo' 'endo'}) && ismember(dp.sf_profile,{'match'})
         ftargfile      = sprintf('../data/bootstrap_samples/%s_ftarg.mat',dp.attn_type);
         ftarg          = load(ftargfile);
         ftarg          = ftarg.ftarg;

         % separate out the fine and coarse textures based on viewing distance in the experiment
            switch dp.attn_type
               case 'exo'
                  fineidx                 = ismember(p.exp_list,{'yc98_exp1' 'yc08_cuesize1' 'ymc08_exp2' 'tc02_lvm'});
                  tmpcenter               = ftarg(fineidx,:);
                  targfreq.fine.center    = nanmedian(tmpcenter(:));
                  targfreq.fine.ci        = [nanmedian(tmpcenter(:))-std(tmpcenter(:)); nanmedian(tmpcenter(:))+std(tmpcenter(:))];
                  targfreq.idx.fine       = fineidx;
                  
                  coarseidx               = ismember(p.exp_list,{'yc98_exp2' 'clh06_baseline'});
                  tmpcenter               = ftarg(coarseidx,:);
                  targfreq.coarse.center  = nanmedian(tmpcenter(:));
                  targfreq.coarse.ci      = [nanmedian(tmpcenter(:))-std(tmpcenter(:)); nanmedian(tmpcenter(:))+std(tmpcenter(:))];
                  targfreq.idx.coarse     = coarseidx;
               case 'endo'
                  fineidx                 = ismember(p.exp_list,{'ymc08_exp3' 'ymc08_exp1' 'bc17_exp1'});
                  tmpcenter               = ftarg(fineidx,:);
                  targfreq.fine.center    = nanmedian(tmpcenter(:));
                  targfreq.fine.ci        = [nanmedian(tmpcenter(:))-std(tmpcenter(:)); nanmedian(tmpcenter(:))+std(tmpcenter(:))];
                  targfreq.idx.fine       = fineidx;

                  coarseidx               = ismember(p.exp_list,{'ymc08_exp4'});
                  tmpcenter               = ftarg(coarseidx,:);
                  targfreq.coarse.center  = nanmedian(tmpcenter(:));
                  targfreq.coarse.ci      = [nanmedian(tmpcenter(:))-std(tmpcenter(:)); nanmedian(tmpcenter(:))+std(tmpcenter(:))];
                  targfreq.idx.coarse     = coarseidx;

            end
            targfreq.all                  = ftarg;
      end


   % prepare variables for output
      out.fit_err       = model_sse;      % matrix of error of model fit
      out.fit_params    = all_fit_params; % matrix of all free parameters fit to the data
      out.model_resp    = all_model_resp; % structure containing all model-derived performance
      out.exp_data      = exp_data;       % behavioral data that were fit
      out.free_params   = free_params;    % structure containing free parameters, formatted into easy-to-view structure
      out.aic           = aic;            % AIC scores
      out.bic           = bic;            % BIC scores
      out.r2            = r2;             % R2 score
      out.sim_params    = p;              % parameters that controlled the simulation
      if exist('targfreq','var')
         out.targfreq   = targfreq;       % target frequency of texture stimuli
      end
