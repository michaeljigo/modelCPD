% Purpose:  Objective function for fit_exps

function [cost model_resp] = fit_exps_objective(stimdrive,supdrive,attn,data,energy,p,params)

%% Re-structure parameters
% re-scale parameters
params = params'.*(p.unscaled_model_bnd(:,2)-p.unscaled_model_bnd(:,1))+p.unscaled_model_bnd(:,1);
[~,~,~,param_vals] = prune_unprune_parameters(p,params,stimdrive,supdrive,attn);

for e = 1:size(param_vals,1)
   tmp_stimdrive{e} = stimdrive(e);
   tmp_supdrive{e} = supdrive(e);
   if ~isempty(attn)
      tmp_attn{e} = attn(e);
   else
      tmp_attn{e} = attn;
   end

   for ii = 1:numel(p.param_list)
      switch p.param_list{ii}
         case {'cg_max' 'cg_slope' 'freq_max' 'freq_slope' 'bw_max'}
            if ~isnan(param_vals(e,ii))
               tmp_stimdrive{e}.(p.param_list{ii}) = param_vals(e,ii);
            else
               tmp_stimdrive{e}.(p.param_list{ii}) = param_vals(p.full_index(e,ii),ii);
            end 
         case {'attn_freq_max' 'attn_freq_slope' 'attn_bw' 'attn_amp_max' 'attn_spread'}
            if ~isempty(tmp_attn)
               if ~isnan(param_vals(e,ii))
                  tmp_attn{e}.(p.param_list{ii}) = param_vals(e,ii);
               else
                  tmp_attn{e}.(p.param_list{ii}) = param_vals(p.full_index(e,ii),ii);
               end 
            end
      end
   end
end


%%%% Compute model response to images
parfor exp_num = 1:numel(data)

   %% Generate parameters
   % extract parameters for this experiment
   current_stimdrive = tmp_stimdrive{exp_num};
   current_supdrive = tmp_supdrive{exp_num};
   if ~isempty(attn)
      current_attn = tmp_attn{exp_num};
   else
      current_attn = [];
   end

   %% Apply parameters to evaluate model
   dprime = nan(size(data(exp_num).dprime,1),numel(data(exp_num).ecc));
   for c = 1:size(data(exp_num).dprime,1)
      % evaluate model
      use_attn = c-1;
      %targresp = run_model(data(exp_num).stim.targ,'energy',energy.targ(exp_num).energy,'ecc',data(exp_num).ecc,'stimdrive',current_stimdrive,'supdrive',current_supdrive,'attn',current_attn,'use_attn',use_attn,'gain_profile',p.gain_profile,'model_variant',p.model_variant);
      %notargresp = run_model(data(exp_num).stim.notarg,'energy',energy.notarg(exp_num).energy,'ecc',data(exp_num).ecc,'stimdrive',current_stimdrive,'supdrive',current_supdrive,'attn',current_attn,'use_attn',use_attn,'gain_profile',p.gain_profile,'model_variant',p.model_variant);
      targresp = imAmodel(data(exp_num).stim.targ,'energy',energy.targ(exp_num).energy,'ecc',data(exp_num).ecc,'stimdrive',current_stimdrive,'supdrive',current_supdrive,'attn',current_attn,'use_attn',use_attn,'sf_profile',p.sf_profile,'model_variant',p.model_variant);
      notargresp = imAmodel(data(exp_num).stim.notarg,'energy',energy.notarg(exp_num).energy,'ecc',data(exp_num).ecc,'stimdrive',current_stimdrive,'supdrive',current_supdrive,'attn',current_attn,'use_attn',use_attn,'sf_profile',p.sf_profile,'model_variant',p.model_variant);
      
      % compute dprime (Eucledian norm)
      dprime(c,:) = sqrt(nansum((targresp(:,:)-notargresp(:,:)).^2,2));
   end
   dprime = (dprime./mean(dprime(1,:))).*mean(data(exp_num).dprime(1,:));
   model_resp(exp_num).dprime = dprime;
end

% compute cost
cost = arrayfun(@(x,y) (x.dprime(:)-y.dprime(:)).^2,model_resp,data,'uniformoutput',0);
cost = nansum(cell2mat(cost(:))); % sum across experiments
