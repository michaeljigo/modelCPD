% Purpose:  Initialize a set of parameters for a collection of experiments. Then prune parameters, based on index table.
% By:       Michael Jigo
%
% Important variables:
%     exp_list       = list of experiments being fit
%     param_list     = list of parameters being manipulated
%     full_index     = full matrix of fixed and free parameters for each experiment
%     fitted_params  = parameters currently being optimized
%     fixed_index    = fixed parameters

function [stimdrive supdrive attn param_vals bnds plaus_bnds fixed_index] = prune_unprune_params(p,fitted_params,stimdrive,supdrive,attn)

if nargin==1
   %% Initialize parameters
   for e = 1:numel(p.exp_list)
      [stimdrive(e) supdrive(e) attn(e)] = init_parameters('sf_profile',p.sf_profile); 

      % set repeated parameters to NaN
      for param = p.param_list
         param = char(param);
         if ~isequal(e,p.full_index(e,ismember(p.param_list,param))) % checks if parameter should follow a different experiment
            switch param
               case {'cg_max' 'cg_slope' 'freq_max' 'freq_slope' 'bw_max' 'bw_slope'}
                  stimdrive(e).(param) = nan;
               case {'attn_freq_max' 'attn_freq_slope' 'attn_bw' 'attn_amp_max' 'attn_spread'}
                  attn(e).(param) = nan;
            end
         end

         % put all parameters in a matrix matching the index matrix
         switch param
            case {'cg_max' 'cg_slope' 'freq_max' 'freq_slope' 'bw_max' 'bw_slope'}
               param_vals(e,ismember(p.param_list,param)) = stimdrive(e).(param);
               param_lo_bnd(e,ismember(p.param_list,param)) = min(stimdrive(e).bnd.(param));
               param_hi_bnd(e,ismember(p.param_list,param)) = max(stimdrive(e).bnd.(param));

               param_plaus_lo(e,ismember(p.param_list,param)) = min(stimdrive(e).plaus_bnd.(param));
               param_plaus_hi(e,ismember(p.param_list,param)) = max(stimdrive(e).plaus_bnd.(param));
            case {'attn_freq_max' 'attn_freq_slope' 'attn_bw' 'attn_amp_max' 'attn_spread'}
               param_vals(e,ismember(p.param_list,param)) = attn(e).(param);
               param_lo_bnd(e,ismember(p.param_list,param)) = min(attn(e).bnd.(param));
               param_hi_bnd(e,ismember(p.param_list,param)) = max(attn(e).bnd.(param));
               
               param_plaus_lo(e,ismember(p.param_list,param)) = min(attn(e).plaus_bnd.(param));
               param_plaus_hi(e,ismember(p.param_list,param)) = max(attn(e).plaus_bnd.(param));
         end
      end
   end

   %% if fitting neutral only, remove attention parameters
   if strcmp(p.attn_type,'neutral')
      attn = [];
   end

   %% collapse all parameters into a vector, then keep track of nan locations
   % collapse
   param_vals = param_vals(:);
   param_lo_bnd = param_lo_bnd(:);
   param_hi_bnd = param_hi_bnd(:);
   param_plaus_lo = param_plaus_lo(:);
   param_plaus_hi = param_plaus_hi(:);
   
   % remove fixed (i.e., nans)
   fixed_index = isnan(param_vals);
   param_vals(fixed_index) = [];
   param_lo_bnd(fixed_index) = [];
   param_hi_bnd(fixed_index) = [];
   param_plaus_lo(fixed_index) = [];
   param_plaus_hi(fixed_index) = [];
   bnds = [param_lo_bnd param_hi_bnd];
   plaus_bnds = [param_plaus_lo param_plaus_hi];

else
   % un-prune paramater values and recreate full matrix of parameters
   param_vals = nan(size(p.full_index));
   param_vals = param_vals(:);
   param_vals(~p.fixed_index) = fitted_params;
   param_vals = reshape(param_vals,size(p.full_index));

   for e = 1:numel(p.exp_list)
      for ii = 1:numel(p.param_list)
         switch p.param_list{ii}
            case {'cg_max' 'cg_slope' 'freq_max' 'freq_slope' 'bw_max' 'bw_slope'}
               if ~isnan(param_vals(e,ii))
                  stimdrive(e).(p.param_list{ii}) = param_vals(e,ii);
               else
                  stimdrive(e).(p.param_list{ii}) = param_vals(p.full_index(e,ii),ii);
               end 
            case {'attn_freq_max' 'attn_freq_slope' 'attn_bw' 'attn_amp_max' 'attn_spread'}
               if ~isnan(param_vals(e,ii))
                  attn(e).(p.param_list{ii}) = param_vals(e,ii);
               else
                  attn(e).(p.param_list{ii}) = param_vals(p.full_index(e,ii),ii);
               end 
         end
      end
   end
end

