% Purpose:  Objective function for fit_thresholds_exp2.m

function [cost model] = fit_thresholds_objective_exp2(stimdrive,supdrive,attn,data,grating,energy,p,params)

%% Re-structure parameters
   % re-scale parameters
   params = params'.*(p.unscaled_model_bnd(:,2)-p.unscaled_model_bnd(:,1))+p.unscaled_model_bnd(:,1);
   
   % re-structure
   [stimdrive,supdrive,attn] = prune_unprune_parameters(p,params,stimdrive,supdrive,attn);

   
   %stimdrive.cg_max = 1.9762;
   %stimdrive.cg_slope = -0.0601;
   %stimdrive.freq_max = 3;
   %stimdrive.freq_slope = -0.0912;
   %stimdrive.bw_max = 6;
   %stimdrive.bw_slope = 0;

   %attn.attn_freq_max = 2.7624;
   %attn.attn_freq_slope = -0.0105;
   %attn.attn_bw = 2.3008;
   %attn.attn_amp_max = 1.6454;
   %attn.attn_spread = 1.1207;
   %attn.baseline = 1;


%% Generate model discriminability to gratings
  % loop through spatial frequencies
  parfor f = 1:numel(data.freq) 
      dprime(:,f,:) = parallel_discrim(p,data,grating,energy,stimdrive,supdrive,attn,f);
  end 


%% Scale model dprime to match data
% compute scalar
   avgdata = squeeze(mean(data.performance,1));
   neutavg = mean(avgdata(1,:));
   scalar = neutavg/mean(dprime(1,:));

% scale
   dprime = dprime*scalar;


%% Separate Neutral performance and Attention effect
   % data
   data_neut = squeeze(data.performance(:,1,:,:));
   data_attn = data.attn_effect;

   % model
   model_neut = squeeze(dprime(1,:,:));
   model_attn = squeeze(diff(dprime,[],1));

   % repeat model to match # of subjects
   nsubj = size(data.performance,1);
   model_neut = repmat(shiftdim(model_neut,-1),[nsubj 1 1]);
   model_attn = repmat(shiftdim(model_attn,-1),[nsubj 1 1]);
   

%% Cost
% sse (summed across Neutral performance and attention effects)
cost_neut = sum((model_neut(:)-data_neut(:)).^2);
cost_attn = sum((model_attn(:)-data_attn(:)).^2);
total_cost = cost_neut+cost_attn;

% scale individual costs to have equal weight to the summed cost
neut_ratio = 1./(cost_neut./total_cost);
attn_ratio = 1./(cost_attn/total_cost);
cost = cost_neut*neut_ratio+cost_attn*attn_ratio;

% put model values into a structure
model.neut = squeeze(model_neut(1,:,:));
model.attn = squeeze(model_attn(1,:,:));
     
   
function dprime = parallel_discrim(p,data,grating,energy,stimdrive,supdrive,attn,f)
   
   % loop through contrasts
     for e = 1:numel(data.ecc)

      % the no-target image doesn't depend on eccentricity
        % loop through cues (neutral, valid)
        for cue = 1:2
           use_attn = cue-1;

           % evaluate model
           targresp = imAmodel(grating.targ(:,:,f,e),'energy',energy.targ(f,e).energy,'ecc',p.im_center,'stimdrive',stimdrive,'supdrive',supdrive,'attn',attn,'use_attn',use_attn,'sf_profile',p.sf_profile,'model_variant',p.model_variant,'px_per_deg',p.px_per_deg,'attended_ecc',data.ecc(e),'channel_freq_bw',0.5);
      
           notargresp = imAmodel(grating.notarg(:,:,f),'energy',energy.notarg(f).energy,'ecc',p.im_center,'stimdrive',stimdrive,'supdrive',supdrive,'attn',attn,'use_attn',use_attn,'sf_profile',p.sf_profile,'model_variant',p.model_variant,'px_per_deg',p.px_per_deg,'attended_ecc',data.ecc(e),'channel_freq_bw',0.5);
      
           % compute dprime (Eucledian norm)
           dprime(cue,e) = sqrt(nansum((targresp(:,:)-notargresp(:,:)).^2,2));
        end
     end
