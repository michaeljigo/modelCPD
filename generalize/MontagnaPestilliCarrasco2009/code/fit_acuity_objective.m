% Purpose:  Objective function for fit_acuity.m

function [cost model] = fit_acuity_objective(stimdrive,supdrive,attn,data,landolt,energy,p,params)

%% Re-structure parameters
   % re-scale parameters
   params = params'.*(p.unscaled_model_bnd(:,2)-p.unscaled_model_bnd(:,1))+p.unscaled_model_bnd(:,1);
   
   % re-structure
   [stimdrive,supdrive,attn] = prune_unprune_parameters(p,params,stimdrive,supdrive,attn);
         
   
         %stimdrive.cg_max = 1.5;
         %stimdrive.cg_slope = -0.0912;
         %stimdrive.freq_max = 3;
         %stimdrive.freq_slope = -0.9;
         %stimdrive.bw_max = 1.1;

         %attn.attn_freq_max = 2;
         %attn.attn_freq_slope = -0.3;
         %attn.attn_bw = 3.5;
         %attn.attn_amp_max = 15;
         %attn.attn_spread = 0.6;


%% Generate model discriminability to landolt squares
  % loop through gap size
   for g = 1:numel(p.gap_size)
      dprime(:,g) = parallel_discrim(p,data,landolt,energy,stimdrive,supdrive,attn,g);
   end


%% Determine gap threshold
   interp_gap = linspace(0,30,1e3);
   for cue = 1:size(data.thresh,2)
      current_dprime = squeeze(dprime(cue,:));
      interp_dprime = interp1(p.gap_size,current_dprime,interp_gap,'pchip');
         
      %if cue==1
         % scale all performance relative to the Neutral condition 
         neut_scalar = 1./max(current_dprime(:)); 
      %end
      d = interp_dprime; d = d.*neut_scalar*2; % scale up so maximum response is 2
      all_d(cue,:) = d;
      [~,thresh_idx] = min(abs(d-1));
      model_thresh(cue) = interp_gap(thresh_idx);
   end


%% Cost
acuity_diff = diff(data.thresh,[],2)./data.thresh(:,1);
model_diff = diff(model_thresh)./model_thresh(1);

% store model output
model.thresh = model_thresh;
model.dprime = all_d;

% repeat model CS to match # of observers
model_thresh  = repmat(model.thresh,[size(data.thresh,1) 1]);
data_thresh   = data.thresh;
model_diff = repmat(model_diff,size(data.thresh,1),1);

% sse
cost_thresh = sum((model_thresh(:)-data_thresh(:)).^2);
cost_diff = sum((model_diff(:)-acuity_diff(:)).^2);
total_cost = cost_thresh+cost_diff;
%cost = cost_thresh;
%cost = cost_diff;

% scale individual costs to have equal weight to the summed cost
thresh_ratio = 1./(cost_thresh./total_cost);
diff_ratio = 1./(cost_diff/total_cost);
cost = cost_thresh*thresh_ratio+cost_diff*diff_ratio;
cost = cost_thresh;

     
   
function dprime = parallel_discrim(p,data,landolt,energy,stimdrive,supdrive,attn,g)
   
      % loop through cues (neutral, valid)
   for cue = 1:size(data.thresh,2)
      use_attn = cue-1;

      % evaluate model
      targresp = imAmodel(landolt.targ(:,:,g),'energy',energy.targ(g).energy,'ecc',data.ecc,'stimdrive',stimdrive,'supdrive',supdrive,'attn',attn,'use_attn',use_attn,'sf_profile',p.sf_profile,'model_variant',p.model_variant,'px_per_deg',p.px_per_deg);
        
      notargresp = imAmodel(landolt.notarg(:,:,g),'energy',energy.notarg(g).energy,'ecc',data.ecc,'stimdrive',stimdrive,'supdrive',supdrive,'attn',attn,'use_attn',use_attn,'sf_profile',p.sf_profile,'model_variant',p.model_variant,'px_per_deg',p.px_per_deg);

      % use highest SF channel only
      targresp = targresp(:,:,1:1,:,:);
      notargresp = notargresp(:,:,1:1,:,:);
   
      % compute dprime (Eucledian norm)
      dprime(cue) = sqrt(nansum((targresp(:,:)-notargresp(:,:)).^2,2));
     end
