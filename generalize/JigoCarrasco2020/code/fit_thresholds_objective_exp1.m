% Purpose:  Objective function for fit_thresholds_exp1.m

function [cost model] = fit_thresholds_objective_exp1(stimdrive,supdrive,attn,data,grating,energy,p,params)

%% Re-structure parameters
   % re-scale parameters
   params = params'.*(p.unscaled_model_bnd(:,2)-p.unscaled_model_bnd(:,1))+p.unscaled_model_bnd(:,1);
   
   % re-structure
   [stimdrive,supdrive,attn] = prune_unprune_parameters(p,params,stimdrive,supdrive,attn);

   
%% Generate model discriminability to gratings
  % loop through spatial frequencies
  parfor f = 1:numel(data.freq)
      dprime(:,f,:,:) = parallel_discrim(p,data,grating,energy,stimdrive,supdrive,attn,f);
  end 


%% Determine contrast threshold at each eccentricity
   interp_contrast = logspace(-3,0,1e3);
   for f = 1:numel(data.freq)
      for e = 1:numel(data.ecc)
         for cue = 1:size(data.crf.thresh,4)
            current_dprime = squeeze(dprime(cue,f,:,e));
            interp_dprime = interp1(log10(p.contrasts),current_dprime,log10(interp_contrast),'spline');
         
            % get thresholds 
            if cue==1 && isequal(data.ecc(e),0)
               % isolate foveal CRFs for this SF
               fovealCRF = dprime(cue,f,:,e);

               % this will serve as the scalar for this SF only
               neut_scalar = 1./max(fovealCRF(:));
            end
            %neut_scalar = 1./max(current_dprime(:)); % get scalar using only the neutral condition at the current eccentricity x SF
            d = interp_dprime; d = d.*neut_scalar*2; % scale up so maximum response is 2
            all_d(cue,f,:,e) = d;
            [~,thresh_idx] = min(abs(d-1));
            modelcs(e,f,cue) = 1./interp_contrast(thresh_idx);
         end
      end
   end
   modelcs(modelcs<1) = 1;



%% Cost

% separate out neutral CS
model_neut = modelcs(:,:,1);
data_neut = 1./data.crf.thresh(:,:,:,1);

% separate out attention effects
model_attn = modelcs(:,:,2)./modelcs(:,:,1);
data_attn = 1./data.crf.thresh;
data_attn = data_attn(:,:,:,2)./data_attn(:,:,:,1);

% save model output
model.neut = model_neut;
model.attn = model_attn;

% repeat model vals to match # of observers
model_neut = repmat(shiftdim(model_neut,-1),[size(data.crf.thresh,1) 1 1]);
model_attn = repmat(shiftdim(model_attn,-1),[size(data.crf.thresh,1) 1 1]);

% log-transform neutral CS
model_neut = log(model_neut);
data_neut = log(data_neut);

% cost
cost_neut = sum((model_neut(:)-data_neut(:)).^2);
cost_attn = sum((model_attn(:)-data_attn(:)).^2);
total_cost = cost_neut+cost_attn;

% scale individual costs to have equal weight to the summed cost
neut_ratio = 1./(cost_neut./total_cost);
attn_ratio = 1./(cost_attn/total_cost);
cost = cost_neut*neut_ratio+cost_attn*attn_ratio;
   

function dprime = parallel_discrim(p,data,grating,energy,stimdrive,supdrive,attn,f)
   
   % loop through contrasts
     for c = 1:numel(p.contrasts)

        % loop through cues (neutral, valid)
        for cue = 1:size(data.crf.thresh,4)
           use_attn = cue-1;

           % evaluate model
           targresp = imAmodel(grating.targ(:,:,f,c),'energy',energy.targ(f,c).energy,'ecc',data.ecc,'stimdrive',stimdrive,'supdrive',supdrive,'attn',attn,'use_attn',use_attn,'sf_profile',p.sf_profile,'model_variant',p.model_variant,'px_per_deg',p.px_per_deg);
           notargresp = imAmodel(grating.notarg(:,:,f,c),'energy',energy.notarg(f,c).energy,'ecc',data.ecc,'stimdrive',stimdrive,'supdrive',supdrive,'attn',attn,'use_attn',use_attn,'sf_profile',p.sf_profile,'model_variant',p.model_variant,'px_per_deg',p.px_per_deg);

           % compute dprime (Eucledian norm)
           dprime(cue,c,:) = sqrt(nansum((targresp(:,:)-notargresp(:,:)).^2,2));
        end
     end
