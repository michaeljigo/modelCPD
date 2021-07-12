% Purpose:  This function will use the best-fitting parameters and evaluate the model at a finer range of SFs.
%
% By:       Michael Jigo
%           06.02.21

function evaluate_exp1_newSF(varargin)

%% Set default parameters
in = {'sf_profile' ...     % 'narrow' or 'broad' or 'space_only'
   'model_variant' ...     % 'main_model', 'minus_sf', 'minus_ori', 'minus_space', 'minus_context', 'minus_sum'
   'attn_type' ...         % 'exo' or 'endo'
   'px_per_deg' ...        % pixels per degree of images
   'im_width' ...          % side length of image in degrees
   'freq_bw' ...           % SF frequency bandwidth of model
   'display_fit'};         % 1=display fits; 0=don't display fits

val = {'narrow' ...        % sf_profile
   'main_model' ...        % model_variant
   'exo' ...               % attn_type
   32 ...                  % px_per_deg
   4 ...                   % im_width
   0.5 ...                 % freq_bw
   1};                     % display_fit

p = parseOptionalInputs(in,val,varargin); 

%% Load best-fitting parameters
   load(sprintf('../data/fitted_parameters/exp1_%s.mat',p.attn_type));


%% Set new SFs
   data = out.data;
   %data.freq = 2.^(-1:0.5:3);
   data.freq = 2.^(-1:1:3);
   data.freq = 2;

   stimdrive = out.stimdrive;
   supdrive = out.supdrive;
   attn = out.attn;
   attn.attn_amp_max = 4;


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
         [~,energy.targ(f,c)] = imAmodel(tmp45,'decompose_only',1,'px_per_deg',p.px_per_deg,'channel_freq_bw',p.freq_bw);
         [~,energy.notarg(f,c)] = imAmodel(tmp135,'decompose_only',1,'px_per_deg',p.px_per_deg,'channel_freq_bw',p.freq_bw);
      end
   end

   
%% Generate model discriminability to gratings
  % loop through spatial frequencies
  for f = 1:numel(data.freq)
      dprime(:,f,:,:) = parallel_discrim(p,data,grating,energy,stimdrive,supdrive,attn,f);
  end 


%% Determine contrast threshold at each eccentricity
   interp_contrast = logspace(-3,0,1e3);
   for f = 1:numel(data.freq)
      for e = 1:numel(data.ecc)
         for cue = 1:size(data.crf.thresh,4)
            current_dprime = squeeze(dprime(cue,f,:,e));
            interp_dprime = interp1(log(p.contrasts),current_dprime,log(interp_contrast),'spline');
         
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


%% Get model output
% separate out neutral CS
model_neut = modelcs(:,:,1);

% separate out attention effects
model_attn = modelcs(:,:,2)./modelcs(:,:,1);

% save model output
model.neut = model_neut;
model.attn = model_attn;
keyboard


function dprime = parallel_discrim(p,data,grating,energy,stimdrive,supdrive,attn,f)
   
   % loop through contrasts
     for c = 1:numel(p.contrasts)

        % loop through cues (neutral, valid)
        for cue = 1:size(data.crf.thresh,4)
           use_attn = cue-1;

           % evaluate model
           targresp = imAmodel(grating.targ(:,:,f,c),'energy',energy.targ(f,c).energy,'ecc',data.ecc,'stimdrive',stimdrive,'supdrive',supdrive,'attn',attn,'use_attn',use_attn,'sf_profile',p.sf_profile,'model_variant',p.model_variant,'px_per_deg',p.px_per_deg,'channel_freq_bw',p.freq_bw);
           notargresp = imAmodel(grating.notarg(:,:,f,c),'energy',energy.notarg(f,c).energy,'ecc',data.ecc,'stimdrive',stimdrive,'supdrive',supdrive,'attn',attn,'use_attn',use_attn,'sf_profile',p.sf_profile,'model_variant',p.model_variant,'px_per_deg',p.px_per_deg,'channel_freq_bw',p.freq_bw);

      
           % compute dprime (Eucledian norm)
           dprime(cue,c,:) = sqrt(nansum((targresp(:,:)-notargresp(:,:)).^2,2));
        end
     end

