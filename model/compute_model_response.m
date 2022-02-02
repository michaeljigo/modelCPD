% Purpose:  Compute model response, given input parameters.
%           Several model alternatives can be evaluated:
%           'main_model'      best-fitting model, use this
%           'minus_ori'       without cross-orientation suppression
%           'minus_sf'        without cross-frequency suppression
%           'minus_space'     without surround suppression
%           'minus_context'   without all contextual modulation
%           'minus_sum'       without spatial summation
%
% Input:    energy         4D matrix of steerable pyramid-derived energy
%           stimdrive      stimulus drive parameters
%           supdrive       suppressive drive parameters
%           attn           attention modulation parameters
%           model_variant  model variant string
%
% Output:   popresp        population response
%
% 
function popresp  = compute_model_response(energy,stimdrive,supdrive,attn,model_variant)

%% Attention modulation
if ~exist('attn','var') || isempty(attn)
   attn = 1;
end

%% Choose model variant to evaluate
switch model_variant
   case 'main_model' 
      popresp = main_model(energy,stimdrive,supdrive,attn);
   case 'minus_context' 
      % no contextual modulation
      popresp = minus_context(energy,stimdrive,supdrive,attn);
   case 'minus_sum'
      % no spatial summation
      popresp = minus_sum(energy,stimdrive,supdrive,attn);
   case 'minus_sf'
      % no cross-frequency suppression
      popresp = minus_sf(energy,stimdrive,supdrive,attn);
   case 'minus_ori'
      % no cross-orientation suppression
      popresp = minus_ori(energy,stimdrive,supdrive,attn);
   case 'minus_space'
      % no surround suppression
      popresp = minus_space(energy,stimdrive,supdrive,attn);
   case 'response_gain'
      % attention scales responses after normalization and before pooling
      popresp = response_gain(energy,stimdrive,supdrive,attn);
   case 'contrast_gain'
      % attention scales the suppressive drive after pooling
      popresp = contrast_gain(energy,stimdrive,supdrive,attn); otherwise
      error('%s is not a model option.',model_variant);
end
return
   
   function popresp = main_model(energy,stimdrive,supdrive,attn)
      
      %% Excitation
      % adjust gain on spatial frequency channels
      energy = stimdrive.sf_gain.*energy;
   
      %% Attention modulation
      exc_energy = attn.*energy;

      %% Lateral suppression
      inh_energy = exc_energy;
      inh_energy = pool_space(inh_energy,supdrive,1); % surround suppression
      inh_energy = pool_feature(inh_energy,supdrive,'ori'); % cross-orientiation
      inh_energy = pool_feature(inh_energy,supdrive,'freq'); % cross-frequency

      %% Normalization
      energy_contrast_gain = supdrive.contrast_gain.^2; % adjust overall sensitivity
      popresp = exc_energy./(energy_contrast_gain+inh_energy);

      %% Spatial pooling
      popresp = pool_space(popresp,supdrive,0);
   return
   
   function popresp = response_gain(energy,stimdrive,supdrive,attn)
      
      %% Excitation
      % adjust gain on spatial frequency channels
      exc_energy = stimdrive.sf_gain.*energy;
   
      %% Lateral suppression
      inh_energy = exc_energy;
      inh_energy = pool_space(inh_energy,supdrive,1); % surround suppression
      inh_energy = pool_feature(inh_energy,supdrive,'ori'); % cross-orientiation
      inh_energy = pool_feature(inh_energy,supdrive,'freq'); % cross-frequency

      %% Normalization
      energy_contrast_gain = supdrive.contrast_gain.^2; % adjust overall sensitivity
      popresp = exc_energy./(energy_contrast_gain+inh_energy);
      
      %% Attention modulation on the population response
      popresp = attn.*popresp;

      %% Spatial pooling
      popresp = pool_space(popresp,supdrive,0);
   return
   
   function popresp = contrast_gain(energy,stimdrive,supdrive,attn)
      
      %% Excitation
      % adjust gain on spatial frequency channels
      exc_energy = stimdrive.sf_gain.*energy;
   
      %% Lateral suppression
      inh_energy = exc_energy;
      inh_energy = pool_space(inh_energy,supdrive,1); % surround suppression
      inh_energy = pool_feature(inh_energy,supdrive,'ori'); % cross-orientiation
      inh_energy = pool_feature(inh_energy,supdrive,'freq'); % cross-frequency

      %% Attention modulation on the lateral suppression component
      %inh_energy = inh_energy./attn;

      %% Normalization
      energy_contrast_gain = supdrive.contrast_gain.^2; % adjust overall sensitivity
      popresp = exc_energy./(energy_contrast_gain./max(attn(:))+inh_energy);
      
      %% Spatial pooling
      popresp = pool_space(popresp,supdrive,0);
   return
   
   
   
   function popresp = minus_sf(energy,stimdrive,supdrive,attn)
   
      %% Excitation
      % adjust gain on spatial frequency channels
      energy = stimdrive.sf_gain.*energy;

      %% Attention modulation
      exc_energy = attn.*energy;

      %% Lateral suppression
      inh_energy = exc_energy;
      inh_energy = pool_space(inh_energy,supdrive,1); % surround suppression
      inh_energy = pool_feature(inh_energy,supdrive,'ori'); % cross-orientiation
   
      %% Normalization
      energy_contrast_gain = supdrive.contrast_gain.^2; % adjust overall sensitivity
      popresp = exc_energy./(energy_contrast_gain+inh_energy);

      %% Spatial pooling
      popresp = pool_space(popresp,supdrive,0);
   return
   
      
   function popresp = minus_ori(energy,stimdrive,supdrive,attn)
   
      %% Excitation
      energy = stimdrive.sf_gain.*energy;

      %% Attention modulation
      exc_energy = attn.*energy;

      %% Lateral suppression
      inh_energy = exc_energy;
      inh_energy = pool_space(inh_energy,supdrive,1); % surround suppression
      inh_energy = pool_feature(inh_energy,supdrive,'freq'); % cross-frequency

      %% Normalization
      energy_contrast_gain = supdrive.contrast_gain.^2; % adjust overall sensitivity
      popresp = exc_energy./(energy_contrast_gain+inh_energy);

      %% Spatial pooling
      popresp = pool_space(popresp,supdrive,0);
   return
   
      
   function popresp = minus_space(energy,stimdrive,supdrive,attn)
      
      %% Excitation
      energy = stimdrive.sf_gain.*energy;

      %% Attention modulation
      exc_energy = attn.*energy;

      %% Lateral suppression
      inh_energy = exc_energy;
      inh_energy = pool_feature(inh_energy,supdrive,'ori'); % cross-orientation
      inh_energy = pool_feature(inh_energy,supdrive,'freq'); % cross-frequency

      %% Normalization
      energy_contrast_gain = supdrive.contrast_gain.^2; % adjust overall sensitivity
      popresp = exc_energy./(energy_contrast_gain+inh_energy);

      %% Spatial pooling
      popresp = pool_space(popresp,supdrive,0);
   return
   
   
   function popresp = minus_context(energy,stimdrive,supdrive,attn)
      
      %% Excitation
      energy = stimdrive.sf_gain.*energy;

      %% Attention modulation
      % apply attention after exponentiation
      exc_energy = attn.*energy;

      %% Local suppression
      inh_energy = exc_energy; % each pixel normalizes itself

      %% Normalization
      energy_contrast_gain = supdrive.contrast_gain.^2; % adjust overall sensitivity
      popresp = exc_energy./(energy_contrast_gain+inh_energy);

      %% Spatial pooling
      popresp = pool_space(popresp,supdrive,0);
   return
   
   
   function popresp = minus_sum(energy,stimdrive,supdrive,attn)
      
      %% Excitation
      energy = stimdrive.sf_gain.*energy;

      %% Attention modulation
      exc_energy = attn.*energy;

      %% Lateral suppression
      inh_energy = exc_energy;
      inh_energy = pool_space(inh_energy,supdrive,1); % surround suppression
      inh_energy = pool_feature(inh_energy,supdrive,'ori'); % cross-orientiation
      inh_energy = pool_feature(inh_energy,supdrive,'freq'); % cross-frequency

      %% Normalization
      energy_contrast_gain = supdrive.contrast_gain.^2; % adjust overall sensitivity
      popresp = exc_energy./(energy_contrast_gain+inh_energy);
   return
   

   
   
   
   
   
%% Pooling functions
function pooled_energy = pool_space(energy,inh,norm_kernel)
   % Usage:    pool = pool_space(energy,inh)
   % Purpose:  Convolve each frequency channel with its corresopnding spatial kernel 
   % Input:    energy   --    energy image that will be convolved
   %           inh      --    structure containing convolution kernels
   %
   % Output:   pooled_energy     --    pooled response across space

   % pre-allocate pool matrix
   pooled_energy = zeros(size(energy));
   switch numel(size(inh.space_kernels))
      case 2 % input is a 1D vector, that will be separably convolved
         for f = 1:size(energy,2)
            % choose the appropriate kernel for this SF
            kernel = inh.space_kernels(f,:);
            kernel = kernel'*kernel;
            if norm_kernel
               kernel = kernel./sum(kernel(:));
            end
            for o = 1:size(energy,1)
               % convolve each orientaiton channel with the kernel
               ori_channel = squeeze(energy(o,f,:,:));
               pooled_energy(o,f,:,:) = convolve2(ori_channel,kernel,'replicate');
            end
         end
      case 3 % input is an image
      for f = 1:size(energy,2)
            kernel = squeeze(inh.space_kernels(f,:,:));
            if norm_kernel
               kernel = kernel./sum(kernel(:));
            end
            for o = 1:size(energy,1)
               % convolve each orientaiton channel with the kernel
               ori_channel = squeeze(energy(o,f,:,:));
               pooled_energy(o,f,:,:) = conv_fft2(ori_channel,kernel,'wrap');
            end
         end
   end
   


function pooled_energy = pool_feature(energy,inh,which_feature)
   % Usage:    pool = pool_feature(energy,inh,which_feature)
   % Purpose:  Weight and sum (i.e., convolve) across bands that are included in the suppressive pool.
   % Input:    energy         --    energy to-be-convolved
   %           inh            --    structure containing convolution kernels
   %           which_feature  --    string defining which feature will be convolved ('freq' or 'ori')
   %
   % Output:   pooled_energy  --    pooled response across feature


   % extract corresponding pool parameters
   switch which_feature
      case 'freq'
         dim = 2;
         kernel = inh.freq_pool;
      case 'ori'
         dim = 1;
         kernel = inh.ori_pool;
   end
   
   % permute matrix such that the desired feature (i.e., dim) is the first dimension
   org_idx = 1:length(size(energy));
   new_idx = [dim setdiff(org_idx,dim)];
   % keep index of original dimension order
   [~,reorder] = ismember(new_idx,org_idx);
   % permute matrix
   energy = permute(energy,new_idx);
   
   % pre-allocate pool matrix
   pooled_energy = zeros(size(energy));
   
   % perform pooling 
   for f = 1:size(energy,1)
      % choose the appropriate kernel for this channel
      k = find(kernel(f,:)>0);
      
      % create weight matrix that will be applied to kernel
      rep_dim = size(energy);
      feature_pool = nansum(bsxfun(@times,energy(k,:,:,:),kernel(f,k)'),1); % faster
      
      % construct pool matrix
      pooled_energy(f,:,:,:) = feature_pool;
   end
   % revert matrix to original ordering
   pooled_energy = permute(pooled_energy,reorder);
