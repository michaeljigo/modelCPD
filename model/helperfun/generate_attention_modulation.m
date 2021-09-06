% Purpose:  Compute attentional modulation across SF, space and eccentricity.
%           Given the inputted attention parameters, this function will adjust filter weights accordingly.
%
% Inputs:   channel        steerable pyramid subband information
%           attn           attention parameters
%           params         image parameters
%           ecc_im         eccentricity indices
%           eccmap         eccentricity map for individual dimensions (i.e., x and y separately)
%
% Outputs:  attn           attention structure updated with appropriate gain modulation
%
% By:       Michael Jigo

function attn = generate_attention_modulation(attn,channel,params,ecc_im,eccmap)

%% Modulate SF sensitivity with cosine function
% amplitude of peak attention effect across SF
amp = attn.attn_amp_max+(attn.attn_amp_slope*ecc_im);

% frequency where attention effects are largest
peak_freq = attn.attn_freq_max+(attn.attn_freq_slope*ecc_im);



%% Set attention gain profile
% create attention effect for each frequency channel at each eccentricity (pixel)
attn.modulation = zeros([size(peak_freq,1) 1 channel.n_freq size(peak_freq,4) size(peak_freq,5)]);
for f = 1:channel.n_freq
   freq_center = log2(channel.freq(f));
   switch params.sf_profile 
      case 'broad'
         % the sum of 3 cosine functions
         [rc rca] = make_cosine_fun(freq_center,squeeze(peak_freq(:,1,f,:,:)),(attn.attn_bw*2),1);
         attn.modulation(:,1,f,:,:) = rc+rca;
      case 'narrow'
         % a single cosine function
         attn.modulation(:,1,f,:,:) = make_cosine_fun(freq_center,squeeze(peak_freq(:,1,f,:,:)),attn.attn_bw*2,1);
      case 'space_only'
         % gain is uniform across SF
         attn.modulation(:,1,f,:,:) = ones(size(squeeze(peak_freq(:,1,f,:,:))));
   end
end


%% Adjust amplitude and spatial spread of attentional modulation
% set amplitude of attentional modulation
amplitude = attn.attn_amp_max+(ecc_im(:,1,:,:,:).*attn.attn_amp_slope);

% determine location of attention field
if isnan(params.attended_ecc)
   attended_ecc = params.ecc;
else
   if numel(params.attended_ecc)==1
      attended_ecc = repmat(params.attended_ecc,1,numel(params.ecc));
   end
end


%% Spatial spread
for e = 1:numel(params.ecc)

   % create map for spatial locus of attention field
   dist_from_center = attended_ecc(e)-params.ecc(e);
   [~,locus] = cart2pol(eccmap.x-dist_from_center,eccmap.y);
   
   % create spatial window
   switch params.spatial_profile
      case 'center'
         % excitatory center modeled as a cosine function
         tmp_spread = make_cosine_fun(locus,0,attn.attn_spread*2,1);
      case 'center_surround2DG' % center-surround profile, defined by 2nd-derivative of Gaussian 
         %tmp_spread = gausderivative2(locus,1,attn.attn_spread,0,0); % y-position always along horizontal meridian
         tmp_spread = gausderivative2(locus,0,attn.attn_spread,1,0); % y-position always along horizontal meridian
      case 'center_surround' 
         % excitatory center with suppressive surround, modeled as the sum of two Gaussians
         pos_amp     = 1;                                       % positive amplitude
         neg_amp     = 1*attn.attn_sup_amp;                     % negative amplitude
         pos_spread  = attn.attn_spread;                        % positive spread
         neg_spread  = attn.attn_spread*attn.attn_sup_spread;   % negative spread
         tmp_spread  = make_difference_of_gaussians(locus,0,pos_amp,neg_amp,pos_spread,neg_spread);
   end
   attn_spatial_spread(e,1,1,:,:) = tmp_spread;
end
% repeat dimensions to match energy matrix
attn_spatial_spread = repmat(attn_spatial_spread,[1 1 channel.n_freq 1 1]);



%% Window attentional modulation
attn.modulation = attn.modulation.*attn_spatial_spread;


%% Adjust amplitude
attn.modulation = (amplitude-attn.attn_baseline).*attn.modulation+attn.attn_baseline;


%% Repeat matrix to account for orientations
attn.modulation = repmat(attn.modulation,[1 channel.n_ori 1 1 1]);
attn.modulation(attn.modulation<0) = 0;
