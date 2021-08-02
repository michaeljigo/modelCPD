% Purpose:  Evaluate log-parabola function for the inputted SFs (or filter center frequencies) at the given eccentricity.
%           
% Input:    Structure containing the eccentricities and frequencies. As well as the parameters controlling the 
%           height and slope of each of the three parameters (amplitude, peak_freq, bandwidth) of the log-parabola.
%
% Output:   Evaluated values of the log-parabola function

function [weights amp peak_freq bw] = evaluate_log_parabola(varargin)

%% Set log-parabola parameters
% compute_sensitivity=1, evaluate the CSF
% compute_sensitivity=0, just output the parameters
in = {'amp_max' 'amp_slope' 'freq_max' 'freq_slope' 'freq_min' 'bw_max' 'bw_slope' 'ecc' 'freq' 'display_function' 'channel_freq' 'channel_ori' 'compute_sensitivity'};
val = {2 -0.0483 2 -0.128 0.5 2 0 0:2:10 2 0 [] [] 1};
params = parseOptionalInputs(in,val,varargin);


%% Evaluate parameter values for the desired eccentricity
% amplitude
amp = 10.^(params.amp_max+(params.amp_slope*params.ecc));

% peak SF 
peak_freq = 2.^(log2(2^params.freq_max-params.freq_min)+(params.freq_slope*params.ecc))+params.freq_min;

% bandwidth
bw = (params.bw_max)*exp(params.bw_slope.*params.ecc);

if ~params.compute_sensitivity
   weights = [];
   return
end

if numel(size(amp))==2
   % input is not an image

   % ensure that values are in a column vector
   if size(amp,2)>1
      amp = amp';
   end
   if size(peak_freq,2)>1
      peak_freq = peak_freq';
   end
   if size(bw,2)>1
      bw = bw';
   end
   % concatenate parameters into single matrix
   logparab_params = [peak_freq bw amp];
end


%% Generate the appropriate contrast sensitivity (weights) for the SFs at each eccentricity
if numel(size(amp))==2
   % input is a vector of eccentricities
   for e = 1:numel(params.ecc)
      weights(e,:) = evalCSF('logparab',params.freq,logparab_params(e,:));
   end
   %weights(weights<1) = 1;
elseif numel(size(amp)==3)
   % input is an eccentricity image
  
   % duplicate channel frequencies to all eccentricities and pixels
   freq_mat = shiftdim(params.channel_freq,-3);
   rep_sizes = size(amp); rep_sizes(3) = [];
   freq_mat = repmat(freq_mat,[rep_sizes 1]);
   freq_mat = permute(freq_mat,[1 2 5 3 4]);
  
   % matrix of log parabola parameter values for each pixel (eccentricity) 
   logparab_params(:,:,:,:,:,1) = peak_freq;
   logparab_params(:,:,:,:,:,2) = bw;
   logparab_params(:,:,:,:,:,3) = amp;

   % generate weights at each pixel
   weights = evalCSF('logparab',freq_mat,logparab_params);
end


%% Display the CSFs at each eccentricity
if params.display_function
   fit_freq = logspace(log10(0.5),log10(128),1e3);
   %colors = linspecer(numel(params.ecc));
   colors = linspecer(11);
   
   figure('Name','Evaluated log-parabola'); 
   for e = 1:numel(params.ecc)
      disp_weights(e,:) = evalCSF('logparab',fit_freq,logparab_params(e,:)); disp_weights(disp_weights<1) = nan;
      loglog(fit_freq,disp_weights(e,:),'-','linewidth',4,'color',colors(e,:)); hold on
   end
   set(gca,'box','off','tickdir','out','xtick',2.^(-3:6),'ytick',[1 10 100 500 1000],'xlim',[0.4 32],'fontname','arial','fontsize',14);
   legend(cellfun(@num2str,num2cell(params.ecc),'uniformoutput',0),'location','southwest');
   xlabel('Spatial frequency (cpd)','fontname','arial','fontsize',16); ylabel('Sensitivity (1/threshold)','fontname','arial','fontsize',16);
end
