% Purpose:  This function will load:
%           1) a parameter set for a given experiment and 
%           2) the experiment's image then
%           Run the image through the model with the specified parameters to recreate a fit displayed in the manuscript.
% 
% Input:    exp_name       desired experiment to be evaluated
%           display_fit    1=display; 0=don't display
%
% Output:   model_dprime   matrix of model dprime [neutral;valid]
%
% By:       Michael Jigo

function model_dprime = generate_best_fit(exp_name,display_fit)
addpath(genpath('./helperFun'));
if ~exist('display_fit','var')
   display_fit = 1;
end

%% Load experiment data and parameters
filename = sprintf('../data/behavior/%s.mat',exp_name);
load(filename);

% add baseline parameter

%params.stimdrive.cg_max = 2.9199;
%params.stimdrive.cg_slope = -0.0579;
%params.stimdrive.freq_max = 1.9364;
%params.stimdrive.freq_slope = -0.5126;
%params.stimdrive.bw_max = 1.523;
%params.attn.attn_freq_max = 2.2955;
%params.attn.attn_freq_slope = -0.1516;
%params.attn.attn_bw = 2.8239;
%params.attn.attn_amp_max = 5.9131;
%params.attn.attn_amp_slope = 0;
params.attn.attn_spread = 4;
params.attn.attn_baseline = 1;
params.attn.attn_sup_amp = 0;
params.attn.attn_sup_spread = 0;



%% Choose appropriate gain profile
switch exp_name
   case {'yc98_exp1' 'yc98_exp2' 'tc02' 'clh06' 'yc08' 'ymc08_exp2'} % involuntary attention experiments
      gain_profile = 'narrow';
   otherwise % voluntary attention experiments
      gain_profile = 'broad';
end


%% Run image through model with specified parameters
% Neutral
neut_targ = imAmodel(data.stim.targ,'ecc',data.ecc,'stimdrive',params.stimdrive,'supdrive',params.supdrive,'attn',params.attn,'use_attn',0);
neut_notarg = imAmodel(data.stim.notarg,'ecc',data.ecc,'stimdrive',params.stimdrive,'supdrive',params.supdrive,'attn',params.attn,'use_attn',0);

% Valid
[valid_targ p] = imAmodel(data.stim.targ,'ecc',data.ecc,'stimdrive',params.stimdrive,'supdrive',params.supdrive,'attn',params.attn,'use_attn',1,'sf_profile',gain_profile);
[valid_notarg p] = imAmodel(data.stim.notarg,'ecc',data.ecc,'stimdrive',params.stimdrive,'supdrive',params.supdrive,'attn',params.attn,'use_attn',1,'sf_profile',gain_profile);

%% Compute d-prime
neut_dprime = sqrt(nansum((neut_targ(:,:)-neut_notarg(:,:)).^2,2));
valid_dprime = sqrt(nansum((valid_targ(:,:)-valid_notarg(:,:)).^2,2));

% scale to match behavioral data
valid_dprime = (valid_dprime./mean(neut_dprime))*mean(data.dprime(1,:));
neut_dprime = (neut_dprime./mean(neut_dprime))*mean(data.dprime(1,:));

% put model dprime into a matrix
model_dprime(1,:) = neut_dprime;
model_dprime(2,:) = valid_dprime;


%% Plot
if display_fit
   figure('name',exp_name);
   for cue = 1:size(data.dprime,1)
      if cue==1 % Neutral
         % observed data
         err = errorbar(data.ecc,data.dprime(cue,:),data.sem(cue,:),'linestyle','none','capsize',0,'linewidth',2,'color','k'); hold on
         set([err.Bar, err.Line],'ColorType','truecoloralpha','colordata',[err.Line.ColorData(1:3); 255*0.25]);
         % group-average
         plot(data.ecc,data.dprime(cue,:),'o','markersize',6,'color','w','markerfacecolor','w'); hold on % white edge
         plot(data.ecc,data.dprime(cue,:),'o','markersize',5,'markerfacecolor','k','markeredgecolor','none') % colored dot;

         % model fit
         leg(cue) = plot(data.ecc,model_dprime(cue,:),'-','linewidth',4,'color','k');  
      else % valid
         % set appropriate colors
         switch exp_name
            case {'yc98_exp1' 'yc98_exp2' 'tc02' 'clh06' 'yc08' 'ymc08_exp2'} % involuntary attention experiments
               col = [5 113 176]./255;
            otherwise % voluntary attention experiments
               col = [202 0 32]./255;
         end

         err = errorbar(data.ecc,data.dprime(cue,:),data.sem(cue,:),'linestyle','none','capsize',0,'linewidth',2,'color',col); hold on
         set([err.Bar, err.Line],'ColorType','truecoloralpha','colordata',[err.Line.ColorData(1:3); 255*0.25]);
         % group-average
         plot(data.ecc,data.dprime(cue,:),'o','markersize',6,'color','w','markerfacecolor','w'); hold on % white edge
         plot(data.ecc,data.dprime(cue,:),'o','markersize',5,'markerfacecolor',col,'markeredgecolor','none') % colored dot;
         
         % model fit
         leg(cue) = plot(data.ecc,model_dprime(cue,:),'-','linewidth',4,'color',col);  
      end
   end

   %% pretty up figure
   % x limits
   xlim(1) = -0.5;
   if max(data.ecc)<8
      xlim(2) = 8;
   elseif max(data.ecc)<=10
      xlim(2) = 10;
   elseif max(data.ecc)<12
      xlim(2) = 12;
   elseif max(data.ecc)<=20
      xlim(2) = 20;
   else
      xlim(2) = 24;
   end

   % y limits
   minval = min(data.dprime(:))-max(data.sem(:)); maxval = max(data.dprime(:))+max(data.sem(:));
   bin = 0.5;
   ymin = floor(minval) + floor( (minval-floor(minval))/bin) * bin; ymin = max(ymin,-0.25);
   ymax = floor(maxval) + ceil( (maxval-floor(maxval))/bin) * bin;
   ylim = [0 ymax];
   ytick = -5:0.5:5;
   yticklabel = cellfun(@num2str,num2cell(ytick),'UniformOutput',0);
   set(gca,'box','off','tickdir','out','linewidth',2,'fontname','arial','fontsize',8,'xtick',0:4:40,'ytick',ytick,'ylim',ylim,'xlim',xlim,'PlotBoxAspectRatio',[1 1 1],'yticklabel',yticklabel);
   xlabel('Eccentricity (deg)','fontname','arial','fontsize',10);
   ylabel('d-prime','fontname','arial','fontsize',10);
   title(exp_name,'interpreter','none');
   legend(leg,{'neutral' 'valid'},'location','northeast');
end
