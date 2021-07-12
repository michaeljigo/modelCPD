% Purpose:  Display fit results of fit_thresholds_exp1.m
%
% By:       Michael Jigo
%           05.25.21

function display_exp1_fits(varargin)

%% Set default parameters
in = {'sf_profile' ...     % 'narrow' or 'broad' or 'space_only'
   'model_variant' ...     % 'main_model', 'minus_sf', 'minus_ori', 'minus_space', 'minus_context', 'minus_sum'
   'attn_type'};           % 'exo' or 'endo'

val = {'narrow' ...        % sf_profile
   'main_model' ...        % model_variant
   'exo'};                 % attn_type

p = parseOptionalInputs(in,val,varargin); 

%% Load path
loaddir = '../data/fitted_parameters/';
filename = sprintf('exp1_%s.mat',p.attn_type);
load([loaddir,filename]);

% load behavior
load(sprintf('../data/behavior/exp1.mat'));
if strcmp(p.attn_type,'exo')
   data = data(1);
else
   data = data(2);
end


%% Compute group-average and SEM
cs = 1./data.crf.thresh;
avg = squeeze(mean(cs,1));
ci = get_bootstrap_ci(cs);

cue_effect = squeeze(cs(:,:,:,2)./cs(:,:,:,1));
avg_effect = squeeze(mean(cue_effect,1));
raw_effect = 1./data.crf.raw_thresh; raw_effect = raw_effect(:,:,:,2)./raw_effect(:,:,:,1);
sem_effect = withinSubjError(raw_effect);
ci_effect = get_bootstrap_ci(raw_effect);

% put model CS (neut + attention) together
modelcs(:,:,1) = out.modelcs.neut;
modelcs(:,:,2) = out.modelcs.neut.*out.modelcs.attn;


%% Get Peak SFs
   % Neutral peak SF from data
   neutSF.avg = exp(mean(log(out.data.csf.peakSF),1));
   neutSF.ci = get_bootstrap_ci(out.data.csf.peakSF);
   neutSF.sem = withinSubjError(out.data.csf.peakSF);

   % attention peak SF from model
   attnSF = out.attn.attn_freq_max+(out.attn.attn_freq_slope.*out.data.ecc);

%% Display fits
switch p.attn_type
   case 'exo'
      colors = [0 0 0; 5 113 176]./255;
   case 'endo'
      colors = [0 0 0; 202 0 32]./255;
end


% Plot contrast sensitivity functions (separate plots for Neutral and Exo)
   figure('name','CSFs (split)','position',[138 322 444 201]);
   interp_freq = logspace(log10(0.5),log10(8),12); % 11
   ecc_colors(:,:,1) = [37 37 37; 99 99 99; 150 150 150; 204 204 204]./255;

   switch p.attn_type
      case 'exo'
         ecc_colors(:,:,2) = [8 81 156; 49 130 189; 107 174 214; 189 215 231]./255;
      case 'endo'
         ecc_colors(:,:,2) = [165 15 21; 222 45 38; 251 106 74; 252 174 145]./255;
   end
   for cue = 1:size(out.data.crf.thresh,4)
      subplot(1,size(out.data.crf.thresh,4),cue);
      for e = 1:numel(out.data.ecc)

         % draw fits
         interp_cs(e,:,cue) = exp(interp1(log10(out.data.freq),log(modelcs(e,:,cue)),log10(interp_freq),'spline'));
         if cue==1
            % draw vertical lines for Neutral peak SF
            line([neutSF.avg(e) neutSF.avg(e)],[1 max(squeeze(interp_cs(e,:,cue)))],'color',ecc_colors(e,:,cue),'linewidth',2); hold on
         end
         leg(cue) = loglog(interp_freq,squeeze(interp_cs(e,:,cue)),'-','linewidth',3,'color',ecc_colors(e,:,cue)); hold on
   
         % draw subject avg 
         loglog(out.data.freq,squeeze(avg(e,1:5,cue)),'o','color','w','markersize',4.5);
         loglog(out.data.freq,squeeze(avg(e,1:5,cue)),'o','markerfacecolor',ecc_colors(e,:,cue),'markersize',4,'markeredgecolor','none');

         % draw confidence interval
         for f = 1:numel(out.data.freq)
            lb = ci(1,e,f,cue);
            ub = ci(2,e,f,cue);
            line([data.freq(f) data.freq(f)],[lb ub],'color',ecc_colors(e,:,cue),'linewidth',1.5); 
         end
      end
      % pretty up subplot
      figureDefaults
      xlabel('Spatial frequency (cpd)','fontname','arial','fontsize',10);
      ylabel('Contrast sensitivity','fontname','arial','fontsize',10);
      set(gca,'xlim',[0.3 16],'xtick',2.^(-1:4),'ylim',[1 300],'ytick',[1 10 100 300],'yscale','log','xscale','log');
   end
   
   % save figure
   savedir = '../figures/exp1/pdfs/';
   if ~exist(savedir,'dir')
      mkdir(savedir);
   end
   filename = sprintf('%s_CSF_SPLIT.pdf',p.attn_type);
   saveas(gcf,[savedir,filename]);


   
% Plot cueing effect (Valid/Neutral)
   figure('name','Cue effect','position',[138 383 561 140]);
   for e = 1:numel(out.data.ecc)
      subplot(1,numel(out.data.ecc),e);

      % draw line at 1
      line([0.125 25],[1 1],'color',[0.5 0.5 0.5],'linewidth',2); hold on
      
      
      % draw vertical lines for Neutral peak SF
      line([neutSF.avg(e) neutSF.avg(e)],[1 1.4],'color','k','linewidth',2); hold on
      lb = neutSF.avg(e)-neutSF.sem(e); ub = neutSF.avg(e)+neutSF.sem(e); % SEM
      %lb = neutSF.ci(1,e); ub = neutSF.ci(2,e);  % CI
      interp_sf = logspace(log10(lb),log10(ub),1e2);
      interp_val_min = ones(1,numel(interp_sf)); interp_val_max = ones(1,numel(interp_sf))+0.4;
      f = fill([interp_sf fliplr(interp_sf)],[interp_val_min fliplr(interp_val_max)],[0.5 0.5 0.5]); hold on
      set(f,'FaceColor',colors(1,:),'FaceAlpha',0.25,'EdgeColor','none');

     
      % draw fits from Jigo & Carrasco 2020
      semilogx(data.groupavgFit.freq,data.groupavgFit.attn(e,:),'-','linewidth',4,'color',[colors(2,:) 0.2]); hold on
      
      
      % draw fits by imAmodel
         % 1-octave spacing
         %interp_fit = interp1(log10(out.data.freq),out.modelcs.attn(e,:),log10(interp_freq),'spline');
         %semilogx(interp_freq,interp_fit,'-','linewidth',3,'color',colors(2,:)); hold on
         
         % 0.5-octave spacing
         interp_fit = interp1((out.finer_SF_sample.freq),out.finer_SF_sample.model.attn(e,:),(interp_freq),'pchip');
         semilogx(interp_freq,interp_fit,'-','linewidth',3,'color',colors(2,:)); hold on

      
      % draw subject avg 
      semilogx(data.freq,avg_effect(e,:),'o','markersize',6,'color','w'); hold on
      semilogx(data.freq,avg_effect(e,:),'o','markersize',5,'markerfacecolor',colors(2,:),'markeredgecolor','none'); hold on
      errorbar(data.freq,avg_effect(e,:),sem_effect(e,:),'.','markersize',1,'color',colors(2,:),'capsize',0,'linestyle','none','linewidth',1.5);

      % pretty up subplot
      figureDefaults
      if e==1
         xlabel('Spatial frequency (cpd)','fontname','arial','fontsize',10);
         ylabel('Attention effect','fontname','arial','fontsize',10);
      end
      set(gca,'xlim',[0.3 16],'xtick',2.^(-1:4),'ytick',1:0.2:5,'xscale','log','yscale','linear','ylim',[0.98 1.3],'ytick',[1:0.15:2]);
      title(sprintf('%i deg',out.data.ecc(e)));
   end
   
   % save figure
   savedir = '../figures/exp1/pdfs/';
   if ~exist(savedir,'dir')
      mkdir(savedir);
   end
   filename = sprintf('%s_attn_effect.pdf',p.attn_type);
   %saveas(gcf,[savedir,filename]);
