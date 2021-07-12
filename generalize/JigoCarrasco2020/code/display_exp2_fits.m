% Purpose:  Display fit results of fit_thresholds_exp2.m
%
% By:       Michael Jigo
%           05.26.21

function display_exp2_fits(varargin)

%% Set default parameters
in = {'sf_profile' ...     % 'narrow' or 'broad' or 'space_only'
   'model_variant' ...     % 'main_model', 'minus_sf', 'minus_ori', 'minus_space', 'minus_context', 'minus_sum'
   'attn_type'};           % 'neutral' or 'involuntary' or 'voluntary'

val = {'narrow' ...        % sf_profile
   'main_model' ...        % model_variant
   'exo'};                 % attn_type

p = parseOptionalInputs(in,val,varargin); 

%% Load path
loaddir = '../data/fitted_parameters/';
filename = sprintf('exp2_%s.mat',p.attn_type);
load([loaddir,filename]);


%% Compute group-average and SEM
avg = squeeze(mean(out.data.performance,1));
sem = withinSubjError(out.data.performance);

cue_effect = squeeze(mean(out.data.attn_effect,1));
sem_cue    = withinSubjError(out.data.attn_effect);


%% Get Peak SFs
   % Neutral peak SF from data
   neutSF.avg = exp(mean(log(out.data.csf.peakSF),1));
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
interpfreq = logspace(log10(0.5),log10(11),30);

% Plot Neutral performance (DATA)
   figure('name','Neutral performance (DATA)','position',[138 272 851 251]);
   for e = 1:numel(out.data.ecc)
      % draw fits
      interp_neut = interp1(log10(out.data.freq),out.model.neut(:,e),log10(interpfreq),'spline');
      semilogx(out.data.freq,out.model.neut(:,e),'-','linewidth',3,'color',colors(e,:),'markersize',20); hold on
      %semilogx(interpfreq,interp_neut,'-','linewidth',3,'color',colors(e,:),'markersize',20); hold on
      
      % draw subject avg with SEM
      semilogx(out.data.freq,squeeze(avg(1,:,e)),'o','markersize',6,'color','w'); hold on
      semilogx(out.data.freq,squeeze(avg(1,:,e)),'o','markersize',5,'markerfacecolor',colors(e,:),'markeredgecolor','none'); hold on
      errorbar(out.data.freq,squeeze(avg(1,:,e)),squeeze(sem(1,:,e)),'.','color',colors(e,:),'markersize',1,'capsize',0,'linestyle','none'); hold on

      % pretty up subplot
      figureDefaults
      xlabel('Spatial frequency (cpd)','fontname','arial','fontsize',10);
      ylabel('Performance (d^{\prime})','fontname','arial','fontsize',10);
      set(gca,'xlim',[0.3 16],'xtick',2.^(-1:4),'ylim',[0 2.5],'ytick',0:0.5:10,'xscale','log');
   end
   legend(cellfun(@num2str,num2cell(out.data.ecc),'uniformoutput',0),'location','southwest');

   % save figure
   savedir = '../figures/exp2/';
   if ~exist(savedir,'dir')
      mkdir(savedir);
   end
   filename = sprintf('%s_neutral_DATA.png',p.attn_type);
   saveas(gcf,[savedir,filename]);

   
% Plot Neutral performance (MODEL)
   %figure('name','Neutral performance (MODEL)','position',[138 272 851 251]);
   %for e = 1:numel(out.data.ecc)
      % draw fits
      %semilogx(out.data.freq,out.model.neut(:,e),'.','linewidth',3,'color',colors(e,:),'markersize',20); hold on
   
      % draw subject avg with SEM
      %errorbar(out.data.freq,squeeze(avg(1,:,e)),squeeze(sem(1,:,e)),'.','color',colors(e,:),'markersize',20,'capsize',0,'linestyle','none');

      % pretty up subplot
      %figureDefaults
      %xlabel('Spatial frequency (cpd)','fontname','arial','fontsize',10);
      %ylabel('Performance (d^{\prime})','fontname','arial','fontsize',10);
      %set(gca,'xlim',[0.3 16],'xtick',2.^(-1:4),'ylim',[0 2.5],'ytick',0:0.5:10,'xscale','log');
   %end
   %legend(cellfun(@num2str,num2cell(out.data.ecc),'uniformoutput',0),'location','southwest');

   % save figure
   %savedir = '../figures/exp2/';
   %if ~exist(savedir,'dir')
      %mkdir(savedir);
   %end
   %filename = sprintf('%s_neutral_MODEL.png',p.attn_type);
   %saveas(gcf,[savedir,filename]);

   
% Plot attention effect (delta d-prime)
   figure('name','Attention effect','position',[138 272 851 251]);
   for e = 1:numel(out.data.ecc)
      subplot(1,numel(out.data.ecc),e);

      % draw line at 0
      line([0.125 25],[0 0],'color',[0.5 0.5 0.5],'linewidth',2); hold on

      % draw vertical lines for Neutral peak SF
      line([neutSF.avg(e) neutSF.avg(e)],[0 1],'color','k','linewidth',2); hold on
      lb = neutSF.avg(e)-neutSF.sem(e); ub = neutSF.avg(e)+neutSF.sem(e); 
      interp_sf = logspace(log10(lb),log10(ub),1e2);
      interp_val_min = zeros(1,numel(interp_sf)); interp_val_max = ones(1,numel(interp_sf));
      f = fill([interp_sf fliplr(interp_sf)],[interp_val_min fliplr(interp_val_max)],[0.5 0.5 0.5]); hold on
      set(f,'FaceColor',colors(1,:),'FaceAlpha',0.25,'EdgeColor','none');
      
      
      % draw vertical lines for Attention peak SF
      if strcmp(p.attn_type,'exo')
         %line([attnSF(e) attnSF(e)],[0 1],'color',colors(2,:),'linewidth',2); hold on
      end
      
      
      % draw fits
      interpfit = interp1(log10(out.data.freq),squeeze(out.model.attn(:,e)),log10(interpfreq),'spline');
      %semilogx(out.data.freq,squeeze(out.model.attn(:,e)),'-','linewidth',3,'color',colors(2,:),'markersize',20); hold on
      semilogx(interpfreq,interpfit,'-','linewidth',3,'color',colors(2,:)); hold on
   
      % draw subject avg and SEM
      semilogx(out.data.freq,cue_effect(:,e),'o','markersize',6,'color','w'); hold on
      semilogx(out.data.freq,cue_effect(:,e),'o','markersize',5,'markerfacecolor',colors(2,:),'markeredgecolor','none'); hold on
      errorbar(out.data.freq,cue_effect(:,e),sem_cue(:,e),'.','markersize',1,'capsize',0,'linewidth',1.5,'linestyle','none','color',colors(2,:));
      
      % pretty up subplot
      figureDefaults
      xlabel('Spatial frequency (cpd)','fontname','arial','fontsize',10);
      ylabel('Attention effect (\Delta d^{\prime})','fontname','arial','fontsize',10);
      set(gca,'xlim',[0.3 16],'xtick',2.^(-1:4),'ylim',[0 0.5],'ytick',0:0.25:2,'xscale','log','yscale','linear');
      title(sprintf('%i deg',out.data.ecc(e)));
   end
   
   % save figure
   savedir = '../figures/exp2/';
   if ~exist(savedir,'dir')
      mkdir(savedir);
   end
   filename = sprintf('%s_attention.png',p.attn_type);
   saveas(gcf,[savedir,filename]);
