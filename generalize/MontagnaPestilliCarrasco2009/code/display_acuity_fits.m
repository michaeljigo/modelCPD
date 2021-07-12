% Purpose:  Display fit results of fit_acuity.m
%
% By:       Michael Jigo
%           05.29.21

function display_acuity_fits(varargin)

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
filename = sprintf('%s.mat',p.attn_type);
load([loaddir,filename]);


%% Compute group-average and SEM
avg_thresh = mean(out.data.thresh);
sem_thresh = std(out.data.thresh,[],1)./sqrt(size(out.data.thresh,1));
ci_thresh = get_bootstrap_ci(out.data.thresh); % 68% confidence interval

thresh_diff = (out.data.thresh(:,1)-out.data.thresh(:,2))./out.data.thresh(:,1)*100;
avg_thresh_diff = mean(thresh_diff);
sem_thresh_diff = std(thresh_diff)./sqrt(size(out.data.thresh,1));


%% Display fits
switch p.attn_type
   case 'exo'
      colors = [0 0 0; 5 113 176]./255;
   case 'endo'
      colors = [0 0 0; 202 0 32]./255;
end

% Plot thresholds from Montagna, Pestilli & Carrasco 2009
   figure('name','DATA thresholds','position',[138 315 207 208]);
   
   % draw the avg 
   xval = [1 1.4];
   for ii = 1:numel(xval)
      bar(xval(ii),avg_thresh(ii),0.3,'edgecolor','none','facecolor',[0.5 0.5 0.5]); hold on
   end

   % draw SEM
   %errorbar(xval,avg_thresh,sem_thresh,'linestyle','none','capsize',0,'color','k','linewidth',2);

   % draw confidence intervals
   for ii = 1:numel(xval)
      lb = ci_thresh(1,ii); ub = ci_thresh(2,ii);
      line([xval(ii) xval(ii)],[lb ub],'color',[0 0 0],'linewidth',3); hold on
   end
   
   % add in model predictions
   plot(xval,out.model.thresh,'o','markersize',6,'markerfacecolor','w','markeredgecolor','none'); hold on
   leg(2) = plot(xval,out.model.thresh,'o','markersize',5,'markerfacecolor',colors(2,:),'markeredgecolor','none'); hold on

   % pretty up plot
   figureDefaults
   ylabel('Acuity threshold (arc min)','fontname','arial','fontsize',10);
   set(gca,'xlim',[min(xval)-0.3 max(xval)+0.3],'xtick',xval,'xticklabel',{'Neutral' 'Valid'},'ylim',[0 15],'ytick',0:5:50);
   %legend(leg,{'data' 'model'});

   % save figure
   savedir = '../figures/pdfs/';
   if ~exist(savedir,'dir')
      mkdir(savedir);
   end
   filename = sprintf('%s_DATA.pdf',p.attn_type);
   %saveas(gcf,[savedir,filename]);
   return


% Display psychometric functions
   figure('name','psychometric function','position',[360 357 317 261]);
   interp_gap = linspace(0,30,1e3);
   for cue = 1:numel(out.model.thresh)

      % draw vertical line at the threshold
      line([out.model.thresh(cue) out.model.thresh(cue)],[0 1],'color',colors(cue,:),'linewidth',1.5); hold on

      % draw horizontal line to threshold
      line([0 out.model.thresh(cue)],[1 1],'color',colors(cue,:),'linewidth',1.5); hold on

      % draw interpolated tuning functions
      interp_dprime = interp1(linspace(0,30,1e3),out.model.dprime(cue,:),interp_gap,'spline');
      leg(cue) = plot(interp_gap,interp_dprime,'-','linewidth',3,'color',colors(cue,:)); hold on
   end
   legend(leg,{'Neutral' 'Valid'},'location','southeast');
   figureDefaults;
   xlabel('Gap size (arc min)','fontname','arial','fontsize',10);
   ylabel('Discriminability (d^\prime)','fontname','arial','fontsize',10);
   set(gca,'ytick',0:0.5:2,'xlim',[0 30],'xtick',0:10:60,'ylim',[0 2]);
   
   
   % save figure
   savedir = '../figures/pdfs/';
   if ~exist(savedir,'dir')
      mkdir(savedir);
   end
   filename = sprintf('%s_PSYCHO_FUN.pdf',p.attn_type);
   %saveas(gcf,[savedir,filename]);
   return


%% Create landolt squares
gap_size = [0 15 30];
for g = 1:numel(gap_size)
   % create target and no-target images
   landolt = make_landolt_square('gap_size',gap_size(g),'im_size',3,'px_per_deg',32,'background_color','white');

   figure; imagesc(landolt); colormap gray; truesize; set(gca,'visible','off');
   
   % save figure
   savedir = '../figures/pdfs/';
   if ~exist(savedir,'dir')
      mkdir(savedir);
   end
   filename = sprintf('LANDOLT_%i.pdf',gap_size(g));
   saveas(gcf,[savedir,filename]);
end
