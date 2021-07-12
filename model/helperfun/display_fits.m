% Purpose:  Display fit results of fit_exps.m
%
% By:       Michael Jigo
%           05.17.21

function display_fits(varargin)

%% Set default parameters
in = {'sf_profile' ...   % 'narrow' or 'broad' or 'space_only'
   'model_variant' ...     % 'main_model', 'minus_sf', 'minus_ori', 'minus_space', 'minus_context', 'minus_sum'
   'attn_type' ...         % 'neutral' or 'involuntary' or 'voluntary'
   'display_fit'};         % 1=display fits; 0=don't display fits

val = {'narrow' ...        % sf_profile
   'main_model' ...        % model_variant
   'voluntary' ...       % attn_type
   1};                     % display_fit

p = parseOptionalInputs(in,val,varargin); 

%% Load path
loaddir = sprintf('../data/fitted_parameters/%s/',p.model_variant);
filename = sprintf('%s_%s.mat',p.attn_type,p.sf_profile);
load([loaddir,filename]);



%% Display fits
figure('position',[207 1402 1375 246]);
for e = 1:numel(out.p.exp_list)
   subplot(1,numel(out.p.exp_list),e);
   for cue = 1:size(out.data(e).dprime,1)
      if cue==1 % Neutral
         % observed data
         err = errorbar(out.data(e).ecc,out.data(e).dprime(cue,:),out.data(e).sem(cue,:),'linestyle','none','capsize',0,'linewidth',2,'color','k'); hold on
         set([err.Bar, err.Line],'ColorType','truecoloralpha','colordata',[err.Line.ColorData(1:3); 255*0.25]);
         % group-average
         plot(out.data(e).ecc,out.data(e).dprime(cue,:),'o','markersize',6,'color','w','markerfacecolor','w'); hold on % white edge
         plot(out.data(e).ecc,out.data(e).dprime(cue,:),'o','markersize',5,'markerfacecolor','k','markeredgecolor','none') % colored dot;

         % model fit
         leg(cue) = plot(out.data(e).ecc,out.model_resp(e).dprime(cue,:),'-','linewidth',4,'color','k');  
      else % valid
         % set appropriate colors
         switch out.p.exp_list{e}
            case {'yc98_exp1' 'yc98_exp2' 'tc02' 'clh06' 'yc08' 'ymc08_exp2'} % involuntary attention experiments
               col = [5 113 176]./255;
            otherwise % voluntary attention experiments
               col = [202 0 32]./255;
         end

         err = errorbar(out.data(e).ecc,out.data(e).dprime(cue,:),out.data(e).sem(cue,:),'linestyle','none','capsize',0,'linewidth',2,'color',col); hold on
         set([err.Bar, err.Line],'ColorType','truecoloralpha','colordata',[err.Line.ColorData(1:3); 255*0.25]);
         % group-average
         plot(out.data(e).ecc,out.data(e).dprime(cue,:),'o','markersize',6,'color','w','markerfacecolor','w'); hold on % white edge
         plot(out.data(e).ecc,out.data(e).dprime(cue,:),'o','markersize',5,'markerfacecolor',col,'markeredgecolor','none') % colored dot;

         % model fit
         leg(cue) = plot(out.data(e).ecc,out.model_resp(e).dprime(cue,:),'-','linewidth',4,'color',col);  
      end
   end

   %% Compute AIC value
   nfree = size(out.p.model_bnd,1);
   nobs = 0;
   for ex = 1:numel(out.p.exp_list)
      nobs = nobs+numel(out.data(ex).dprime); 
   end
   aic = nobs*log(out.sse./nobs)+(2*nfree);
   aicc = aic+((2*nfree*(nfree+1))./(nobs-nfree-1));

  
   %% pretty up figure
   % x limits
   xlim(1) = -0.5;
   if max(out.data(e).ecc)<8
      xlim(2) = 8;
   elseif max(out.data(e).ecc)<=10
      xlim(2) = 10;
   elseif max(out.data(e).ecc)<12
      xlim(2) = 12;
   elseif max(out.data(e).ecc)<=20
      xlim(2) = 20;
   else
      xlim(2) = 24;
   end

   % y limits
   minval = min(out.data(e).dprime(:))-max(out.data(e).sem(:)); maxval = max(out.data(e).dprime(:))+max(out.data(e).sem(:));
   bin = 0.5;
   ymin = floor(minval) + floor( (minval-floor(minval))/bin) * bin; ymin = max(ymin,-0.25);
   ymax = floor(maxval) + ceil( (maxval-floor(maxval))/bin) * bin;
   ylim = [0 ymax];
   ytick = -5:0.5:5;
   yticklabel = cellfun(@num2str,num2cell(ytick),'UniformOutput',0);
   set(gca,'box','off','tickdir','out','linewidth',2,'fontname','arial','fontsize',8,'xtick',0:4:40,'ytick',ytick,'ylim',ylim,'xlim',xlim,'PlotBoxAspectRatio',[1 1 1],'yticklabel',yticklabel);
   xlabel('Eccentricity (deg)','fontname','arial','fontsize',10);
   ylabel('d-prime','fontname','arial','fontsize',10);

   if e==1
      title(sprintf('%s; AICc=%.2f',out.p.exp_list{e},aicc),'interpreter','none');
      %legend(leg,{'neutral' 'valid'},'location','best');
   else
      title(out.p.exp_list{e},'interpreter','none');
   end
end

savedir = '../figures/fits/';
if ~exist(savedir,'dir')
   mkdir(savedir);
end
filename = sprintf('%s_%s_%s.png',out.p.model_variant,out.p.attn_type,out.p.sf_profile);
saveas(gcf,[savedir,filename]);
