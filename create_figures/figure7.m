% Purpose:  Recreate Figure 7 showing model comparison results of attention.
%
% By:       Michael Jigo

function figure7
addpath(genpath('../../modelCPD'));

   % load bootstrap samples for all models fit to exogenous and endogenous attention experiments
      attntype = {'exo' 'endo'};
      models   = {'match' 'mismatch' 'space_only'};
      for a = 1:numel(attntype)
         for m = 1:numel(models)
            boot.(attntype{a}).(models{m}) = combine_bootstrap_iterations('attn_type',attntype{a},'sf_profile',models{m});
         end
      end

   % compare models
      bestmodel   = 'match';
      for a = 1:numel(attntype)
         % arrange output as [narrow broad space] models
            if strcmp(attntype{a},'exo')
               permuteidx = [1 2 3];
            else
               permuteidx = [2 1 3];
            end
         for m = 1:numel(models)
            % get proper index
               n = permuteidx(m);

            % distribution
               aic.(attntype{a}).raw(:,n) = boot.(attntype{a}).(models{m}).aic.boot-boot.(attntype{a}).(bestmodel).aic.boot;
               bic.(attntype{a}).raw(:,n) = boot.(attntype{a}).(models{m}).bic.boot-boot.(attntype{a}).(bestmodel).bic.boot;

            % median
               aic.(attntype{a}).center(n) = median(boot.(attntype{a}).(models{m}).aic.boot-boot.(attntype{a}).(bestmodel).aic.boot);
               bic.(attntype{a}).center(n) = median(boot.(attntype{a}).(models{m}).bic.boot-boot.(attntype{a}).(bestmodel).bic.boot);

            % ci
               aic.(attntype{a}).ci(:,n) = quantile(boot.(attntype{a}).(models{m}).aic.boot-boot.(attntype{a}).(bestmodel).aic.boot,[0.025 1-0.025]);
               bic.(attntype{a}).ci(:,n) = quantile(boot.(attntype{a}).(models{m}).bic.boot-boot.(attntype{a}).(bestmodel).bic.boot,[0.025 1-0.025]);
         end
      end


   % Plot
      figure('name','Figure 7');
      for a = 1:numel(attntype)
         subplot(1,2,a)
         aic_xval  = [0 0.5 1];
         bic_xval  = [0 0.5 1]+0.05;
         xval(:,1) = aic_xval; 
         xval(:,2) = bic_xval; 
         xval      = mean(xval,2);

         % aic
            for ii = 1:numel(aic_xval)
               line([aic_xval(ii) aic_xval(ii)],[aic.(attntype{a}).ci(1,ii) aic.(attntype{a}).ci(2,ii)],'color',[0.5 0.5 0.5],'linewidth',1.5); hold on
            end
            leg(1) = plot(aic_xval,aic.(attntype{a}).center,'o','markersize',4,'markerfacecolor',[0.5 0.5 0.5],'markeredgecolor','k');

         % bic
            for ii = 1:numel(bic_xval)
               line([bic_xval(ii) bic_xval(ii)],[bic.(attntype{a}).ci(1,ii) bic.(attntype{a}).ci(2,ii)],'color','k','linewidth',1.5);
            end
            leg(2) = plot(bic_xval,bic.(attntype{a}).center,'o','markersize',4,'markerfacecolor','w','markeredgecolor','k');

         % pretty it up
            figureDefaults
            title(attntype{a},'fontname','arial','fontsize',10);
            if strcmp(attntype{a},'exo')
               legend(leg,{'\Delta AIC' '\Delta BIC'},'location','southeast');
            else
               legend(leg,{'\Delta AIC' '\Delta BIC'},'location','northeast');
            end
            set(gca,'ticklength',[0.025 0.05],'xtick',xval,'xticklabel',{'narrow' 'broad' 'space'},'ylim',[-2 50],'ytick',0:20:100,'xlim',[-0.25 1.2]);
            ylabel('Model performance','fontname','arial','fontsize',10);
      end


   % Save figure
      saveas(gcf,'./figure7.pdf');
