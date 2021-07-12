% Purpose:  Recreate Figure 4 showing Neutral performance and statistical analyses.
%
% By:       Michael Jigo
% Edited:   07.06.21

function figure4
addpath(genpath('../../modelCPD_v4'));

%% Preprocess and format data before displaying
   % visualization parameters
      disp.ci     = 0.16;     % 68% confidence interval 
      disp.extrap = 'pchip';  % extrapolation method
      disp.colors = [0 0 0]; 
      

   % load bootstrap samples for all models fit to neutral data
      models = {'main_model' 'minus_ori' 'minus_space' 'minus_sf' 'minus_context' 'minus_sum'};
      for m = 1:numel(models)
         boot.(models{m}) = combine_bootstrap_iterations('model_variant',models{m},'attn_type','neutral');
      end


   % get indices of experiments to display
      explist  = boot.(models{1}).sim_params.exp_list;
      exps     = {'yc98_exp1' 'yc98_exp2'};
      expidx   = find(ismember(explist,exps));


   % interpolate median and 68% confidence interval of model responses
      for e = 1:numel(expidx)
         expdata              = boot.main_model.exp_data(expidx(e));
         modelresp            = boot.main_model.model_resp(expidx(e)).dprime;
         medianresp           = squeeze(nanmedian(modelresp,3));
         modelci              = squeeze(quantile(modelresp,[disp.ci 1-disp.ci],3));
         modeldisp(e).ecc     = 0:(1/3):max(expdata.ecc);
         modeldisp(e).center  = interp1(expdata.ecc,medianresp,modeldisp(e).ecc,disp.extrap);
         modeldisp(e).ci      = interp1(expdata.ecc,modelci,modeldisp(e).ecc,disp.extrap)';

      % add in error for each eccentricity and validity condition
         err = load(sprintf('../data/behavior/%s.mat',exps{e}));
         sem                          = err.data.sem;

         % compute upper and lower bound based on SEM
         y = boot.main_model.exp_data(expidx(e)).dprime;
         modeldisp(e).err.lb          = y-sem(1,:);
         modeldisp(e).err.ub          = y+sem(1,:);
      end


   % concatenate data for scatter plot
      scatterdata.obs   = [];
      scatterdata.model = [];
      for e = 1:numel(boot.main_model.exp_data)
         % get index for best-fitting parameters for this experiment (for visualization purposes)
            [~,bestidx] = min(boot.main_model.fit_err(:,e));
            bestmodel   = boot.main_model.model_resp(e).dprime(:,:,bestidx);

         % store the data
            scatterdata.obs   = cat(1,scatterdata.obs,boot.main_model.exp_data(e).dprime(:));
            scatterdata.model = cat(1,scatterdata.model,bestmodel(:));
      end
      % best-fit regression line
         [scatterdata.fit, scatterdata.ci] = regress(scatterdata.model,[ones(size(scatterdata.obs)) scatterdata.obs]);


   % compute delta AIC and BIC scores for model comparison figure
      for m = 1:numel(models)
         % compute comparison metrics
         aicval         = boot.(models{m}).aic.boot-boot.main_model.aic.boot;
         bicval         = boot.(models{m}).bic.boot-boot.main_model.bic.boot;

         % store median
         aic.center(m)  = median(aicval);
         bic.center(m)  = median(bicval);

         % store 95% confidence interval
         aic.ci(:,m)    = quantile(aicval,[0.025 1-0.025]);
         bic.ci(:,m)    = quantile(bicval,[0.025 1-0.025]);
      end


   % generate matrix images showing SF selectivity
      for e = 1:numel(expidx)
         % load texture stimulus for simulation
            stim        = load(sprintf('../data/behavior/%s.mat',exps{e}));         
            stim        = stim.data.stim;

         % get index for best-fitting parameters for this experiment
            [~,bestidx] = min(boot.main_model.fit_err(:,expidx(e)));

         % get parameters for simulation
            ecc         = 0:1:20;
            params      = {'cg_max' 'cg_slope' 'freq_max' 'freq_slope' 'bw_max'};
            [stimdrive, supdrive] = init_parameters;
            fitparams   = boot.main_model.free_params(expidx(e));
            for p = 1:numel(params)
               stimdrive.(params{p}) = fitparams.(params{p})(bestidx);
            end

         % run the simulation
            [targresp tmp] = imAmodel(stim.targ,'ecc',ecc,'stimdrive',stimdrive,'supdrive',supdrive,'use_attn',0);
            notargresp     = imAmodel(stim.notarg,'ecc',ecc,'stimdrive',stimdrive,'supdrive',supdrive,'use_attn',0);

         % compute absolute differences between target and no-target
            dresp       = abs(targresp-notargresp);
            % average across orientation and space
            dresp       = squeeze(mean(mean(mean(dresp,2),4),5));

         % store matrix for display
            freqmat.vals(:,:,e)  = flip(dresp,2)';
            freqmat.ecc          = ecc;
            freqmat.freq         = tmp.energy.channel.freq;
      end


%% Display
figure('name','Figure 4','position',[680 950 572 748]);
   % CPD fits
      subplotidx  = [1 3];
      xlims       = [0 12; 0 24];
      xtick       = 0:4:24; 
      ylim        = [0 2];
      ytick       = 0:0.5:2;
      titles      = {'fine-scale' 'coarse-scale'};

      for e = 1:numel(modeldisp)
         subplot(3,2,subplotidx(e));
         % draw model fits (error)
            x = [modeldisp(e).ecc, fliplr(modeldisp(e).ecc)];
            y = [modeldisp(e).ci(1,:) fliplr(modeldisp(e).ci(2,:))];
            f = fill(x,y,disp.colors);
            set(f,'facecolor',disp.colors,'facealpha',0.25,'edgecolor','none'); hold on

         % draw model fits (median)
            plot(modeldisp(e).ecc,modeldisp(e).center,'-','linewidth',3,'color',disp.colors);

         % draw observed data
            x = boot.main_model.exp_data(expidx(e)).ecc;
            y = boot.main_model.exp_data(expidx(e)).dprime;
            plot(x,y,'o','markersize',5,'color','w','markerfacecolor','w');
            plot(x,y,'o','markersize',4,'markerfacecolor',disp.colors,'markeredgecolor','none');

         % draw error bars
            for ec = 1:size(y,2)
               lb = modeldisp(e).err.lb(ec);
               ub = modeldisp(e).err.ub(ec);
               line([x(ec) x(ec)],[lb ub],'color',[disp.colors 0.5],'linewidth',2); 
            end

         % pretty up subplot
            figureDefaults
            set(gca,'xlim',xlims(e,:),'xtick',xtick,'ylim',ylim,'ytick',ytick,'ticklength',[0.025 0.05]);
            xlabel('Eccentricity (^o)','fontname','arial','fontsize',10); 
            ylabel('Performance (d^{\prime})','fontname','arial','fontsize',10); 
            title(titles{e},'fontname','arial','fontsize',10);
      end

   
   % Matrix of absolute differences
      subplotidx  = [2 4];
      yticklabel  = cellfun(@num2str,num2cell(fliplr(freqmat.freq)),'uniformoutput',0);
      for e = 1:size(freqmat.vals,3)
         subplot(3,2,subplotidx(e));
         imagesc(freqmat.ecc,1:size(freqmat.vals,1),freqmat.vals(:,:,e)); colormap gray; figureDefaults;

         % pretty up figure
            set(gca,'yscale','linear','ytick',1:5,'yticklabel',yticklabel,'xtick',0:4:24,'ydir','normal','ticklength',[0.025 0.05]);
            xlabel('Eccentricity (^o)','fontname','arial','fontsize',10);
            ylabel('Spatial frequency (cpd)','fontname','arial','fontsize',10);
      end


   % Scatter plot
      subplot(3,2,5);
      xlim = [0 2];
      % unity line
         line([0 2],[0 2],'color',disp.colors,'linestyle','--','linewidth',1.5); hold on
      
      % regression line
         % confidence interval
         x  = linspace(0,2,1e3);
         lb = scatterdata.ci(2,1)*x+scatterdata.ci(1,1);
         ub = scatterdata.ci(2,2)*x+scatterdata.ci(1,2);
         f  = fill([x fliplr(x)],[lb fliplr(ub)],disp.colors); hold on
         set(f,'facecolor',disp.colors,'facealpha',0.15,'edgecolor','none');
         % regression line
         y  = scatterdata.fit(2)*x+scatterdata.fit(1);
         plot(x,y,'-','linewidth',3,'color',disp.colors);

      % data points
         plot(scatterdata.obs,scatterdata.model,'o','markersize',5,'color','w','markerfacecolor','w');
         plot(scatterdata.obs,scatterdata.model,'o','markersize',4,'markerfacecolor',disp.colors,'markeredgecolor','none');

      % pretty it up
         figureDefaults
         set(gca,'xlim',xlim,'ylim',ylim,'ticklength',[0.025 0.05]);
         xlabel('Measured d^{\prime}','fontname','arial','fontsize',10);
         ylabel('Predicted d^{\prime}','fontname','arial','fontsize',10);


   % Model comparisons
      subplot(3,2,6)
      aic_xval  = [0 0.5 0.75 1 1.25 1.75];
      bic_xval  = [0 0.5 0.75 1 1.25 1.75]+0.1;
      xval(:,1) = aic_xval; 
      xval(:,2) = bic_xval; 
      xval      = mean(xval,2);

      % aic
         for ii = 1:numel(aic_xval)
            line([aic_xval(ii) aic_xval(ii)],[aic.ci(1,ii) aic.ci(2,ii)],'color',[0.5 0.5 0.5],'linewidth',1.5); hold on
         end
         leg(1) = plot(aic_xval,aic.center,'o','markersize',4,'markerfacecolor',[0.5 0.5 0.5],'markeredgecolor','k');

      % bic
         for ii = 1:numel(bic_xval)
            line([bic_xval(ii) bic_xval(ii)],[bic.ci(1,ii) bic.ci(2,ii)],'color','k','linewidth',1.5);
         end
         leg(2) = plot(bic_xval,bic.center,'o','markersize',4,'markerfacecolor','w','markeredgecolor','k');

      % pretty it up
         figureDefaults
         legend(leg,{'\Delta AIC' '\Delta BIC'},'location','northwest');
         set(gca,'ticklength',[0.025 0.05],'xtick',xval,'xticklabel',{'full' '-\theta' '-x,y' '-f' '-all' '-sum'},'ylim',[-2 50],'ytick',0:10:50,'xlim',[-0.25 2]);
         ylabel('Model performance','fontname','arial','fontsize',10);

saveas(gcf,'./fig4.png');
