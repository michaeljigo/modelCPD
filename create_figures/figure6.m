% Purpose:  Recreate Figure 6 showing Endogenous attention performance and statistical analyses.
%
% By:       Michael Jigo

function figure6
addpath(genpath('../../modelCPD'));

%% Preprocess and format data before displaying
   % visualization parameters
      disp.ci     = 0.16;     % 68% confidence interval 
      disp.extrap = 'pchip';  % extrapolation method
      disp.colors = [0 0 0; 202 0 32]./255; 


   % load bootstrap samples for exogenous attention fits
      boot = combine_bootstrap_iterations('model_variant','main_model','attn_type','endo');


   % get indices of experiments to display
      explist  = boot.sim_params.exp_list;
      exps     = {'ymc08_exp3' 'ymc08_exp4'};
      expidx   = find(ismember(explist,exps));

   % concatenate data for scatter plot
      scatterdata.all.obs   = [];
      scatterdata.all.model = [];
      for e = 1:numel(boot.exp_data)
         % get index for best-fitting parameters for this experiment (for visualization purposes)
            [~,bestidx] = min(boot.fit_err(:,e));
            bestmodel   = boot.model_resp(e).dprime(:,:,bestidx);

         % store the data
            for c = 1:size(boot.exp_data(e).dprime,1)
               scatterdata.obs{e}(c,:)   = boot.exp_data(e).dprime(c,:);
               scatterdata.model{e}(c,:) = bestmodel(c,:);
            end
            scatterdata.all.obs   = cat(1,scatterdata.all.obs,boot.exp_data(e).dprime(:));
            scatterdata.all.model = cat(1,scatterdata.all.model,bestmodel(:));
      end
      % best-fit regression line
         [scatterdata.fit, scatterdata.ci] = regress(scatterdata.all.model(:),...
            [ones(size(scatterdata.all.obs(:))) scatterdata.all.obs(:)]);



   % interpolate median and 68% confidence interval of model responses
      for e = 1:numel(expidx)
         expdata              = boot.exp_data(expidx(e));
         modelresp            = boot.model_resp(expidx(e)).dprime;
         medianresp           = squeeze(nanmedian(modelresp,3));
         modelci              = squeeze(quantile(modelresp,[disp.ci 1-disp.ci],3));
         modeldisp(e).ecc     = linspace(0,max(expdata.ecc),20);

         for c = 1:size(expdata.dprime,1)
            modeldisp(e).center(c,:)  = interp1(expdata.ecc,medianresp(c,:),modeldisp(e).ecc,disp.extrap);
            modeldisp(e).ci(c,:,:)    = interp1(expdata.ecc,squeeze(modelci(c,:,:)),modeldisp(e).ecc,disp.extrap)';
         end

      % add in error for each eccentricity and validity condition
         err = load(sprintf('../data/behavior/%s.mat',exps{e}));
         sem                          = err.data.sem;

         % compute upper and lower bound based on SEM
         y = boot.exp_data(expidx(e)).dprime;
         modeldisp(e).err.lb          = y-sem;
         modeldisp(e).err.ub          = y+sem;
      end



   % compute attention decline across eccentricity
      ecc = linspace(0,20,50);
      for b = 1:size(boot.fit_err,1)
         for e = 1:size(boot.fit_err,2)
            freqmax     = boot.free_params(e).attn_freq_max(b);
            freqslope   = boot.free_params(e).attn_freq_slope(b);
            attndisp.all_center(:,e,b) = freqmax+(freqslope*ecc);
         end
      end
      attndisp.ecc      = ecc;
      attndisp.ci       = quantile(2.^permute(attndisp.all_center(:,:),[2 1]),[0.16 1-0.16]);
      attndisp.center   = squeeze(median(2.^attndisp.all_center(:,:),2))';



   % compute attentional gain on the target 
      % comparison will be made to the average target frequency
      exp_targfreq = mean(boot.targfreq.all,2);
      gain_on_targ.ecc              = linspace(0,20,50); 
      for b = 1:size(boot.fit_err,1) % loop through bootstrap samples
         for e = 1:size(boot.fit_err,2) % ...and experiments
            for ec = 1:numel(gain_on_targ.ecc) % and compute the gain on the target frequency at each eccentricity
               targfreq          = exp_targfreq(e);               % target SF
               attnfreq          = attndisp.all_center(ec,e,b);   % attentional center SF
               attnbw            = boot.free_params(e).attn_bw(b);% attentional bandwidth
               
               % make SF profile
               [rc rca] = make_cosine_fun(log2(targfreq),attnfreq,attnbw*2/3,1);
               gaintarg(ec,e,b)  = rc+rca;
            end
         end
      end
      % fine textures
         gain_on_targ.fine.all      = gaintarg(:,boot.targfreq.idx.fine,:); 
         gain_on_targ.fine.center   = mean(gain_on_targ.fine.all(:,:),2)';
         gain_on_targ.fine.ci       = [gain_on_targ.fine.center-std(gain_on_targ.fine.all(:,:),[],2)'; ...
                                       gain_on_targ.fine.center+std(gain_on_targ.fine.all(:,:),[],2)']; 

      % coarse textures
         gain_on_targ.coarse.all    = gaintarg(:,boot.targfreq.idx.coarse,:); 
         gain_on_targ.coarse.center = mean(gain_on_targ.coarse.all(:,:),2)';
         gain_on_targ.coarse.ci     = [gain_on_targ.coarse.center-std(gain_on_targ.coarse.all(:,:),[],2)'; ...
                                       gain_on_targ.coarse.center+std(gain_on_targ.coarse.all(:,:),[],2)']; 


   % generate matrix images showing SF selectivity of attention
      % load texture stimulus for simulation
         stim        = load(sprintf('../data/behavior/%s.mat','ymc08_exp1'));         
         modelp      = stim.params;
         stim        = stim.data.stim;

      % get index for best-fitting parameters for this experiment
         [~,bestidx] = min(boot.fit_err(:,1));
         bestidx = 58;

      % get parameters for simulation
         ecc         = 0:1:20;
         params      = {'cg_max' 'cg_slope' 'freq_max' 'freq_slope' 'bw_max' 'attn_freq_max' ...
                        'attn_freq_slope' 'attn_bw' 'attn_amp_max'};
         [stimdrive, supdrive, attn] = init_parameters;
         fitparams   = boot.free_params(1);
         stimdrive = modelp.stimdrive;
         attn      = modelp.attn;

            [stimdrive, supdrive, attn] = init_parameters;
            fitparams   = boot.free_params(1);
            for p = 1:numel(params)
               if ismember(params{p},fieldnames(stimdrive))
                  stimdrive.(params{p}) = fitparams.(params{p})(bestidx);
               end
               if ismember(params{p},fieldnames(attn))
                  attn.(params{p})      = fitparams.(params{p})(bestidx);
               end
            end
            attn.attn_baseline = 1;
            attn.attn_spread   = 4;

      % run the simulation
         [~, tmp] = imAmodel(stim.targ,'ecc',ecc,'stimdrive',stimdrive,'supdrive',supdrive,'attn',attn,'use_attn',1,'sf_profile','broad','spatial_profile','center');
      % get the model's attention modulation across eccentricity
         modulation           = tmp.attn.modulation;
         % average across orientation and space
         modulation           = squeeze(mean(mean(mean(modulation,2),4),5));

      % store matrix for display
         modmat.vals          = flip(modulation,2)';
         modmat.ecc           = ecc;
         modmat.freq          = tmp.energy.channel.freq;
         

%% Plot
figure('name','Figure 6','position',[680 950 572 748]);
   % CPD fits
      subplotidx  = [1 2];
      xlims       = [0 12; 0 24];
      xtick       = 0:4:24; 
      ylim        = [0 2.5];
      ytick       = 0:0.5:5;
      titles      = {'fine-scale' 'coarse-scale'};

      for e = 1:numel(modeldisp)
         subplot(3,2,subplotidx(e));
         for c = 1:size(modeldisp(e).center,1)
            % draw model fits (error)
               x = [modeldisp(e).ecc, fliplr(modeldisp(e).ecc)];
               y = [squeeze(modeldisp(e).ci(c,1,:))' fliplr(squeeze(modeldisp(e).ci(c,2,:))')];
               f = fill(x,y,disp.colors(c,:));
               set(f,'facecolor',disp.colors(c,:),'facealpha',0.25,'edgecolor','none'); hold on

            % draw model fits (median)
               plot(modeldisp(e).ecc,modeldisp(e).center(c,:),'-','linewidth',3,'color',disp.colors(c,:));
   
            % draw observed data
               x = boot.exp_data(expidx(e)).ecc;
               y = boot.exp_data(expidx(e)).dprime(c,:);
               plot(x,y,'o','markersize',5,'color','w','markerfacecolor','w');
               plot(x,y,'o','markersize',4,'markerfacecolor',disp.colors(c,:),'markeredgecolor','none');

            % draw error bars
               for ec = 1:size(y,2)
                  lb = modeldisp(e).err.lb(c,ec);
                  ub = modeldisp(e).err.ub(c,ec);
                  line([x(ec) x(ec)],[lb ub],'color',[disp.colors(c,:) 0.5],'linewidth',2); 
               end
         end
         % pretty up subplot
            figureDefaults
            set(gca,'xlim',xlims(e,:),'xtick',xtick,'ylim',ylim,'ytick',ytick,'ticklength',[0.025 0.05]);
            xlabel('Eccentricity (^o)','fontname','arial','fontsize',10); 
            ylabel('Performance (d^{\prime})','fontname','arial','fontsize',10); 
            title(titles{e},'fontname','arial','fontsize',10);
      end



   % Matrix of absolute differences
      yticklabel  = cellfun(@num2str,num2cell(fliplr(modmat.freq)),'uniformoutput',0);
         subplot(3,2,3);
         imagesc(modmat.ecc,1:size(modmat.vals,1),modmat.vals); colormap gray; figureDefaults;

         % pretty up figure
            set(gca,'yscale','linear','ytick',1:5,'yticklabel',yticklabel,'xtick',0:4:24,'ydir','normal','ticklength',[0.025 0.05]);
            xlabel('Eccentricity (^o)','fontname','arial','fontsize',10);
            ylabel('Spatial frequency (cpd)','fontname','arial','fontsize',10);



   % Scatter plot
      clear leg
      subplot(3,2,4);
      xlim = [0 2.5];
      xtick = 0:0.5:5;
      % unity line
         line(xlim,xlim,'color',disp.colors(1,:),'linestyle','--','linewidth',1.5); hold on
      
      % regression line
         % confidence interval
         x  = linspace(0,max(xlim),1e3);
         lb = scatterdata.ci(2,1)*x+scatterdata.ci(1,1);
         ub = scatterdata.ci(2,2)*x+scatterdata.ci(1,2);
         f  = fill([x fliplr(x)],[lb fliplr(ub)],disp.colors(1,:)); hold on
         set(f,'facecolor',disp.colors(1,:),'facealpha',0.15,'edgecolor','none');
         % regression line
         y  = scatterdata.fit(2)*x+scatterdata.fit(1);
         plot(x,y,'-','linewidth',3,'color',disp.colors(1,:));

      % data points
         for e = 1:numel(scatterdata.obs)
            for c = 1:size(scatterdata.obs{e},1)
               plot(scatterdata.obs{e}(c,:),scatterdata.model{e}(c,:),'o','markersize',5,'color','w','markerfacecolor','w');
               leg(c) = plot(scatterdata.obs{e}(c,:),scatterdata.model{e}(c,:),'o','markersize',4,'markerfacecolor',disp.colors(c,:),'markeredgecolor','none');
            end
         end

      % pretty it up
         figureDefaults
         set(gca,'xlim',xlim,'ylim',xlim,'ticklength',[0.025 0.05],'xtick',xtick,'ytick',xtick);
         xlabel('Measured d^{\prime}','fontname','arial','fontsize',10);
         ylabel('Predicted d^{\prime}','fontname','arial','fontsize',10);
         legend(leg,{'neutral' 'central'},'fontname','arial','fontsize',6,'location','northwest');



   % SF preference of attention
      subplot(3,2,5)
      xlim              = [0 20];
      xtick             = 0:4:20;
      ylim              = [0.25 8];
      ytick             = [0.25 0.5 1 2 4 8];
      targfreq_colors   = [0 0 0; 150 150 150]./255;

      % add in the target SFs
         freqtypes = {'fine' 'coarse'};
         for f = 1:numel(freqtypes)
            center = boot.targfreq.(freqtypes{f}).center;
            leg(f) = line(xlim,[center center],'color',targfreq_colors(f,:),'linewidth',2); hold on

            % error region
            x = attndisp.ecc;
            y = zeros(2,numel(x))+boot.targfreq.(freqtypes{f}).ci;
            x = [x, fliplr(x)];
            y = [y(1,:) fliplr(y(2,:))];
            g = fill(x,y,targfreq_colors(f,:));
            set(g,'facecolor',targfreq_colors(f,:),'facealpha',0.25,'edgecolor','none'); hold on
         end

      % add attention SF preference
         leg(3) = semilogy(attndisp.ecc,attndisp.center,'-','linewidth',2.5,'color',disp.colors(2,:)); hold on

         % error region
         x = [attndisp.ecc, fliplr(attndisp.ecc)];
         y = [attndisp.ci(1,:) fliplr(attndisp.ci(2,:))];
         g = fill(x,y,targfreq_colors(f,:));
         set(g,'facecolor',disp.colors(c,:),'facealpha',0.25,'edgecolor','none'); hold on

      % pretty up figure
         figureDefaults
         set(gca,'yscale','log','ticklength',[0.025 0.05],'ytick',ytick,'xlim',xlim,'ylim',ylim);
         xlabel('Eccentricity (^o)','fontname','arial','fontsize',10);
         ylabel('Spatial frequency (cpd)','fontname','arial','fontsize',10);
         legend(leg,{'f_{fine}' 'f_{coarse}' 'f_{narrow}'},'location','northeast','fontname','arial','fontsize',8);



   % Gain overlap on target
      subplot(3,2,6)
      clear leg
      xlim              = [0 20];
      xtick             = 0:4:20;
      ylim              = [0 1];
      ytick             = [0 0.5 1];
      targfreq_colors   = [0 0 0; 150 150 150]./255;

      % plot curves
         freqtypes = {'fine' 'coarse'};
         for f = 1:numel(freqtypes)
            % error region
            x = gain_on_targ.ecc;
            y = gain_on_targ.(freqtypes{f}).ci;
            x = [x, fliplr(x)];
            y = [y(1,:) fliplr(y(2,:))];
            g = fill(x,y,targfreq_colors(f,:));
            set(g,'facecolor',targfreq_colors(f,:),'facealpha',0.25,'edgecolor','none'); hold on

            % center
            center = gain_on_targ.(freqtypes{f}).center;
            leg(f) = plot(gain_on_targ.ecc,center,'-','linewidth',2,'color',targfreq_colors(f,:)); hold on
         end

      % pretty up figure
         figureDefaults
         set(gca,'xlim',xlim,'ylim',ylim,'ticklength',[0.025 0.05],'xtick',xtick,'ytick',ytick);
         legend(leg,freqtypes,'fontname','arial','fontsize',8,'location','southeast');
         xlabel('Eccentricity (^o)','fontname','arial','fontsize',10); 
         ylabel('Gain on target','fontname','arial','fontsize',10); 


   % Save figure
      saveas(gcf,'./figure6.pdf');
