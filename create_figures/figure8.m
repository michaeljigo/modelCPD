% Purpose:  Recreate Figure 8 showing fits to texture segmentation experiments.
%
% By:       Michael Jigo
% Edited:   07.08.21

function figure8
addpath(genpath('../../modelCPD_v4'));

%% Preprocess and format data before displaying
   % visualization parameters
      disp.ci           = 0.16;     % 68% confidence interval 
      disp.extrap       = 'pchip';  % extrapolation method
      disp.endo_colors  = [0 0 0; 202 0 32]./255; 
      disp.exo_colors   = [0 0 0; 5 113 176]./255; 


   % load bootstrap samples for exogenous and endogenous attention fits
      attntype = {'exo' 'endo'};
      for a = 1:numel(attntype)
         boot.(attntype{a}) = combine_bootstrap_iterations('model_variant','main_model','attn_type',attntype{a});
      end


   % get indices of experiments to display
      exps     = {'tc02_lvm' 'clh06_baseline' 'yc08_cuesize1' 'ymc08_exp2' 'ymc08_exp1' 'bc17_exp1'};
      expfiles = {'tc02' 'clh06' 'yc08' 'ymc08_exp2' 'ymc08_exp1' 'bc17'};


   % interpolate median and 68% confidence interval of model responses
      for e = 1:numel(exps)
         % find the experiment of interest and create color matrix
            if ismember(exps{e},boot.exo.sim_params.exp_list)
               data     = boot.exo;
               disp.colors(:,:,e) = disp.exo_colors;
            else
               data  = boot.endo;
               disp.colors(:,:,e) = disp.endo_colors;
            end
            expidx   = find(ismember(data.sim_params.exp_list,exps{e}));

         expdata              = data.exp_data(expidx);
         modelresp            = data.model_resp(expidx).dprime;
         medianresp           = squeeze(nanmedian(modelresp,3));
         modelci              = squeeze(quantile(modelresp,[disp.ci 1-disp.ci],3));
         modeldisp(e).ecc     = linspace(0,max(expdata.ecc),20);

         for c = 1:size(expdata.dprime,1)
            modeldisp(e).center(c,:)  = interp1(expdata.ecc,medianresp(c,:),modeldisp(e).ecc,disp.extrap);
            modeldisp(e).ci(c,:,:)    = interp1(expdata.ecc,squeeze(modelci(c,:,:)),modeldisp(e).ecc,disp.extrap)';
         end

      % add in error for each eccentricity and validity condition
         err = load(sprintf('../data/behavior/%s.mat',expfiles{e}));
         sem                          = err.data.sem;

         % compute upper and lower bound based on SEM
         y = data.exp_data(expidx).dprime;
         modeldisp(e).err.lb          = y-sem;
         modeldisp(e).err.ub          = y+sem;

      % store observed data
         modeldisp(e).obs.ecc         = expdata.ecc;
         modeldisp(e).obs.dprime      = expdata.dprime;
      end


%% Plot
      figure('name','Figure 8','position',[360 61 428 557]);
      xlims       = [0 12; 0 24; 0 12; 0 12; 0 12; 0 8];
      xtick       = 0:4:24; 
      ylim        = [0.75 2.25; 0 2.25; 0.5 2; 0.4 1.6; 0.4 2; 0 3];
      ytick       = 0:0.5:5;

      for e = 1:numel(modeldisp)
         subplot(3,2,e);
         for c = 1:size(modeldisp(e).center,1)
            % draw model fits (error)
               x = [modeldisp(e).ecc, fliplr(modeldisp(e).ecc)];
               y = [squeeze(modeldisp(e).ci(c,1,:))' fliplr(squeeze(modeldisp(e).ci(c,2,:))')];
               f = fill(x,y,disp.colors(c,:,e));
               set(f,'facecolor',disp.colors(c,:,e),'facealpha',0.25,'edgecolor','none'); hold on

            % draw model fits (median)
               plot(modeldisp(e).ecc,modeldisp(e).center(c,:),'-','linewidth',3,'color',disp.colors(c,:,e));
   
            % draw observed data
               x = modeldisp(e).obs.ecc;
               y = modeldisp(e).obs.dprime(c,:);
               plot(x,y,'o','markersize',5,'color','w','markerfacecolor','w');
               plot(x,y,'o','markersize',4,'markerfacecolor',disp.colors(c,:,e),'markeredgecolor','none');

            % draw error bars
               for ec = 1:size(y,2)
                  lb = modeldisp(e).err.lb(c,ec);
                  ub = modeldisp(e).err.ub(c,ec);
                  line([x(ec) x(ec)],[lb ub],'color',[disp.colors(c,:,e) 0.5],'linewidth',2); 
               end
         end
         % pretty up subplot
            figureDefaults
            set(gca,'xlim',xlims(e,:),'xtick',xtick,'ylim',ylim(e,:),'ytick',ytick,'ticklength',[0.025 0.05]);
            if e==(numel(modeldisp)-1)
               xlabel('Eccentricity (^o)','fontname','arial','fontsize',10); 
               ylabel('Performance (d^{\prime})','fontname','arial','fontsize',10); 
            end
            title(expfiles{e},'fontname','arial','fontsize',10);
      end
