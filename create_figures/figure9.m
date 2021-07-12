% Purpose:  Recreate Figure 9 showing fits to acuity and texture segmentation experiments
%
% By:       Michael Jigo
% Edited:   07.08.21

function figure9
addpath(genpath('../../modelCPD_v4'));

   % display parameters
      colors = [5 113 176; 202 0 32]./255;


   % load best-fitting parameters for
      attntypes = {'exo' 'endo'};
      % acuity
      acuitydir = '../generalize/MontagnaPestilliCarrasco2009/data/fitted_parameters/';
      for a = 1:numel(attntypes)
         tmp  = load([acuitydir,attntypes{a},'.mat']);
         acuity.(attntypes{a}) = tmp.out;

         % compute group-average and SEM of behavioral acuity thresholds
         acuity.(attntypes{a}).avg = mean(acuity.(attntypes{a}).data.thresh);
         acuity.(attntypes{a}).sem = std(acuity.(attntypes{a}).data.thresh)./sqrt(size(acuity.(attntypes{a}).data.thresh,1));
      end
      
      % contrast sensitivity
      csmodel_dir = '../generalize/JigoCarrasco2020/data/fitted_parameters/';
      csbehav     = load('../generalize/JigoCarrasco2020/data/behavior/exp1.mat','data');
      csbehav     = csbehav.data;

      for a = 1:numel(attntypes)
         tmp = load([csmodel_dir,'exp1_',attntypes{a},'.mat']);
         cs.(attntypes{a})  = tmp.out;

         % compute group-average and SEM of behavioral cueing effects
            sensitivity = 1./csbehav(a).crf.thresh; % convert from thresholds to sensitivity
            cue_effect  = squeeze(sensitivity(:,:,:,2)./sensitivity(:,:,:,1));
            cs.(attntypes{a}).avg = squeeze(mean(cue_effect));
            err_effect  = 1./csbehav(a).crf.raw_thresh; err_effect = err_effect(:,:,:,2)./err_effect(:,:,:,1);
            cs.(attntypes{a}).sem = squeeze(std(err_effect)./sqrt(size(err_effect,1)));

         % compute Neutral peak SF 
            cs.(attntypes{a}).neutSF = exp(mean(log(cs.(attntypes{a}).data.csf.peakSF),1));
      end


   % Display acuity
      figure('name','Figure 9','position',[360 78 235 540]);
      xval = [1 1.4];
      for a = 1:numel(attntypes)
         subplot(5,2,a)

         % draw the average
            avg = acuity.(attntypes{a}).avg;
            sem = acuity.(attntypes{a}).sem;
            for ii = 1:numel(xval)
               leg(1) = bar(xval(ii),avg(ii),0.3,'edgecolor','none','facecolor',colors(a,:)); hold on
               set(leg(1),'facealpha',0.5);
            end

         % draw SEM
            errorbar(xval,avg,sem,'linestyle','none','capsize',0,'color',[0.5 0.5 0.5],'linewidth',2);

         % add in model predictions
            plot(xval,acuity.(attntypes{a}).model.thresh,'o','markersize',6,'markerfacecolor','w','markeredgecolor','none'); hold on
            leg(2) = plot(xval,acuity.(attntypes{a}).model.thresh,'o','markersize',5,'markerfacecolor',colors(a,:),'markeredgecolor','none'); hold on

         % pretty up plot
            figureDefaults
            ylabel('Threshold (arc min)','fontname','arial','fontsize',10);
            set(gca,'xlim',[min(xval)-0.3 max(xval)+0.3],'xtick',xval,'xticklabel',{'Neutral' 'Valid'},'ylim',[0 15],'ytick',0:5:50,'fontsize',6);
            title(attntypes{a},'fontname','arial','fontsize',10);
      end


   % Display contrast sensitivity
      subplotidx = [3:2:9;4:2:10];
      interp_freq = logspace(log10(0.5),log10(8),12); % 11

      for a = 1:numel(attntypes)
         for e = 1:numel(cs.exo.data.ecc)
            subplot(5,2,subplotidx(a,e));

            % draw line at 1
               line([0.125 25],[1 1],'color',[0.5 0.5 0.5],'linewidth',2); hold on
      
            % draw vertical lines for Neutral peak SF
               peakSF = cs.(attntypes{a}).neutSF(e);
               line([peakSF peakSF],[1 1.4],'color','k','linewidth',2); hold on

            % draw fits from Jigo & Carrasco 2020
               jc20fit.freq = csbehav(a).groupavgFit.freq;
               jc20fit.attn = csbehav(a).groupavgFit.attn(e,:);
               semilogx(jc20fit.freq,jc20fit.attn,'-','linewidth',4,'color',[colors(a,:) 0.2]); hold on

            % draw fits by model
               x = cs.(attntypes{a}).finer_SF_sample.freq;
               y = cs.(attntypes{a}).finer_SF_sample.model.attn(e,:);
               interp_fit = interp1(x,y,interp_freq,'pchip');
               semilogx(interp_freq,interp_fit,'-','linewidth',3,'color',colors(a,:)); hold on

            % draw subject avg 
               freq   = csbehav(a).freq;
               effect = cs.(attntypes{a}).avg(e,:);
               err    = cs.(attntypes{a}).sem(e,:);
               semilogx(freq,effect,'o','markersize',6,'color','w'); hold on
               semilogx(freq,effect,'o','markersize',5,'markerfacecolor',colors(a,:),'markeredgecolor','none'); hold on
               errorbar(freq,effect,err,'.','markersize',0.01,'color',colors(a,:),'capsize',0,'linestyle','none','linewidth',1.5);

            % pretty up subplot
               figureDefaults
               if e==4 && a==1
                  xlabel('Spatial frequency (cpd)','fontname','arial','fontsize',6);
                  ylabel('Attention effect','fontname','arial','fontsize',6);
               end
               set(gca,'xlim',[0.3 16],'xtick',2.^(-1:4),'ytick',1:0.2:5,'xscale','log','yscale','linear','ylim',[0.98 1.3],'ytick',[1:0.15:2],'fontsize',6);
         end
      end
