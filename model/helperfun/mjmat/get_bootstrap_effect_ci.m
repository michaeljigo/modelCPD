% Purpose:  Compute bootstrapped confidence intervals of the main effects and interactions among the factors in the inputted data matrix.
%           The first dimension of matrix must contain separate subjects.
%
% Syntax:   [ci pval raw] = get_bootstrap_effect_ci(datamat,nboot,factornames,levelvals)
%
% By:       Michael Jigo
% Edited:   07.02.21

function [ci pval raw] = get_bootstrap_effect_ci(datamat,nboot,factornames,levelvals)
   
   % defaults
      if nargin<2
         nboot = 1e3;
      end
      if nargin<3
         factornames = [];
      end
      if nargin<4
         levelvals = [];
      end


   % main effects
      n_factors = numel(size(datamat))-1;
      for n = 1:n_factors
         if isempty(factornames)
            factorname  = sprintf('factor%i',n);
         else
            factorname  = factornames{n};
         end
         factoridx   = n+1;
         n_levels    = size(datamat,factoridx);

               
         % average across other factors
            avgidx      = setdiff(1:numel(size(datamat)),[1 factoridx]);
            if numel(avgidx>1)
               factoravg = datamat;
               obs       = datamat;
               for a = 1:numel(avgidx)
                  factoravg = nanmean(factoravg,avgidx(a));
                  obs       = nanmean(obs,avgidx(a));
               end
            else
               factoravg   = nanmean(datamat,avgidx);
               obs         = nanmean(obs,avgidx);
            end
            factoravg = squeeze(factoravg);


            % store factor average
            raw.(factorname) = nanmean(factoravg,1);
         

         % get x-values for the levels in the current factor
            if isempty(levelvals) || any(isnan(levelvals{n}))
               xval = 0:size(factoravg,2)-1;
            else
               xval = levelvals{n};
            end
               


         % compute confidence interval of contrasts between levels of factor
         teststat   = [];
         contrasts   = combnk(1:n_levels,2); % contrasts to-be-computed
         switch size(contrasts,1)
            case 1 % difference between 2 levels in a factor
               for c = 1:size(contrasts,1)
                  teststat(c,:) = factoravg(:,contrasts(c,1))-factoravg(:,contrasts(c,2)); 
               end
               % collapse across contrasts
                  teststat   = teststat(:);
            
               % get 95% confidence interval of group-average test statistic
                  [ci.(factorname) bootdist] = get_bootstrap_ci(teststat,0.025);
                  ci.(factorname)            = ci.(factorname)';
               
               % compute two-tailed p-value...only doing this is reviewer asks for p-value
                  lefttail                   = sum(bootdist<0)./numel(bootdist);
                  righttail                  = sum(bootdist>0)./numel(bootdist);
                  pval.(factorname)          = min(lefttail,righttail)*2;

            otherwise % fit a slope across >= 3 levels in a factor
               % bootstrap group-average estimates
                  [~, bootdist] = get_bootstrap_ci(factoravg,0.025);
                  teststat = [];
                  for b = 1:size(bootdist,1)
                     % remove nans from estimation
                        yval        = bootdist(b,:);
                        nanidx      = ~isnan(yval);

                     % fit linear function to data
                        lineparams  = polyfit(xval(nanidx),yval(nanidx),1);
                        teststat(b) = lineparams(1); % extract slope
                  end

               % get 95% confidence interval of group-average test statistic
                  ci.(factorname)      = quantile(teststat,[0.025 1-0.025]);

               % compute two-tailed p-value...only doing this is reviewer asks for p-value
                  lefttail                   = sum(teststat<0)./numel(bootdist);
                  righttail                  = sum(teststat>0)./numel(bootdist);
                  pval.(factorname)          = min(lefttail,righttail)*2;
               end
         end


   %% interaction
   switch n_factors
      case {2 3 4} % two, three, and four-way interactions
         % get # of two-way interactions
            interactions = combnk(1:n_factors,2);

         % compute interactions
         for ii = 1:size(interactions,1)
            % factors of interest
               factoridx   = interactions(ii,:);
               factorname  = sprintf('%s_x_%s',factornames{factoridx(1)},factornames{factoridx(2)});
            
            % average across non-interest factors
               if n_factors>2
                  nointerest  = setdiff(1:n_factors,factoridx)+1;
                  factoravg   = squeeze(nanmean(datamat,nointerest));
               else
                  factoravg   =  datamat;
               end

            % create separate bootstrap samples 
               [~,bootdist] = get_bootstrap_ci(factoravg);

            % fit separate slope to each level within 1st dimension
               xval = levelvals{factoridx(end)}; % select xvals of 2nd dimension for fits
               for b = 1:size(bootdist,1)
                  slope = [];
                  for d = 1:size(bootdist,2)
                     % remove nans from estimation
                        yval        = squeeze(bootdist(b,d,:))';
                        nanidx      = ~isnan(yval);

                     % fit linear function to data
                        lineparams  = polyfit(xval(nanidx),yval(nanidx),1);
                        slope(d)    = lineparams(1);
                  end
                  % test statistic is difference between slopes
                     teststat(b) = mean(diff(slope));

                  % store slopes
                     raw.(factorname)(b,:) = slope;
               end
               raw.(factorname) = nanmean(raw.(factorname),1);

            % get confidence intervals
               ci.(factorname) = quantile(teststat,[0.025 1-0.025]);
               
            % compute two-tailed p-value...only doing this is reviewer asks for p-value
               lefttail                   = sum(teststat<0)./numel(bootdist);
               righttail                  = sum(teststat>0)./numel(bootdist);
               pval.(factorname)          = min(lefttail,righttail)*2;
         end

         % add 3-way interactions
            % # of 3-way interactions
               interactions = combnk(1:n_factors,3);

            for ii = 1:size(interactions,1)
               factoridx   = interactions(ii,:);
               factorname  = sprintf('%s_x_%s_x_%s',factornames{factoridx(1)},factornames{factoridx(2)},factornames{factoridx(3)});
            
               % average across non-interest factor
                  if n_factors>3
                     nointerest  = setdiff(1:n_factors,factoridx)+1;
                     factoravg   = squeeze(nanmean(datamat,nointerest));
                  else
                     factoravg   =  datamat;
                  end

               % compute the difference across dimension that has 2 levels 
                  nlevels     = size(factoravg);
                  factordiff  = squeeze(mean(diff(factoravg,[],2),2));

               % now analyze as if a two-way interaction, as above
                  % create separate bootstrap samples of group-average within factors
                     [~,bootdist] = get_bootstrap_ci(factordiff);
   
                  % fit separate slope to each level within 1st dimension
                     xval = levelvals{factoridx(end)}; % select xvals of 2nd dimension for fits
                     for b = 1:size(bootdist,1)
                        slope = [];
                        for d = 1:size(bootdist,2)
                           % remove nans from estimation
                              yval        = squeeze(bootdist(b,d,:))';
                              nanidx      = ~isnan(yval);

                           % fit linear function to data
                              lineparams  = polyfit(xval(nanidx),yval(nanidx),1);
                              slope(d)    = lineparams(1);
                        end
                        % test statistic is average difference across slopes
                           teststat(b) = mean(diff(slope));
      
                        % store slopes
                           raw.(factorname)(b,:) = slope;
                     end
                     raw.(factorname) = nanmean(raw.(factorname),1);
   
                  % get confidence intervals
                     ci.(factorname) = quantile(teststat,[0.025 1-0.025]);
                  
                  % compute two-tailed p-value...only doing this is reviewer asks for p-value
                     lefttail                   = sum(teststat<0)./numel(bootdist);
                     righttail                  = sum(teststat>0)./numel(bootdist);
                     pval.(factorname)          = min(lefttail,righttail)*2;
               end


            % add 4-way interaction
               if n_factors<4
                  return
               end
               interactions = combnk(1:n_factors,4);
               factorname  = sprintf('%s_x_%s_x_%s_x_%s',factornames{1},factornames{2},factornames{3},factornames{4});

               % take difference from 2 of those dimensions, squeeze, then analyze as if a 2-way interaction
                  factordiff  = mean(diff(datamat,[],2),2); 
                  factordiff  = squeeze(mean(diff(factordiff,[],3),3));
               
               % now analyze as if a two-way interaction, as above
                  % create separate bootstrap samples of group-average within factors
                     [~,bootdist] = get_bootstrap_ci(factordiff);
   
                  % fit separate slope to each level within 1st dimension
                     for b = 1:size(bootdist,1)
                        slope = [];
                        for d = 1:size(bootdist,2)
                           lineparams  = polyfit(xval,squeeze(bootdist(b,d,:))',1);
                           slope(d)    = lineparams(1);
                        end
                        % test statistic is difference between slopes
                           teststat(b) = mean(diff(slope));
      
                        % store slopes
                           raw.(factorname)(b,:) = slope;
                     end
                     raw.(factorname) = nanmean(raw.(factorname),1);
   
                  % get confidence intervals
                     ci.(factorname) = quantile(teststat,[0.025 1-0.025]);
                  
                  % compute two-tailed p-value...only doing this is reviewer asks for p-value
                     lefttail                   = sum(teststat<0)./numel(bootdist);
                     righttail                  = sum(teststat>0)./numel(bootdist);
                     pval.(factorname)          = min(lefttail,righttail)*2;

      otherwise
         error('ERROR: Only two and three-way interactions are supported.');
   end
