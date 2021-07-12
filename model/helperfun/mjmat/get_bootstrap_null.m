% Purpose:  Compute bootstrapped null distribution, as in Efron 1979 (Bootstrap methods, Another look at the jackknife)
%           The first dimension of data matrix must contain separate subjects. Other dimensions contain separate independent factors for analysis.
%
%           The full matrix will be shuffled and resampled with replacement. 
%           Contrasts will be taken to build null distributions of main effects and interactions among factors.
%
%           Currently works for 2-factor designs only.
%        
% By:       Michael Jigo
% Edited:   07.02.21

function null = get_bootstrap_null(datamat,nboot,factornames)

   % defaults
      if nargin<2
         nboot = 1e3;
      end
      if nargin<3
         factornames = [];
      end

   % sample, with replacement, and shuffle group-average measures
      rng('shuffle');
      alldata     = datamat(:);
      shuffleidx  = randi(numel(alldata),[nboot numel(alldata)]);
      for n = 1:nboot
         null(n,:) = alldata(shuffleidx(n,:));
      end

     
   % reshape into original form
      matsize = size(datamat);
      null = reshape(null,[nboot matsize]);

   
   % form group-level null distribution
      null = squeeze(nanmean(null,2));

   

   % create null distribution of main effects
      n_factors = numel(size(null))-1;
      for n = 1:n_factors
         if isempty(factornames)
            factorname  = sprintf('factor%i',n);
         else
            factorname  = factornames{n};
         end
         factoridx   = n+1;
         n_levels    = size(null,factoridx);


         % average across other factors
            avgidx      = setdiff(1:numel(size(null)),[1 factoridx]);
            if numel(avgidx>1)
               factoravg = null;
               obs       = datamat;
               for a = 1:numel(avgidx)
                  factoravg = nanmean(factoravg,avgidx(a));
                  obs       = nanmean(obs,avgidx(a));
               end
            else
               factoravg   = nanmean(null,avgidx);
               obs         = nanmean(obs,avgidx);
            end
            factoravg = squeeze(factoravg);
            obs       = squeeze(nanmean(obs,1));


         % compute contrasts between levels of factor
         contrasts   = combnk(1:n_levels,2); % contrasts to-be-computed
         for c = 1:size(contrasts,1)
            test.null.(factorname)(c,:) = factoravg(:,contrasts(c,1))-factoravg(:,contrasts(c,2)); 
            test.obs.(factorname)(c,:)  = obs(contrasts(c,1))-obs(contrasts(c,2));
         end
         test.null.(factorname)  = mean(test.null.(factorname),1);
         test.obs.(factorname)   = mean(test.obs.(factorname));
         
         % two-sided p-value
         lefttail                = sum(test.obs.(factorname)<test.null.(factorname))./nboot;
         righttail               = sum(test.obs.(factorname)>test.null.(factorname))./nboot;
         pval.(factorname)       = min(lefttail,righttail)*2;
      end
      keyboard
