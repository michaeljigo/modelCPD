% Purpose:  Do a permutation F-test with ANOVA.
%           Algorithm follows Manly, 2007.
%           Manly, B. F. J. (2007) Randomization, Bootstrap, and Monte Carlo Methods 
%           in Biology (3rd ed.), London: Chapman & Hall.
%           
%           Only set up for repeated-measures designs.

function out = permutation_ftest(datamat,factors,varargin)

   nperm = 1e3;   % default # permutation samples
   verbose = 0;   % print iterations
   if nargin>2 
      nperm = varargin{1};
      if nargin>3
         verbose = varargin{2};
      end
   end

   if verbose
      init = tic;
   end
   rng('shuffle')
   nsubj = size(datamat,1);
   for n = 1:nperm
      % collapse across factors
      permmat = datamat(:,:);

      % randomize labels
      randidx = randi(size(permmat,2),[1 size(permmat,2)]);
      for s = 1:nsubj
         permmat(s,:) = permmat(s,randidx);
      end
      %permmat = permmat(:,randidx);

      % reshape
      permmat = reshape(permmat,size(datamat));

      % repeated-measures ANOVA
      tbl = simple_mixed_anova(permmat,[],factors);

      % store f-ratio 
      f(n,:) = tbl{3:2:end,4};

      if verbose
         fprintf('%i ',n);
      end
   end
   if verbose
      time = toc(init);
      fprintf('\n Took %.2f seconds\n',time);
   end
   % store ANOVA test labels
   labels = tbl.Properties.RowNames;
   labels = labels(3:2:end);
   labels = cellfun(@(x) x(13:end),labels,'uniformoutput',0);


   % observed F-ratio
   tbl = simple_mixed_anova(datamat,[],factors);
   obs.tbl = tbl;
   obs.f = tbl{3:2:end,4}';
   
   % permutation-based p-val
   pval = sum(obs.f<f)./nperm;
   
   % format output structure
   out.null_fratio = f;
   out.tests = labels;
   out.pval = pval;
   out.obs = obs;
