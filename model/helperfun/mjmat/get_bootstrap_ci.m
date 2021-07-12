% Purpose:  Compute bootstrapped confidence intervals given an input matrix of data.
%           The first dimension of data matrix must contain separate subjects.
%
% Syntax:   [ci boot boot_subj] = get_bootstrap_ci(datamat,cival,nboot)
%
% By:       Michael Jigo
% Edited:   07.02.21

function [ci boot boot_subj] = get_bootstrap_ci(datamat,cival,nboot)

   % defaults
      if nargin<3
         nboot = 1e3;
      end
      if nargin<2
         cival = 0.16; % 68% confidence interval
      end
      if nargin<1
         error('Need to input a matrix');
      end

   % create resampling matrix
      rng('shuffle');
      nsubj = size(datamat,1);
      randidx = randi(nsubj,nboot,nsubj);
      for n = 1:nboot
         sample = datamat(randidx(n,:),:);
         boot(n,:) = nanmean(sample,1);
         boot_subj(n,:,:) = sample;
      end

   % reshape into original form
      matsize = size(datamat);
      matsize = matsize(2:end);
      boot = reshape(boot,[nboot matsize]);
      boot_subj = reshape(boot_subj,[nboot size(datamat)]);


   % compute confience intervals
      ci = quantile(boot,[cival 1-cival],1);
