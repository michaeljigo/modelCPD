% Purpose:  Correct dprime using Hautas adjustment (1995).
%           Input is the output from condParser with the target-presence variable in the first dimension.

function dprime = hautas_adjustment(dprime)
                  
   fa_hits = cellfun(@nansum, dprime.raw)+0.5;
   noise_sig = cellfun(@numel, dprime.raw)+1;

   % find zero-trial conditions
   numtrials = cellfun(@numel, dprime.raw);
   zerotrial = squeeze(any(numtrials==0,1));

   dprime.perf = fa_hits./noise_sig;
   dprime.perf = squeeze(diff(norminv(dprime.perf),[],1));
   dprime.perf(zerotrial) = nan;
   dprime.numtrials = squeeze(sum(numtrials));
