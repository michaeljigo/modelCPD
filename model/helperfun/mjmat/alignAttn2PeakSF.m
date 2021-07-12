% Purpose:  This function will align the attention effects with respect to the peak SF in the neutral condition.
%           In addition, the function will place attention effects in a matrix whose size is determined by the 
%           full range of octaves spanned in the entire dataset.

function aligned = alignAttn(fullOct,subjOct,subjAttn)

if numel(size(subjAttn))==2
   subjAttn = shiftdim(subjAttn,-1);
   subjOct = shiftdim(subjOct,-1);
end

% initialize matrix to hold aligned attention effect 
aligned = nan(size(subjAttn,1),size(subjAttn,3),numel(fullOct));

% loop through subjects (or resampled data) and eccentricity and align attention effect
for s = 1:size(subjAttn,1)
   for e = 1:size(subjAttn,3)
      [preAlign,postAlign] = ismember(squeeze(subjOct(s,:,e)),fullOct);
      aligned(s,e,postAlign) = subjAttn(s,preAlign,e);
   end
end
