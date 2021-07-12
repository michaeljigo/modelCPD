% Purpose:  This function will align the attention effects with respect to the peak SF in the neutral condition.
%           In addition, the function will place attention effects in a matrix whose size is determined by the 
%           full range of octaves spanned in the entire dataset.

function [aligned fullOct] = alignAttn(fullOct,subjOct,subjAttn,bin)

if ~exist('bin','var') || isempty(bin)
   bin = 1;
end

if numel(size(subjAttn))==2
   subjAttn = shiftdim(subjAttn,-1);
   subjOct = shiftdim(subjOct,-1);
end

% bin octaves
if bin
   binOct = round(subjOct*1/bin)./(1/bin);
   fullOct = unique(round(fullOct*1/bin)./(1/bin));
   %binOct = subjOct;
   %fullOct = fullOct;
else
   binOct = subjOct;
end

% initialize matrix to hold aligned attention effect 
aligned = nan(size(subjAttn,1),size(subjAttn,3),numel(fullOct));

% loop through subjects (or resampled data) and eccentricity and align attention effect
newAttn = subjAttn;
for s = 1:size(subjAttn,1)
   for e = 1:size(subjAttn,3)
      % check for octaves that were binned together, for these values average the binned values together
      for point = 1:numel(subjOct(s,:,e))
         if bin
            %% bin to range around point
            binRange = [-bin/2 bin/2]+subjOct(s,point,e);
            binIdx = subjOct(s,:,e)>=min(binRange) & subjOct(s,:,e)<=max(binRange);
            if sum(binIdx)>1
               binAvg = mean(subjAttn(s,binIdx,e));
               newAttn(s,binIdx,e) = binAvg;
            end

            %% bin to nearest octave (upward); one value cannot belong to two bins
            %binIdx = binOct(s,point,e)==binOct(s,:,e);
            %if sum(binIdx)>1
               %newAttn(s,point,e) = mean(subjAttn(s,binIdx,e));
            %end
         end
      end
      subjAttn(s,:,e) = newAttn(s,:,e);

      % do alignment
      [preAlign,postAlign] = ismember(squeeze(binOct(s,:,e)),fullOct);
      aligned(s,e,postAlign) = subjAttn(s,preAlign,e);
   end
end
