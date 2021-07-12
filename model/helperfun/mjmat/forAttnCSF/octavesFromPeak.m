% Purpose:  This function will compute the number of octaves on either side of the peak SF during the neutral condition.
%           This computation is a pre-requisite for computing the aligned attention effect.

function [octaves peakSF] = octavesFromPeak(thresh,sfVal,fit_peakSF);

if ~exist('useFit','var')
   useFit = 1;
end   

% get peak SF
if ~exist('fit_peakSF','var') || isempty(fit_peakSF)
   % get peak based directly from data (non-parameteric)
   [~,peakSF] = min(thresh,[],2);
   peakSF = squeeze(peakSF(1,:,:));
   peakSF = sfVal(peakSF);
else
   peakSF = fit_peakSF;
end

% compute octaves from peak
for e = 1:length(peakSF)
   octaves(:,e) = log2(sfVal./peakSF(e));
end 
octaves = round(octaves*1e2)./1e2;
