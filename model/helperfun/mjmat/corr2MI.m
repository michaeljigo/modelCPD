function mi = corr2mi(r)
% transform correlation coefficients to mutual information

mi = -0.5.*log(1-r.^2);