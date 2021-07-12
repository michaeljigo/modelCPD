% Purpose:  This function will change the range of a vector or matrix to span a specified range.
% By:       Michael Jigo

function out = change_range(in,new_min,new_max,org_min,org_max);

if nargin==3
   % get original range from inputted matrix
   org_min = min(in(:));
   org_max = max(in(:));
else
   % use user-inputted min/max as the original range
end


% change to new range
out = (new_max-new_min)*((in-org_min)./(org_max-org_min))+new_min;
