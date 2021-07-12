function m = nangeomean(x,dim)
% NANGEOMEAN Geometric mean that ignores nans
% see help of geomean for more info.

if nargin < 2 || isempty(dim)
    % Figure out which dimension sum will work along.
    dim = find(size(x) ~= 1, 1);
    if isempty(dim), dim = 1; end
end
m = geomean(x(~isnan(x)),dim);
