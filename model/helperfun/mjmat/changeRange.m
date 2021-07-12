% usage:    out = changeRange(in,newRange)
% by:       Michael Jigo
% date:     03/07/18
% purpose:  Adjust range of vector, matrix, or cell array to arbitrary range. 
%           If input is a matrix or cell array, all elements or cells will be 
%           normalized to the same same range.
%
% INPUTS:
% in        input vector, matrix, or cell array
% newRange  nx2 vector corresponding to the min and max (default: [0 1])

function normOut = changeRange(in,newRange)

% make single vector/matrix into a cell to make the function work
if ~iscell(in)
   in = {in};
end

if ~exist('newRange','var')
   newRange = [0 1];
else
   newRange = sort(newRange);
end
newMin = newRange(1); 
newMax = newRange(2);

% normalize each cell
for i = 1:length(in)
   minI = min(in{i}(:));
   maxI = max(in{i}(:));
   I = in{i};
   normOut{i} = (newMax-newMin)*(I-minI)./(maxI-minI)+newMin;
end
if i==1
   normOut = cell2mat(normOut);
end
