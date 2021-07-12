% calculates the lengths of the inputted vectors and their angular
% separation

function [len, ang] = vecLenAngle(v1,v2,dispOut)

if ieNotDefined('dispOut')
    dispOut = 0;
end

% make any row vectors into column vectors
M = {v1 v2};
[row, ~] = cellfun(@size,M);
row = find(row==1); % index row vectors
M(row) = cellfun(@transpose,M(row),'UniformOutput',false); % transpose row vectors
try
    len = sqrt(sum(cell2mat(M).^2)); % calculate lengths
    ang = acosd(dot(M{1},M{2})/prod(len)); % angle difference in degrees
    if dispOut
        fprintf('vector lengths = [%.2f %.2f]; \n angle = %.2f degrees \n',len,ang);
    end
catch
    error('Input vectors should be the same length');
end

end