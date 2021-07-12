% build a data matrix for 4-way repeated measures ANOVA. The input should
% be organized as follows: data{subject}(IV1,IV2,IV3,IV4)

% 03/22/17: Added functionality to use matrix as data input. Data input
% should be organied as follows: data(s,IV1,IV2,IV3,IV4)

function x = buildrmaov4Matrix(data)

if iscell(data)
    nSubj = length(data);
    nIV1 = size(data{1},1);
    nIV2 = size(data{1},2);
    nIV3 = size(data{1},3);
    nIV4 = size(data{1},4);
else
    if length(size(data))~=5
        error('Input matrix should be nx5');
    end
    [nSubj, nIV1, nIV2, nIV3, nIV4] = size(data);
end
matrixSize = nSubj*nIV1*nIV2*nIV3*nIV4;
x = []; IV1 = 1; IV2 = 1; IV3 = 1; IV4 = 1;

while ~all(size(x)==[matrixSize 6])
    for iSubj = 1:nSubj
        if iscell(data)
            x = vertcat(x,[data{iSubj}(IV1,IV2,IV3,IV4),IV1,IV2,IV3,IV4,iSubj]);
        else
            x = vertcat(x,[data(iSubj,IV1,IV2,IV3,IV4),IV1,IV2,IV3,IV4,iSubj]);
        end
    end
    if IV4==nIV4
       IV4 = 0;
       if IV3~=nIV3
          IV3 = IV3+1;
       else
         IV3 = 1;
         if IV2~=nIV2
            IV2 = IV2+1;
         else
            IV2 = 1;
            if IV1~=nIV1
               IV1 = IV1+1;
            end
         end
       end
    end
    IV4 = IV4+1;
end

