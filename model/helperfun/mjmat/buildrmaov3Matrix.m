% build a data matrix for 3-way repeated measures ANOVA. The input should
% be organized as follows: data{subject}(IV1,IV2,IV3)

% 03/22/17: Added functionality to use matrix as data input. Data input
% should be organied as follows: data(s,IV1,IV2,IV3)

function x = buildrmaov3Matrix(data)

if iscell(data)
    nSubj = length(data);
    nIV1 = size(data{1},1);
    nIV2 = size(data{1},2);
    nIV3 = size(data{1},3);
else
    if length(size(data))~=4
        error('Input matrix should be nx4');
    end
    [nSubj, nIV1, nIV2, nIV3] = size(data);
end
matrixSize = nSubj*nIV1*nIV2*nIV3;
x = []; IV1 = 1; IV2 = 1; IV3 = 1;

while ~all(size(x)==[matrixSize 5])
    for iSubj = 1:nSubj
        if iscell(data)
            x = vertcat(x,[data{iSubj}(IV1,IV2,IV3),IV1,IV2,IV3,iSubj]);
        else
            x = vertcat(x,[data(iSubj,IV1,IV2,IV3),IV1,IV2,IV3,iSubj]);
        end
    end
    if IV3==nIV3
        IV3 = 0;
        if IV2~=nIV2
            IV2 = IV2+1;
        else
            IV2 = 1;
            if IV1~=nIV1
                IV1 = IV1+1;
            end
        end
    end
    IV3 = IV3+1;
end