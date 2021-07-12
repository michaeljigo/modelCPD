% written by: Michael Jigo
% build a data matrix for 2-way repeated measures ANOVA (rmaov2). The input should
% be organized as follows: data{subject}(IV1,IV2)

% 03/22/2017: Added in functionality to make anova matrix from a matrix
% input. Thus, if the data input is organized as follows:
% data(subject,IV1,IV2)
% The anova matrix can be made

function x = buildrmaov2Matrix(data)

if iscell(data)
    nsubj = length(data);
    nIV1 = size(data{1},1);
    nIV2 = size(data{1},2);
else
    if length(size(data))~=3
        error('Input matrix should be nx3');
    end
    [nsubj, nIV1, nIV2] = size(data);
end
matrixSize = nIV2*nIV1*nsubj;
x = []; IV2 = 1; IV1 = 1;

while ~all(size(x)==[matrixSize 4])
    for iSubj = 1:nsubj
        if iscell(data)
            x = vertcat(x,[data{iSubj}(IV1,IV2),IV1,IV2,iSubj]);
        else
            x = vertcat(x,[data(iSubj,IV1,IV2),IV1,IV2,iSubj]);
        end
    end
    if mod(IV2,nIV2)==0
        IV1 = IV1+1;
        IV2 = 0;
    end
    IV2 = IV2+1;
end