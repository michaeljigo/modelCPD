% build a data matrix for 1-way repeated measures ANOVA. The input should
% be organized as follows: data{subject}(IV1)

% 03/22/17: Added functionality to use a matrix as data input. The matrix
% input should be organized as follows: data(subject,IV1)

function x = buildrmaov1Matrix(data)

if iscell(data)
    nsubj = length(data);
    nIV1 = length(data{1});
else
    if length(size(data))~=2
        error('Input matrix should be nx2'];
    end
    [nsubj, nIV1] = size(data);
end
matrixSize = nIV1*nsubj;
x = []; IV1 = 1;

while ~all(size(x)==[matrixSize 3])
    for iSubj = 1:nsubj
        if iscell(data)
            x = vertcat(x,[data{iSubj}(IV1),IV1,iSubj]);
        else
            x = vertcat(x,[data(iSubj,IV1),IV1,iSubj]);
        end
    end
    IV1 = IV1+1;
end