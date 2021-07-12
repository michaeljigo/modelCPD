% This script does a permutation test on the input 1xn or 2xn cell array, 
% each row is a condition and each column is a level in a condition. The
% test statistic is defined by the string input with c1 and c2
% representing condition 1 and condition 2.

% Outputs the percentile of the data test statistic and the null
% distribution

% selectSubset is an integer input that specifies, how many values in the
% data input will be placed in each condition. this input is necessary when
% a 1xn cell array is inputted as data.

% In the case where a 1xn cell array is the data input, the function will
% only output a null distribution. The p-value and t-statistic will be
% non-meaningful.

function [p, tStat, nullDist] = permutationTest(data,nPerm,testStat,selectSubset)

% make sure format of data is compatible with this code
correctFormatIdx = cellfun(@(x) size(x,1)==1,data);
correctFormat = cellfun(@(x) transpose(x),data(~correctFormatIdx),'UniformOutput',false);
data(~correctFormatIdx) = correctFormat;

% ensure the correct variables are inputted
if size(data,1)==1
    if ieNotDefined('selectSubset')
        selectSubset = 0;
        warning('Each condition will comprise of half of the inputted data.');
    end
else
    selectSubset = 0;
end

% initialize
rng('shuffle');
nullDist = nan(nPerm,size(data,2));

% loop through each level of the conditions
p = nan(1,size(data,2));

for c = 1:size(data,2)
    % collapse across cells 
    subset = cell2mat(data(:,c)');
    
    % determine the number of values in each condition.
    c2End = length(subset);
    if selectSubset
        c1Trials = selectSubset;
        c2End = c1Trials+selectSubset;
    elseif ~selectSubset
        c1Trials = round(length(subset)/2);
    else % use number of values in each cell as the number of trials
        % get number of trials for the conditions (c2Trials = length(nullMat)-c1Trials)
        c1Trials = size(data{1,c},2);
    end
    
    % calculate test statistic for the data
    c1 = subset(1:c1Trials); c2 = subset(c1Trials+1:c2End);
    tStat(c) = eval(testStat);
    
    % concatenate trials of both conditions
    nullMat = repmat(subset,nPerm,1);
    
    % create nPermxn matrix of trial number permutations in order to
    % quickly index all the desired permutations without a for loop
    [~, permMat] = sort(rand(nPerm,size(nullMat,2)),2);
    
    % change trial number index to linear indices that can index nullMat
    shuffleIdx = bsxfun(@plus,(permMat-1)*nPerm,(1:nPerm)');
    nullMat = nullMat(shuffleIdx);
    c1 = nullMat(:,1:c1Trials);
    c2 = nullMat(:,c1Trials+1:c2End);
    
    % make the null distribution of test statistics for this level
    nullDist(:,c) = eval(testStat);
    
    % get p-value of data test statistic
    p(c) = sum(tStat(c)>=nullDist(:,c))/length(nullDist(:,c)); 
end