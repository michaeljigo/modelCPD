% train svm and classify test instances
% assumes inputs are already a matrix of instances

function classRes = doSVM(train,test,preprocess)
% only use a linear kernel to classify
kernelfun = 'linear';

trainInstances = train;
testInstances = test;

% preprocess instances
% set up pSettings structure to work with edited preprocessing coee
% (preprocessInst)
pSettings.demean = 0;
pSettings.linearScale0_1 = 0;
pSettings.zscore = 0;
if ~ieNotDefined('preprocess')
    switch preprocess
        case 'demean=1'
            pSettings.demean = 1;
        case 'zscore=1'
            pSettings.zscore = 1;
            pSettings.demean = 1;
        case 'linearScale0_1'
            pSettings.linearScale0_1 = 1;
    end
    [trainInstances, ~] = preprocessInst(trainInstances,pSettings); % right now, preprocessing is kept separate for train and test
    [testInstances, ~] = preprocessInst(testInstances,pSettings);
else
    [trainInstances, pSettings] = MJ_preprocessInstances(trainInstances,s.whichPreProcess);
    [testInstances, ~] = MJ_preprocessInstances(testInstances,pSettings);
end

% make instances libsvm-compatible
trainGrps = []; testGrps = []; numClasses = length(trainInstances);
for iClass = 1:numClasses
    trainGrps = [trainGrps; iClass*ones(size(trainInstances{iClass},1),1)];
    testGrps = [testGrps; iClass*ones(size(testInstances{iClass},1),1)];
end

% calculate c
bestcv = 0; bestc = 1; bestg = 2;
if strcmp(kernelfun,'linear')
    for log2c = -8:2:2
        opt = ['-t 0 -v ',num2str(6),' -c ',num2str(2^log2c),' -q'];
        cv = svmtrain(cell2mat(trainInstances'), trainGrps, opt);
        if cv > bestcv
            bestcv = cv; bestc = 2^log2c; bestg = nan;
        end
    end
end

% train a svm with best C from cross-validation search
if strcmp(kernelfun,'linear')
    opt = ['-t 0 -c ',num2str(bestc), ' -q'];
end
model = svmtrain(cell2mat(trainInstances'),trainGrps, opt);

% use this model to classify test instances
[predictedClass,~,~,] = svmpredict(testGrps,cell2mat(testInstances'),model);
classRes(iROI) = sum(predictedClass==testGrps)/length(testGrps);