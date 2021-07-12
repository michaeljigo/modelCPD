% identical to the original preprocessInstance code. i just added linear
% scaling as a preprocessing step --MJ

% preprocessInstances.m
%
%        $Id:$ 
%      usage: instances = preprocessInstances(instances)
%         by: justin gardner
%       date: 10/25/13
%    purpose: Runs various pre-processing steps on instances. You can:
%             remove the mean across all instances (regardless of type):
%               instances = preprocessInstances(instances,'demean=1');
%             z-score the instances across all instances
%               instances = preprocessInstances(instances,'zecore=1');
%             do pca and keep only components that account for a certain amount of variance
%             across instances. For example to keep components accounting for 60% of variance
%               [instances pSettings] = preprocessInstances(instances,'pca=0.6');
%             note that the variable pSettings will contain the amount of variance accounted 
%             for by each component pSettings.v and the pc components pSettings.pc. The value
%             pSettings.nPC tells you how many PC components were kept (this calls getInstancesPCA)
%             You can also keep a fixed number of pc's. for exmaple to keep 5:
%               [instances pSettings] = preprocessInstances(instances,'pca=5');
%
%             This function is usually called by buildClassifier and classifyInstance. The settings
%             used by buildClassifier are remembered and applied to classifyInstance. So for example:
%             In build the call could look like this:
%               [instances pSettings] = preprocessInstances(instances,'pca=5');
%             and in test:
%               [instances] = preprocessInstances(instances,pSettings);
%             Note that this will apply the PC xform found in the build to the test (i.e. it doesn't
%             compute a new PC basis for the test instances).
%             
%             
function [instances, pSettings] = preprocessInst(instances,pSettings)

% demean and zscore processing
if pSettings.demean || pSettings.zscore || pSettings.linearScale0_1
  % make data matrix which should be one row for each instance, one column for each voxel
  % so d = n x k where n is number of instances and k is number of voxels
  d = [];nClasses = length(instances);
  for iInstance = 1:nClasses
      if isempty(instances{iInstance})
          nInstances(iInstance) = 0;
          continue
      end
    d(end+1:end+size(instances{iInstance},1),:) = instances{iInstance};
    % remember how many instances there originaly was
    nInstances(iInstance) = size(instances{iInstance},1);
  end
  instances = {};

  % demean
  if pSettings.demean
    if ~isfield(pSettings,'demeanInfo')
      % remove mean across instances (i.e. each voxels response will be 0 across all instances)
      meand = mean(d);
      pSettings.demeanInfo.meand=meand;
    else
      meand=pSettings.demeanInfo.meand;
    end
    d = bsxfun(@minus,d,meand);
  end
    

  % z-score
  if pSettings.zscore
    % remove std across instances (i.e. each voxels response will have std=1 across all instances)
    if ~isfield(pSettings,'zscoreInfo')
      stdd = std(d);
      pSettings.zscoreInfo.stdd=stdd;
    else
      stdd=pSettings.zscoreInfo.stdd;
    end
    d = bsxfun(@rdivide,d,stdd);
  end
  
  % scale values so that all values are between 0 and 1
  if pSettings.linearScale0_1
      mind=min(d); maxd=max(d);
      d = bsxfun(@minus,d,mind);
      d = bsxfun(@rdivide,d,maxd-mind);
      pSettings.linearScaleInfo.mind = mind;
      pSettings.linearScaleInfo.maxd = maxd;
  end

  % now sort back into instances
  thisRow = 1;
  for iInstance = 1:nClasses
    thisEndRow = thisRow + nInstances(iInstance) - 1;
    instances{iInstance} = d(thisRow:thisEndRow,:);
    thisRow = thisEndRow+1;
  end
end