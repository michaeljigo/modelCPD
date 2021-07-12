% Purpose:  Objective function for cross-validation of different functional shapes, which are all based on an asymmetric
%           generalized Gaussian.

function [p trainErr testErr] = performFit2AsymGaus(train,test,fixedParams,bounds,seedParams)
% fitting options
options = optimoptions(@fmincon,'Algorithm','active-set','Display','none','MaxIter',1e6,'functionTolerance',1e-6,'StepTolerance',1e-14,'MaxfunctionEvaluations',1e50);
ms = MultiStart('XTolerance',1e-6,'Display','off','StartPointsToRun','bounds','UseParallel',0);

% repeat bounds to account for all eccentricities
lb = repmat(bounds(1,:),size(seedParams,1),1);
ub = repmat(bounds(2,:),size(seedParams,1),1);

% create fitting anonymous function
fitFun = @(params)fit_asymGauss(train,params,fixedParams);
modelProb = createOptimProblem('fmincon','objective',fitFun,'x0',seedParams(:),'lb',lb(:),'ub',ub(:),'options',options);
[p trainErr] = run(ms,modelProb,50);

% reshape output parameters
p = reshape(p,size(train.attn,1),numel(fixedParams));

% if fixing any parameters, do fixing now
for e = 1:size(p,1)
   for fp = 1:numel(fixedParams)
      if isnan(fixedParams{fp})
         continue
      else
         if ischar(fixedParams{fp})
            p(e,fp) = eval(fixedParams{fp});
         else
            p(e,fp) = fixedParams{fp};
         end
      end
   end
end

% compute error of test set
for e = 1:size(test.oct,1)
   if numel(p(e,:))==7
      % two peak Gaussian
      expected(e,:) = asymGaussian_twoPeak(test.oct(e,:),p(e,:));
   else
      % one peak asym Gaussian
      expected(e,:) = asymGaussian(test.oct(e,:),p(e,:));
   end
end
% remove 0s from attention effects
%zeroIdx = ~isnan(test.attn) & test.attn==0;
%test.attn(zeroIdx) = nan;
%test.nsubj = double(test.nsubj);
%test.nsubj(zeroIdx) = nan;

testErr = sqrt(nanmean(test.nsubj(:).*((expected(:)-test.attn(:)).^2)));
%testErr = testErr./nanmedian(test.attn(:).*test.nsubj(:));



function [cost,model] = fit_asymGauss(data,p,fixedParams)
% parameters for the asymmetric gaussian are
% p(1) = mu (center)
% p(2) = sigma left (spread to the left of center)
% p(3) = sigma right (spread to the right of center)
% p(4) = exponent (shape)
% p(5) = amplitude

% reshape parameters to have e x 6 matrix
p = reshape(p,size(data.attn,1),numel(fixedParams));

% if fixing any parameters, do fixing now
for e = 1:size(p,1)
   for fp = 1:numel(fixedParams)
      if isnan(fixedParams{fp})
         continue
      else
         if ischar(fixedParams{fp})
            p(e,fp) = eval(fixedParams{fp});
         else
            p(e,fp) = fixedParams{fp};
         end
      end
   end
end

% evaluate model
for e = 1:size(data.attn,1)
   if numel(p(e,:))==7
      % two peak Gaussian
      model(e,:) = asymGaussian_twoPeak(data.oct(e,:),p(e,:));
   else
      % one peak asym Gaussian
      model(e,:) = asymGaussian(data.oct(e,:),p(e,:));
   end
end

% remove 0s from data
%zeroIdx = ~isnan(data.attn) & data.attn==0;
%data.nsubj = double(data.nsubj);
%data.nsubj(zeroIdx) = nan;
%data.attn(zeroIdx) = nan;
%model(zeroIdx) = nan;

%data.attn(data.attn<0) = 0;

% compute cost of function
cost = sqrt(nanmean(data.nsubj(:).*((model(:)-data.attn(:)).^2)));
