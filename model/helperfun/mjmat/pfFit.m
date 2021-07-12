% fit PF with inputted PF model
% data is a matrix where each row contains data that needs to be fit
%
% Fits a model to psychometric function with maximum likelihood estimation.
%
% Usage: fits =
% pfFit(fitModel,data,levels,trialsPerLevel,startParams,paramLB,paramUB)
%
% INPUT:
% fitModel can be:
% 'weibull'
% 'probit' or 'cumulative normal'
% 'log normal'
% 'logit' or 'logistic'
% 'quick'
%
% data = vector or matrix consisting of the probability correct at each
% stimulus level. Each row is a new PF and each column is a stimulus level.
%
% levels = the stimulus levels that the PF were evaluated at
%
% trialsPerLevel = scalar or vector of the # of trials per level. If input
% is a scalar, each level will be assumed to have the same number of trials
%
% startParams = initial guesses for fitting. for fixed parameters, the
% initial guess will be used as the final parameter.
%
% paramLB = lower bound of each parameter
%
% paramUB = upper bound of each parameter
%
% if any element in paramLB==paramUB, that parameter will be fixed at that value
%
% OUTPUT:
% fits is a structure with the following fields:
% params = nxlevels matrix where n is equal to the number of rows in data.
% This contains the best-fitting parameters for each PF in each row of
% data.
%
% model = the name of the model used to fit data

function fits = pfFit(fitModel,data,levels,trialsPerLevel,startParams,paramLB,paramUB)

%% initialize parameters
% if strcmp(fitModel,'custom')
%     if ieNotDefined('customFun')
%         error('Input the function that will be used to fit.');
%     end
% else
%     customFun = [];
% end

options = optimset('Display','off','MaxIter',1e5,'TolX',1e-40,'TolFun',1e-40);

if strcmp(fitModel,'exp')
    fits.params = nan(size(data,1),2);
else
    fits.params = nan(size(data,1),length(startParams));
end

%% convert proportion correct to number of trials
if length(trialsPerLevel)==1
    trialsPerLevel = repmat(trialsPerLevel,1,size(data,2));
    nCorr = data.*repmat(trialsPerLevel,size(data,1),1);
elseif length(trialsPerLevel)==size(data,2)
    nCorr = data.*repmat(trialsPerLevel,size(data,1),1);
else
    error(['trialsPerLevel should be a scalar or a vector whose length is ',...
        'equal to the number of columns in data']);
end

%% fit using maximum likelihood estimation

% REMOVE LATER
% disppercent(-inf,'Fitting...');
for i = 1:size(data,1)
    switch lower(fitModel)
        case 'weibull'
            fun = @(p)weibullFit(levels,nCorr(i,:),trialsPerLevel,p);
        case {'probit' 'cumulative normal' 'log normal'}
            fun = @(p)cumNormalFit(levels,nCorr(i,:),trialsPerLevel,p);
        case {'logit' 'logistic'}
            fun = @(p)logitFit(levels,nCorr(i,:),trialsPerLevel,p);
        case {'quick'}
            fun = @(p)quickFit(levels,nCorr(i,:),trialsPerLevel,p);
        case {'gumbel'}
            fun = @(p)gumbelFit(levels,nCorr(i,:),trialsPerLevel,p);
        case {'exp'}
            fun = @(p)exponentialFit(levels,data(i,:),p);
        case {'adaptattnmodel'}
            %             fun = @(p)eval([customFun,'(levels,nCorr(i,:),trialsPerLevel,p);']);
            fun = @(p)adaptAttnModel(levels,nCorr(i,:),trialsPerLevel,p);
    end
    fits.params(i,:) = fminsearchbnd(fun,startParams,paramLB,paramUB,options);
    
    % compute and store likelihood
    F = evalPF(fitModel,levels,fits.params(i,:));
    fits.logLikelihood(i,:) = computeLL(F,nCorr(i,:),trialsPerLevel);
    fits.fitPF(i,:) = F;
    
    % REMOVE LATER
%     disppercent(i/size(data,1));
end
% REMOVE LATER
% disppercent(inf,'Fit.');
fits.model = fitModel;

function F = exponentialFit(x,y,params)
if any(isinf(y))
    x(isinf(y)) = [];
    y(isinf(y)) = [];
end

% least squares
F = sum((y-(params(1).*exp(params(2).*x))).^2);

function F = gumbelFit(x,nCorr,trialsPerLevel,params)

F = params(3)+(1-params(3)-params(4)).*(1-exp(-10.^(params(2).*(x-params(1)))));
F = -computeLL(F,nCorr,trialsPerLevel);

function F = weibullFit(x,nCorr,trialsPerLevel,params)

F = params(3)+(1-params(3)-params(4)).*(1-exp(-(x./params(1)).^params(2)));
F = -computeLL(F,nCorr,trialsPerLevel);

function F = cumNormalFit(x,nCorr,trialsPerLevel,params)

F = params(3)+(1-params(3)-params(4)).*normcdf(x,params(1),1/params(2));
F = -computeLL(F,nCorr,trialsPerLevel);

function F = logitFit(x,nCorr,trialsPerLevel,params)

F = params(3)+(1-params(3)-params(4)).*(1./(1+exp(-params(2)*(x-params(1)))));
F = -computeLL(F,nCorr,trialsPerLevel);

function F = quickFit(x,nCorr,trialsPerLevel,params)

F = params(3)+(1-params(3)-params(4)).*(1-2.^(-(x./params(1)).^params(2)));
F = -computeLL(F,nCorr,trialsPerLevel);
