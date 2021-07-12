% usage:    modelParams = fitCSF(model,sfs,cs,<lowerBound>,<upperBound>,<dispFig>)
% by:       Michael Jigo
% date:     07/08/18
% purpose:  Fit a contrast senstivitiy function (CSF) to contrast sensitivity data.
%
% INPUTS:
% model        string specifying which model to fit
%  'emg'       exponential minus gaussian
%  'hmg'       hyperbolic secant minus gaussian
%  'hpmg'      hyperbolic secant raised to power minus gaussian
%  'ms'        Mannos & Sakrison, 1974
%  'yqm'       Yang, Qi, & Makous, 1995
%  'dexp'      double exponential (following Rohaly & Owsley, 1993; Equation 3)
%  'apf'       asymmetric parabolic function (Chung & Legge, 2016)
%
%  sfs         vector of tested spatial frequencies
%
%  cs          vector of contrast sensitivtiy (1/threshold) at each spatial frequency
%
% lowerBound   scalar or vector specifying the lower bound of each model parameter  
%
% upperBound   scalar or vector specifying the upper bound of each model parameter
%
% dispFig      display the fitted model or not (default=0)

function modelParams = fitCSF(model,sfs,cs,lowerBound,upperBound,dispFig)

%% Defaults
% ensure that cs is a row vector
if ~isvector(cs)
   error('Contrast sensitivities must be in a row vector');
elseif size(cs,1)>1
   cs = cs';
end
if ~isvector(sfs)
   error('Spatial frequencies must be in a row vector');
elseif size(sfs,1)>1
   sfs = sfs';
end

if ~exist('dispFig','var')
   dispFig = 0;
end

% store original spatial frequencies, in case transformations are required
orgSF = sfs;

%% Select model and create objective function
switch lower(model)
   % NOTES: 
   % f = spatial frequencies
   % f0 = high-SF scalar
   % f1 = low-SF scalar
   % a = scalar on high or low SF (depending on form)
   % g = uniform gain (common to all forms)
   case 'emg'
      % exponential minus gaussian
      paramNames = {'f0' 'f1' 'a' 'g'};
      eqn = @(p,sfs) p(4)*(exp(-sfs/p(1))-p(3)*exp(-(sfs/p(2)).^2)); 
   case 'hmg'
      % hyperbolic secant minus gaussian
      paramNames = {'f0' 'f1' 'a' 'g'};
      eqn = @(p,sfs) p(4)*(sech(sfs/p(1))-p(3)*exp(-(sfs/p(2)).^2));
   case 'hpmg'
      % hyperbolic secant, raised to power, minus gaussian
      paramNames = {'f0' 'f1' 'a' 'p' 'g'};
      eqn = @(p,sfs) p(5)*(sech((sfs/p(1)).^p(4))-p(3)*exp(-(sfs/p(2)).^2));
   case 'hmh'
      % hyperbolic secant minus hyperbolic secant
      paramNames = {'f0' 'f1' 'a' 'g'};
      eqn = @(p,sfs) p(4)*(sech(sfs/p(1))-p(3)*sech(sfs/p(2)));
   case 'hpmh'
      % hyperbolic secant, raised to power, minus hyperbolic secant
      paramNames = {'f0' 'f1' 'a' 'p' 'g'};
      eqn = @(p,sfs) p(5)*(sech((sfs/p(1)).^p(4))-p(3)*sech(sfs/p(2)));
   case 'ms'
      % Mannos & Sakrison, 1974 
      paramNames = {'f0' 'a' 'p' 'g'};
      eqn = @(p,sfs) 10.^(p(4)*((1-p(2)+sfs/p(1)).*exp(-(sfs/p(1)).^p(3))));
   case 'yqm'
      % Yang, Qi, & Makous, 1995
      paramNames = {'f0' 'f1' 'a' 'g'};
      eqn = @(p,sfs) p(4)*(exp(-sfs/p(1))./(1+(p(3)./(1+(sfs/p(2)).^2))));
   case 'dexp'
      % double exponential (following Rohaly & Owsley, 1993; Equation 3)
      paramNames = {'f0' 'p' 'g'};
      eqn = @(p,sfs) p(3)*sfs.^p(2).*exp(-p(1)*sfs);
   case 'apf'
      % asymmetric parabolic function (Chung & Legge, 2016)
      paramNames = {'a' 'fMax' 'wLow' 'wHigh'};
      sfs = log10(sfs);
      eqn = @(p,sfs) 10.^[p(1)-(sfs(sfs<p(2))-p(2)).^2.*p(3)^2, ...
         p(1)-(sfs(sfs>=p(2))-p(2)).^2.*p(4)^2];
end
initP = rand(1,length(paramNames));

% bounds
if ~exist('upperBound') || isempty(upperBound)
   upperBound = inf(1,length(initP));
elseif isscalar(upperBound)
   upperBound = repmat(upperBound,1,length(initP));
elseif length(upperBound)~=length(initP)
   error('Length of bound must be %i',length(initP));
end
if ~exist('lowerBound') || isempty(lowerBound)
   lowerBound = [-inf(1,length(initP))];
elseif isscalar(lowerBound)
   lowerBound = [repmat(lowerBound,1,length(initP))] ;
elseif length(lowerBound)~=length(initP)
   error('Length of bound must be %i',length(initP));
end

%% Perform fitting
fitType = 'nonlinear';
switch fitType
   case 'linear'
      % objective function (least-square regression)
      options = optimset('Display','off','MaxFunEvals',1e5,'TolFun',1e-14);
      fun = @(p) sum((cs-eqn(p,sfs)).^2);
      [p, fval] = fminsearchbnd(fun,initP,lowerBound,upperBound,options);
   case 'nonlinear'
      options = optimoptions('lsqcurvefit','Display','off',...
         'MaxFunctionEvaluations',1e5,'FunctionTolerance',1e-5);
      [p, fval] = lsqcurvefit(eqn,initP,sfs,cs,lowerBound,upperBound,options);
end

%% Create output structure
modelParams.model =  model;
modelParams.paramNames = paramNames;
modelParams.params = p;
modelParams.minErr = fval;
modelParams.fitType = fitType;

% add in high-SF cutoff
switch model
   case 'apf'
      highSFCutoff = 10.^(sqrt(p(1)./(p(4)^2))+p(2));
      modelParams.highSFCutoff = highSFCutoff;
end


%% Plot data and fit
if dispFig
   fitX = linspace(min(orgSF)/2,max(orgSF)*2,5e2);
   fitY = evalCSF(modelParams,fitX);
   loglog(orgSF,cs,'ks'); hold on
   loglog(fitX,fitY,'k-');
   xlabel('Spatial frequency');
   ylabel('Sensitivity');
   set(gca,'yLim',[1 150]);
end
