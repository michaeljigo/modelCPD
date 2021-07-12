% Simulates a subject based on PF model.
%
% Usage: subj = simulateExp(pf,pfParams,levels,trials,nSim)
%
% INPUTS:
% pf = model PF; 
% can be:
% 'weibull'
% 'probit' or 'cumulative normal'
% 'log normal'
% 'logit' or 'logistic'
% 'quick'
%
% pfParams is a vector of parameters for each model PF; vectors should be
% arranged as follows: [alpha beta gamma lambda]
% if gamma or lambda are not given, they are assumed to be 0.5 and 0,
% respectively
%
% levels is a scalar correponding to the number of stimulus levels you want
% to simulate
%
% trials is a scalar representing the number of trials at each level
%
% nSim is a scalar representing how many simulated PFs will be produced
%
% OUTPUT:
% subj contains the following fields:
% modelPF = the model used to generate the data
% pfParams = the params used to generate the date
% levels = stimulus levels being simulated
% trials = # of trials at each level
% simPF = simulated PFs

function subj = simulateExp(pf,pfParams,levels,trials,nSim)

%% Create subject model
switch lower(pf)
    case 'weibull'
        levels = linspace(0,ceil(pfParams(1))+ceil(pfParams(1)),levels);
        
        % pfParams = [alpha beta gamma lambda]
        if length(pfParams)==2
            pfParams = [pfParams 0.5 0];
        end
        
        pfFun = pfParams(3) + ...
            (1-pfParams(3)-pfParams(4)).*(1-exp(-(levels./pfParams(1)).^pfParams(2)));
        
    case {'probit' 'cumulative normal' 'log normal'}
        if strcmpi(pf,'log normal')
            levels = log(linspace(1,exp(pfParams(1))+10,levels));
        else
            % levels = linspace(pfParams(1)-3,pfParams(1)+3,levels);
            levels = linspace(-1,1,levels);
        end
        
        % pfParams = [alpha beta gamma lambda]
        % for log normal, alpha must be in log units
        if length(pfParams)==2
            pfParams = [pfParams 0.5 0];
        end
        pfFun = pfParams(3)+(1-pfParams(3)-pfParams(4)).*normcdf(levels,pfParams(1),1/pfParams(2));
    case {'logit' 'logistic'}
        levels = log(linspace(1,exp(pfParams(1))+10,levels));
        % pfParams = [alpha beta gamma lambda]
        % alpha needs to be in log units
        if length(pfParams)==2
            pfParams = [pfParams 0.5 0];
        end
        pfFun = pfParams(3)+(1-pfParams(3)-pfParams(4)).*(1./(1+exp(-pfParams(2)*(levels-pfParams(1)))));
    case 'quick'
        levels = linspace(0,ceil(pfParams(1))+ceil(pfParams(1)),levels);
        % pfParams = [alpha beta gamma lambda]
        if length(pfParams)==2
            pfParams = [pfParams 0.5 0];
        end
        pfFun = pfParams(3)+(1-pfParams(3)-pfParams(4)).*(1-2.^(-(levels./pfParams(1)).^pfParams(2)));
end

%% Create trials for each simulation
for n = 1:nSim
    % trial responses at each level
    trialResp = rand(trials,length(levels));
    % subject response at each level (1=correct; 0=incorrect)
    subjResp = trialResp<repmat(pfFun,trials,1);
    % create simulated PF
    simPF(n,:) = mean(subjResp);
end

% Output structure
subj.modelPF = pf;
subj.pfParams = pfParams;
subj.levels = levels;
subj.trials = trials;
subj.simPF = simPF;
