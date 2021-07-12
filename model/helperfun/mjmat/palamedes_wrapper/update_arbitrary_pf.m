% ----------------------------------------------------------------------
% RF=update_arbitrary_pf(RF,amplitude,response,threshPerformance)
% ----------------------------------------------------------------------
% Purpose:  Update an arbitrary psychometric function (PF). 
%           The PF is arbitrary to the extent that its "base" depends on the 
%           inputted threshPerformance variable (see more info on 
%           threshPerformance below). Such a flexible "base" allows the 
%           adaptive procedure to target any performance level an experimenter
%           desires.
%
%           This function is adapted from PAL_AMRF_updateRF.
%           Thus, it updates the structure which contains the settings for
%           and results of a running fit adaptive method.
%
% ----------------------------------------------------------------------
% Input(s)
% RF                 : structure containing settings for running fit, 
%                      as initialized by PAL_AMRF_setupRF
% amplitude          : tested stimulus value (e.g., contrast value) on current trial
% response           : observer response (1=correct; 0=incorrect) to tested stimulus value
% threshPerformance  : desired convergence performance level 
%
% ----------------------------------------------------------------------
% Output(s)
% RF                 : updated structure containing settings for running fit
% ----------------------------------------------------------------------
% Function adapted by Michael Jigo from the Palamedes toolbox
% Last update: 02.19.21
% ----------------------------------------------------------------------

function RF = update_arbitrary_pf(RF, amplitude, response, threshPerformance)

trial = length(RF.response)+1;
RF.x(trial) = amplitude;
RF.response(trial) = response;

if trial == 1
    RF.xStaircase(trial) = RF.x(trial);
    if response == 1
        RF.direction = -1;        
    else
        RF.direction = 1;
    end
end

%% Update posterior distribution
% unpack parameter space for family of psychometric functions
params.alpha = RF.priorAlphaRange;
params.beta = RF.beta;
params.gamma = RF.gamma;
params.lambda = RF.lambda;

% evaluate probability of correct response at given stimulus amplitude for each psychometric function
p = RF.PF(params,amplitude,threshPerformance);

% update posterior probability distribution
if response == 1
    RF.pdf = RF.pdf.*p;
else
    RF.pdf = RF.pdf.*(1 - p);
end
RF.pdf = RF.pdf./sum(RF.pdf);

% obtain descriptive stats of posterior
[RF.mode, RF.mean, RF.sd] = PAL_AMRF_pdfDescriptives(RF.pdf, RF.priorAlphaRange);

% choose next contrast level
if strcmpi(RF.meanmode,'mean')
    if (RF.mean > RF.xCurrent && RF.direction == 1) || (RF.mean < RF.xCurrent && RF.direction == -1)
        RF.reversal(trial) = 0;
    end
    if RF.mean > RF.xCurrent && RF.direction == -1
        RF.reversal(trial) = sum(RF.reversal~=0) + 1;
        RF.direction = 1;
    end
    if RF.mean < RF.xCurrent && RF.direction == 1
        RF.reversal(trial) = sum(RF.reversal~=0) + 1;
        RF.direction = -1;
    end
    RF.xCurrent = RF.mean;
end
if strcmpi(RF.meanmode,'mode')
    if (RF.mode > RF.xCurrent && RF.direction == 1) || (RF.mode < RF.xCurrent && RF.direction == -1)
        RF.reversal(trial) = 0;
    end
    if RF.mode > RF.xCurrent && RF.direction == -1
        RF.reversal(trial) = sum(RF.reversal~=0) + 1;
        RF.direction = 1;
    end
    if RF.mode < RF.xCurrent && RF.direction == 1
        RF.reversal(trial) = sum(RF.reversal~=0) + 1;
        RF.direction = -1;
    end
    RF.xCurrent = RF.mode;
end

% add in functionality for delaying the update for a given number of trials
if isfield(RF,'updateAfterTrial')
   if numel(RF.x)<RF.updateAfterTrial
      % until the trial-update-threshold has been reached, continue presenting the custom contrast levels
      RF.xCurrent = RF.preUpdateLevels(numel(RF.x)+1);
   end
end

RF.xStaircase(trial+1) = RF.xCurrent;

if (strncmpi(RF.stopCriterion,'reversals',4) && sum(RF.reversal~=0) == RF.stopRule)||(strncmpi(RF.stopCriterion,'trials',4) && trial == RF.stopRule)
    RF.stop = 1;
    [RF.modeUniformPrior, RF.meanUniformPrior, RF.sdUniformPrior] = PAL_AMRF_pdfDescriptives(RF.pdf./RF.prior, RF.priorAlphaRange);
end
