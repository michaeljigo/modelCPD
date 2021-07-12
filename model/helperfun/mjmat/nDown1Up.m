% usage:    params = nDown1Up(params,correct)
% by:       Michael Jigo
% date:     05/04/18
% purpose:  Initialize and update an n-down 1-up staircase protocol*.
%
%           *Step size cannot be adaptively changed (e.g., PEST) in current form.
%
% INPUTS:
% params    structure containing the following fields:
%     n                 # of correct responses needed to step down (default (d) = 1)
%     stepDown          units covered in a single downward step (d = 0.23)
%     stepUp            units covered in a single upward step (d = 0.7)
%     startLevel        stimulus level where the staircase will start (d = 0)
%     startAtReversal   the n-down staircase procedure will start at the 1st reversal
%     accuracy          binary vector of responses
%     reversal          binary vector of reversal trials
%     allLevels         vector with all presented stimulus levels
%     currentLevel      stimulus level that will be displayed next
%     lastDirection     scalar stating whether the last step was down (-1) or up (1)
%     stepTilReversal   step size that will be used until the 1st reversal if
%                       startAtReversal is used
%     nBoundTrials      scalar specifying the "catch" or "lapse" trials in the
%                       staircase
%     boundLevel        scalar specifying the stimulus level that will be presented
%                       on each "catch" or "lapse" trial
%     stayWithinBounds  specifies whether the staircase will present values within
%                       some inputted range. For example, an input of "nan" will 
%                       make the range unrestricted. A vector input (e.g., [0 1]) will
%                       prevent the staircase from showing stimulus levels less than 0
%                       or greater than 1.
%
% correct   scalar specifying if the observer was correct (1) or not (0)
%
% OUTPUTS:
% params    same as above

function params = nDown1Up(params,correct)

%% Initialize/validate params
validParams = {'n' 'stepDown' 'stepUp' 'startLevel' 'startAtReversal' 'accuracy' ...
   'reversal' 'allLevels' 'currentLevel' 'lastDirection' 'nTrials' ...
   'stepTilReversal' 'nBoundTrials' 'boundLevel' 'stayWithinBounds'};
defaultParams = {1 0.23 0.7 0 0 [] [] [] 0 nan 50 0.1 0 nan nan};
switch nargin
   case 0
      % use default parameters
      for p = 1:length(validParams)
         params.(validParams{p}) = defaultParams{p};
      end
      params.currentLevel = params.startLevel;
      params.boundTrials = nan(1,params.nTrials);
   case 1
      % get the inputted parameters and values
      inParams = fieldnames(params);
      for i = 1:length(inParams)
         inVal{i} = params.(inParams{i});
      end

      % if any parameters were undefined, use the defaults
      missingParams = ~ismember(validParams,inParams);
      allParams = [inParams; validParams(missingParams)'];
      allVal = [inVal'; defaultParams(missingParams)'];
      for i = 1:length(allParams)
         params.(allParams{i}) = allVal{i};
      end
      params.currentLevel = params.startLevel;

      % insert the catch trial vector
      params.boundTrials = nan(1,params.nTrials);
      if params.nBoundTrials>0
         params.boundTrials(randsample(params.nTrials,params.nBoundTrials)) = ...
            params.boundLevel;
      end
end
params.type = 'downup';

%% Update
steps = [params.stepDown params.stepUp];
updateStair = 0;
reversal = 0;
direction = params.lastDirection;
switch nargin
   case 2
      % accuracy
      params.accuracy(end+1) = correct;
      params.allLevels(end+1) = params.currentLevel;
      if correct, c = [-1 0]; else, c = [0 1]; end

      % trial index
      trialIdx = length(params.accuracy);

      % remove boundary trials, before assessing whether to update staircase
      acc = params.accuracy(isnan(params.boundTrials(1:trialIdx)));
      level = params.allLevels(isnan(params.boundTrials(1:trialIdx)));

      % check if a step should occur and its direction
      if  isnan(params.boundTrials(trialIdx)) && ~acc(end)
         updateStair = 1;
         direction = 1; % upward
      else
         % check if correctness satisfies step rule
         if length(acc)>=params.n && isnan(params.boundTrials(trialIdx))
            acc = acc(end-(params.n-1):end);
            level = level(end-(params.n-1):end);
            level = round(level*1e3)/1e3;
            if sum(acc)==params.n && length(unique(level))==1
               updateStair = 1;
               % expected direction of step
               direction = -1; % downard
            end
         elseif ~isnan(params.boundTrials(trialIdx))
            params.allLevels(end) = params.boundTrials(trialIdx);
      end
   end

   % reversals
   if ~isnan(params.lastDirection) && direction~=params.lastDirection
      reversal = max(params.reversal)+1;
   end
   params.reversal(end+1) = reversal;

   % determine update based on reversal rule
   if params.startAtReversal && ~any(params.reversal==1) && ...
         isnan(params.boundTrials(trialIdx))
      % staircase procedure (i.e., contingency on consecutive correct responses)
      % will not start until the first reversal occurs
      updateStair = 1;
      direction = c(c~=0);
      steps = repmat(params.stepTilReversal,1,2); % use arb-defined step size
   end
   % update last step direction
   params.lastDirection = direction;

   % now update stair based on accuracy
   if updateStair
      nextLevel = sum(c.*steps);
      params.currentLevel = params.currentLevel+nextLevel;
   else
      params.currentLevel = params.currentLevel;
   end

   % determine whether the staircase should stay within the bounds
   if all(~isnan(params.stayWithinBounds))
      if params.currentLevel<min(params.stayWithinBounds)
         params.currentLevel = min(params.stayWithinBounds);
      elseif params.currentLevel>max(params.stayWithinBounds)
         params.currentLevel = max(params.stayWithinBounds);
      end
   end
end
