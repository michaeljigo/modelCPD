% ----------------------------------------------------------------------
% [staircase, stairParams]=usePalamedesStaircase(stairParams,response)
% ----------------------------------------------------------------------
% Purpose:  Initialize and update adaptive procedures. 
%           This function simply calls the Palamedes toolbox, but makes
%           it easier to set and change parameters.
%
%           You can pass in a stairParams structure.
%           e.g., [staircase stairParams] = usePalamedesStaircase(stairParams);
%
%           If run with two inputs, the first input must be 'staircase'
%           (i.e., the output of initialization) and the second input must
%           be either 1 or 0 for correct and incorrect responses, respectively.
%           This will update the adaptive method.
%           e.g., staircase = usePalamedesStaircase(staircase,1);
% 
% ----------------------------------------------------------------------
% Input(s)
% stairParams           : structure containing parameters (listed below) to control adaptive method
%     whichStair        : 1=best PEST; 2=QUEST     (default=1)
%     alphaRange        : vector of possible threshold stimulus values  (default=0.01:0.01:1)
%     fitBeta           : slope of underlying psychometric function     (default=2)
%     fitLambda         : lapse rate (i.e., 1-upper asymptote) of psychometric function (default=0.01)
%     fitGamma          : guess rate (i.e., lower asymptote) of psychometric function   (default=0.5)
%     threshPerformance : target threhsold performance (must be specified if using arbWeibull, see PF)
%     lastPosterior     : posterior distribution from earlier run. when inputted, adaptive method will continue where previous run stopped (default=[])
%     PF                : shape of underlying psychometric function. can be: PAL_Weibull, PAL_Quick, PAL_Gumbel, PAL_HyperbolicSecant, PAL_Logistic, or arbWeibull(default)
%     updateAfterTrial  : if >0, the adaptive method will not use posterior-based estimates until the trial number matches the input for this variable (default=0)
%     preUpdateLevels   : these are the stimulus levels that will be tested before the adaptive method is updated (only works if updateAfterTrial>1)
% response              : 1=correct; 0=incorrect
%
% ----------------------------------------------------------------------
% Output(s)
% staircase             : structure controlling Palamedes adaptive method
% stairParams           : structure of parameters used to initialize Palamedes
% ----------------------------------------------------------------------
% Function created by Michael Jigo
% Last update : Feb. 17 2021
% ----------------------------------------------------------------------

function [staircase, stairParams] = usePalamedesStaircase(stairParams,response)

if nargin==0
   %% Use default staircase parameters
   paramnames = {'whichStair' 'alphaRange' 'fitBeta' 'fitLambda' 'fitGamma' 'threshPerformance' 'lastPosterior' 'PF' 'updateAfterTrial' 'preUpdateLevels'};
   for p = 1:numel(paramnames)
      stairParams.(paramnames{p}) = [];
   end

   % set defaults
   stairParams = validateParams(stairParams);

   % initialize the staircase
   staircase = initStaircase(stairParams);

elseif nargin==1
   % Use this if a structre of parameters are passed in
   % validate parameters
   stairParams = validateParams(stairParams);

   % initialize the staircase
   staircase = initStaircase(stairParams);

elseif nargin==2
   % When there are two inputs to this function, the first is assumed to
   % be the initialized staircase structure (NOT the stairParams structure) 
   % and the second is assumed to be the response (0=incorrect or 1=correct) to a trial.

   % Update the staircase
   if ismember(func2str(stairParams.PF),{'PAL_Weibull' 'PAL_Quick' 'PAL_Gumbel' ...
         'PAL_HyperbolicSecant' 'PAL_Logistic'})
      staircase = PAL_AMRF_updateRF(stairParams,stairParams.xCurrent,response);
   else
      staircase = update_arbitrary_pf(stairParams,stairParams.xCurrent,response,stairParams.threshPerformance);
   end
end


%%%% Helper functions %%%%   
function staircase = initStaircase(stairParams)
% Initialize the staircase
switch stairParams.whichStair
   case 1 % best PEST
      % uniform prior with the mode selected as xCurrent
      meanMode = 'mode';
      prior = ones(1,length(stairParams.alphaRange));
      prior = prior/sum(prior); % make a uniform prior
   case 2 % QUEST
      % normally-distributed prior with the mean selected as xCurrent
      meanMode = 'mean';
      prior = PAL_pdfNormal(stairParams.alphaRange,stairParams.questMean,stairParams.questSD);
end

% set prior to be last posterior if it is provided
if isfield(stairParams,'lastPosterior') && ~isempty(stairParams.lastPosterior)
   prior = stairParams.lastPosterior;
end

staircase = PAL_AMRF_setupRF('priorAlphaRange',stairParams.alphaRange,...
   'stopcriterion','trials','stoprule',inf,'beta',stairParams.fitBeta,...
   'lambda',stairParams.fitLambda,'gamma',stairParams.fitGamma,...
   'meanmode',meanMode,'PF',stairParams.PF,'prior',prior);
staircase.threshPerformance = stairParams.threshPerformance;
staircase.type = 'maxlikelihood';
staircase.updateAfterTrial = stairParams.updateAfterTrial;
staircase.preUpdateLevels = stairParams.preUpdateLevels;
if staircase.updateAfterTrial>0
   staircase.xCurrent = staircase.preUpdateLevels(1);
end


   
function stairParams = validateParams(stairParams)
% check if all parameters were set, if not, then set to the default
setParams = fieldnames(stairParams);
setParams = sort(setParams);

for i = 1:length(setParams)
   switch setParams{i}

      % setting possible threshold estimates
      case 'alphaRange'
         if ~isfield(stairParams,'alphaRange') || isempty(stairParams.alphaRange)
            stairParams.alphaRange = 0.01:0.01:1;
            fprintf('ALPHA RANGE: Set to 0.01:0.01:1 (DEFAULT)\n');
         end

         % check that inputtted PF is accurate
      case 'PF'
         if ~isfield(stairParams,'PF') || isempty(stairParams.PF)
            stairParams.PF = @arbWeibull;
            fprintf('PSYCHOMETRIC function: Set to arbWeibull (DEFAULT)\n');
         else
            if ~isa(stairParams.PF,'function_handle')
               stairParams.PF = eval(['@',stairParams.PF]);
            elseif ~ismember(func2str(stairParams.PF),{'arbWeibull' 'arbLogistic' 'PAL_Weibull' 'PAL_Quick' 'PAL_Gumbel' ...
                  'PAL_HyperbolicSecant' 'PAL_Logistic'})
            error('PSYCHOMETRIC function: Inputted PF is not in the possible list of options.');
            end
         end

         % check for proper beta
      case 'fitBeta'
         if ~isfield(stairParams,'fitBeta') || isempty(stairParams.fitBeta) 
            stairParams.fitBeta = 2;
            fprintf('BETA: Set to 2 (DEFAULT)\n');
         end

         % check for proper gamma
      case 'fitGamma'
         if ~isfield(stairParams,'fitGamma') || isempty(stairParams.fitGamma) 
            stairParams.fitGamma = 0.5;
            fprintf('GAMMA: Set to 0.5 (DEFAULT)\n');
         end

         % check for proper lambda
      case 'fitLambda'
         if ~isfield(stairParams,'fitLambda') || isempty(stairParams.fitLambda) 
            stairParams.fitLambda = 0.01;
            fprintf('LAMBDA: Set to 0.01 (DEFAULT)\n');
         end

         % check if experimenter wants to continue staircase from
         % previous run
      case 'lastPosterior'
         if isempty(stairParams.lastPosterior)
            stairParams.lastPosterior = [];
         elseif ~isfield(stairParams,'lastPosterior')
            stairParams.lastPosterior = [];
         else
            fprintf(['POSTERIOR: Staircase will continue from '...
               'previously computed posterior.\n']);
         end


         % check if experimenter wants to use an arbritrary threshold
         % performance
      case 'threshPerformance'
         if isempty(stairParams.threshPerformance)
            stairParams.threshPerformance = 0.75;
         end
         % if using arbWeibull, force user to input a threshold
         % performance
         if ~isempty(strfind(func2str(stairParams.PF),'arb')) && isempty(stairParams.threshPerformance)
            error(['EXPECTED PERFORMANCE: If using an arbitrary (modified) performance level you '...
               'need to input that expected performance level.']);
         end

         % choosing staircases
      case 'whichStair'
         if ismember(stairParams.whichStair,[1 2])
            if stairParams.whichStair==2
               % check for mean and sd input of prior
               questInput = {'questMean' 'questSD'};
               questDefault = {'stairParams.alphaRange(round(length(stairParams.alphaRange)/2))' '1'};
               idx = ismember(questInput,setParams);
               if any(idx)
                  % check if values are empty
                  questVal = {stairParams.(questInput{1}) ...
                     stairParams.(questInput{2})};
               emptyIdx = find(cellfun('isempty',questVal));
               if ~all(isnumeric([stairParams.questMean stairParams.questSD])) ...
                     && isempty(emptyIdx)
                  % check that values inputted are numbers
                  error('QUEST: mean and standard deviation should be numbers');
               else
                  % set empty values to the default
                  for q = 1:length(emptyIdx)
                     fprintf(['Setting ',questInput{emptyIdx(q)}, ' to ',...
                        questDefault{emptyIdx(q)}, ' (DEFAULT)\n']);
                  stairParams.(questInput{emptyIdx(q)}) = ...
                     eval(questDefault{emptyIdx(q)});
                  end
               end
               else
                  % set the empty or the non-set input to the
                  % default
                  fprintf(['Mean and SD of QUEST prior were not set. '...
                     'Setting them to default.']);
                  for q = 1:length(questInput)
                     stairParams.(questInput{q}) = eval(questDefault{q});
                  end
               end
            end
         else
            fprintf('whichStair: Setting to use best PEST (DEFAULT)\n');
            stairParams.whichStair = 1;
         end

      % trial number, after which, staircase will update
      case 'updateAfterTrial'
         if ~isfield(stairParams,'updateAfterTrial') || isempty(stairParams.updateAfterTrial) 
            stairParams.updateAfterTrial = 0;
         end

      % levels to display before the staircase gets updated
      case 'preUpdateLevels'
         if ~isfield(stairParams,'preUpdateLevels') || isempty(stairParams.preUpdateLevels) 
            if stairParams.updateAfterTrial>0
               stairParams.preUpdateLevels = repmat(median(stairParams.alphaRange),1,stairParams.updateAfterTrial);
            else
               stairParams.preUpdateLevels = [];
            end
         elseif numel(stairParams.preUpdateLevels)~=stairParams.updateAfterTrial
            warning('# of pre-update levels need to match number of pre-update trials. Using default (median of alpha range)');
            stairParams.preUpdateLevels = repmat(median(stairParams.alphaRange),1,stairParams.updateAfterTrial);
         end
   end
end
