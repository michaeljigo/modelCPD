% Use to initialize and update best PEST or QUEST staircase as implemented
% in the Palamedes toolbox.

% Will output both the initialized staircase and the parameters used to
% create that staircase

% Sample stairParams are:

function [staircase, stairParams] = usePalamedesStaircase(stairParams,response)

%% Initialize staircase
% Initialize the parameters
if nargin==0
    % if nothing is passed in, ask experimenter for each of the parameters
    stairParams.whichStair = input(['Which staircase protocol do you want to use ',...
        '(1=best PEST; 2=QUEST; default=1)? ']);
    % check that whichStair is valid
    if ~ismember(stairParams.whichStair,[1 2])
        warning('Input can only be 1 or 2. Setting to default of 1.');
        stairParams.whichStair = [];
    end
    
    % if using the QUEST staircase, specify MEAN and SD of prior (gaussian
    % distribution)
    if stairParams.whichStair==2
        stairParams.questMean = input('Enter the mean of prior (default=median): ');
        stairParams.questSD = input('Enter the standard deviation of prior (default=1): ');
    end
    
    stairParams.alphaRange = input(['Enter the range, and spacing, of the ',...
        'stimulus values (default=0.01:0.01:1): ']);
    % check that alphaRange is valid
    if length(stairParams.alphaRange)==1
        warning('You must input a range of values. Setting to default.')
        stairParams.alphaRange = [];
    end
    
    stairParams.fitBeta = input(['Assumed slope (beta) of the ',...
        'underlying psychometric function (PF; default=2): ']);
    stairParams.fitLambda = input(['Assumped upper asymptote (lambda) of the ',...
        'underlying PF (default=0.01): ']);
    stairParams.fitGamma = input('Assumed lower asymptote (gamma) of the PF (default=0.5): ');
    
    stairParams.threshPerformance = input('Enter the expected performance level: ');
    
    % ask whether experimenter wants to use posterior probability from a
    % previous staircase run
    stairParams.lastPosterior = input(['Use a previous posterior probability? Input ',...
        'the probability distribution (default=[]): ']);
    
    % ask which function the experimenter wants to use
    stairParams.PF = input('Input the function you want to fit (default arbWeibull): ');
    
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
    
    %% UPDATE THE STAIRCASE
elseif nargin==2
    % When there are two inputs to this function, the first is assumed to
    % be the initialized the staircase and the second is assumed to be the
    % response (0=incorrect or 1=correct) to a trial.
    
    % Update the staircase
    if ismember(func2str(stairParams.PF),{'PAL_Weibull' 'PAL_Quick' 'PAL_Gumbel' ...
            'PAL_HyperbolicSecant' 'PAL_Logistic'})
        staircase = PAL_AMRF_updateRF(stairParams,stairParams.xCurrent,response);
    else
        staircase = updateArbWeibull(stairParams,stairParams.xCurrent,response,stairParams.threshPerformance);
    end
end

function staircase = initStaircase(stairParams)
% Initialize the staircase
switch stairParams.whichStair
    case 1 % best PEST
        meanMode = 'mode';
        prior = ones(1,length(stairParams.alphaRange));
        prior = prior/sum(prior); % make a uniform prior
    case 2 % QUEST
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

function stairParams = validateParams(stairParams)

fprintf('\n');
fprintf('\n');

% check if all parameters were set, if not, then set to the default
setParams = fieldnames(stairParams);
setParams = sort(setParams);

for i = 1:length(setParams)
    switch setParams{i}
        
        % setting possible threshold estimates
        case 'alphaRange'
            if isempty(stairParams.alphaRange)
                stairParams.alphaRange = 0.01:0.01:1;
                fprintf('ALPHA RANGE: Set to 0.01:0.01:1 (DEFAULT)\n');
            elseif any(stairParams.alphaRange<0)
                error('ALPHA RANGE: Cannot give negative values\n');
            end
            
            % check that inputtted PF is accurate
        case 'PF'
            if isempty(stairParams.PF)
                stairParams.PF = @arbWeibull;
                fprintf('PSYCHOMETRIC FUNCTION: Set to arbWeibull (DEFAULT)\n');
            else
                if ~isa(stairParams.PF,'function_handle')
                    stairParams.PF = eval(['@',stairParams.PF]);
                elseif ~ismember(func2str(stairParams.PF),{'arbWeibull' 'PAL_Weibull' 'PAL_Quick' 'PAL_Gumbel' ...
                        'PAL_HyperbolicSecant' 'PAL_Logistic'})
                    error('PSYCHOMETRIC FUNCTION: Inputted PF is not in the possible list of options.');
                end
            end
            
            % check for proper beta
        case 'fitBeta'
            if isempty(stairParams.fitBeta)
                stairParams.fitBeta = 2;
                fprintf('BETA: Set to 2 (DEFAULT)\n');
            end
            
            % check for proper gamma
        case 'fitGamma'
            if isempty(stairParams.fitGamma)
                stairParams.fitGamma = 0.5;
                fprintf('GAMMA: Set to 0.5 (DEFAULT)\n');
            end
            
            % check for proper lambda
        case 'fitLambda'
            if isempty(stairParams.fitLambda)
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
                stairParams.threshPerformance = [];
            end
            % if using arbWeibull, force user to input a threshold
            % performance
            if strcmp(func2str(stairParams.PF),'arbWeibull') && isempty(stairParams.threshPerformance)
                error(['EXPECTED PERFORMANCE: If using arbWeibull you '...
                    'need to input an expected performance level.']);
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
    end
end