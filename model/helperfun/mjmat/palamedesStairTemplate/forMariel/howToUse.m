%% Initialize
% The parameters of the staircase have to be set. This can be done
% in two ways:

%%%%% 1) Manually create a stairParams variable (PREFERRED METHOD)

% stairParams is a structure with fields holding the parameter values. For
% instance, a best PEST staircase controlling the tilt of a gabor in a 2AFC
% could be set up as follows:

stairParams.whichStair = 1; % 1=best PEST; 2=QUEST
stairParams.alphaRange = 0.5:0.5:45; % these are the possible Gabor tilts in degrees(*)
stairParams.fitBeta = 2; % slope parameter; a slope of 2 works well
stairParams.fitLambda = 0.01; % lapse parameter
stairParams.fitGamma = 0.5; % lower bound/chance parameter
stairParams.threshPerformance = 0.75; % expected performance level @ threshold
stairParams.PF = 'arbWeibull'; % (**)

% Then, initialize the staircase (myStair) with these parameters by invoking
% usePalamedesStaircase
myStair = usePalamedesStaircase(stairParams);

% myStair is a strucutre with various fields. We'll only focus on the
% response and xCurrent fields for now. Response contains the subject's
% accuracy (should be empty) and xCurrent contains the current threshold
% estimate (by default, this is the median of alphaRange).

% (*) Note that when using best PEST, the inputted possible Gabor tilts
% will be the only tilts that can be shown to the subject. If tilts are
% coarsely sampled (e.g., 0.5:15:45), there will only be a few threshold
% estimates to choose from and the correct threshold will be almost 
% impossible to estimate. Therefore, a general rule is that the range of
% alphaRange should closely match the expected range of stimulus values 
% (e.g., tilts) in the subject's psychometric function.

% (**) The advantage of best PEST and QUEST is that you can arbitrarily set
% the level of performance. The arbWeibull function is a custom PF that
% will ensure that the estimated threshold will converge at a performance
% level within ~5% of threshPerformance.

%%%%% 2) Follow the prompts specified by usePalamedesStaircase

% You can also run usePalamedesStaircase without any input, and prompts will
% be printed in the terminal that allow you to specify the parameters. The
% ouptut will be the stairParams variable. For example,

myStair = usePalamedesStaircase;

%% Update

% When the staircase needs to be updated, simply use the subject's accuracy
% (0 or 1, for incorrect or correct respectively) as the second input to
% usePalamedesStaircase. For example, if the subject is correct do:

myStair = usePalamedesStaircase(myStair,1);

% When run like this, the function assumes that the subject was presented
% with the threshold estimate contained in myStair.xCurrent. After updating
% the staircase, you will notice that xCurrent has changed to the newest
% threshold estimate. Additionally, response will contain values
% corresponding to the accuracy of the subject.

%% Carryover

% If you want to resume the staircase from a previous block, only the
% posterior of the staircase needs to be carried over. To this:

% 1) Take the posterior from the previous block and add it as a parameter
% in stairParams

stairParams.lastPosterior = myStair.pdf; %(***)

%(***) This assumes that stairParams was defined previously and that it is
%still an active variable in your workspace. If not, you have to re-define
%stairParams as shown in "Initialize"

% 2) Re-initialize the staircase to start from where the previous block
% ended

myStair = usePalamedesStaircase(stairParams);

% You will now see the message: "Staircase will continue from previously
% computed posterior."