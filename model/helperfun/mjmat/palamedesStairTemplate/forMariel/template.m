% A short template of using the staircase

%% Initialize
clear

stairParams.whichStair = 1;
stairParams.alphaRange = 0.5:0.5:45;
stairParams.fitBeta = 2;
stairParams.fitLambda = 0.01;
stairParams.fitGamma = 0.5;
stairParams.threshPerformance = 0.75;
stairParams.PF = 'arbWeibull';

myStair = usePalamedesStaircase(stairParams);

% I will be simulating an ideal observer to show that the staircase is
% working:
subjPF = @PAL_Weibull;
subjParams = [15 2 0.5 0.15];

%% Run 1st block
nTrials = 30;
for trial = 1:nTrials
    % show the threshold estimate to the subject
    targetTilt = myStair.xCurrent;
    
    % now update the staircase based on accuracy
    accuracy = rand(1)<subjPF(subjParams,targetTilt);
    myStair = usePalamedesStaircase(myStair,accuracy);  
end

% Plotting the staircase will show that it is working.
figure;
subplot(1,2,1)
scatter(1:nTrials,myStair.x,50,'k','filled'); hold on
plot(1:nTrials,myStair.x,'k-');
xlabel('trial #'); ylabel('tilt');
title('First block');
set(gca,'YLim',[0 45]);

%% Carryover staircase between blocks

stairParams.lastPosterior = myStair.pdf;
myStair = usePalamedesStaircase(stairParams);

%% No carryover
% myStair = usePalamedesStaircase(stairParams);

%% Run 2nd block
for trial = 1:nTrials
    % show the threshold estimate to the subject
    targetTilt = myStair.xCurrent;
    
    % now update the staircase based on accuracy
    accuracy = rand(1)<subjPF(subjParams,targetTilt);
    myStair = usePalamedesStaircase(myStair,accuracy);  
end

subplot(1,2,2)
scatter(1:nTrials,myStair.x,50,'k','filled'); hold on
plot(1:nTrials,myStair.x,'k-');
xlabel('trial #'); ylabel('tilt');
title('Second block');
set(gca,'YLim',[0 45]);

% For a test, try commenting out the carryover and see how the staircase
% will act differently on the 2nd block.
