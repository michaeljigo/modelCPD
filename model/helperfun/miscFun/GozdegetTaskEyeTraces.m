% GozdegetTaskEyeTraces.m
%
%        $Id:$ 
%      usage: e = GozdegetTaskEyeTraces(stimfileName,<taskNum=1>,<phaseNum=1>,<dispFig=0>,<dataPad=3>,<removeBlink=1>)
%         by: justin gardner
%         partially by: Gözde ?entürk
%       date: 13/10/14
%    purpose: same as getTaskEyeTraces
%             Additionally,
%             This code can retrieve the exact timing of certain eye
%             events: saccades, blinks, and fixations. There are some cases
%             that the timing of the experiment and eyelink timer do not
%             exactly overlap. To prevent any minimal timing error, the eye
%             event times are returned back in every milisecond. A note to
%             user: I am a novice in programming. Please do not hesitate to
%             correct my mistakes and improve the code.
%             
%    e = getTaskEyeTraces('100616_stim01','taskNum=1','phaseNum=1','dataPad=3');



function e = GozdegetTaskEyeTraces(stimfileName,varargin)

e = [];

% check arguments
if nargin == 0
  help getTaskEyeTraces
  return
end

taskNum=[];phaseNum=[];dispFig=[];dataPad=[];removeBlink=[];
if exist('getArgs') == 2
  getArgs(varargin,{'taskNum=1','phaseNum=1','dispFig=0','dataPad=3','removeBlink=0','segNum=1'});
else
  disp(sprintf('(getTaskEyeTraces) To run this program you need functions from the mrTools distribution. \nSee here: http://gru.brain.riken.jp/doku.php/mgl/gettingStarted#initial_setup'));
  return
end

% if we are passed in the name of a file then assume it is a stimfile and load it through
% getTaskParameters
if isstr(stimfileName)
  [e stimfile] = getTaskParameters(stimfileName);
  keyboard
  if isempty(e)
    return
  end
  % make sure e/task are the correct form of cell array
  e = cellArray(e);
  stimfile.task = cellArray(stimfile.task,2);
else
  disp(sprintf('(getTaskEyeTraces) Must pass in stimfile name'));
  return
end

% get the correct task and phase
if taskNum > length(e)
  disp(sprintf('(getTaskEyeTraces) taskNum=%i out of range for this stimfile (1:%i)',taskNum,length(e)));
  return
end
if phaseNum > length(e{taskNum})
  disp(sprintf('(getTaskEyeTraces) phaseNum=%i out of range for this task (1:%i)',phaseNum,length(e{taskNum})));
  return
end

e = e{taskNum}(phaseNum);

% keep stimfile
e.stimfile = stimfile;
e.stimfile.taskNum = taskNum;
e.stimfile.phaseNum = phaseNum;

% check for taskID
if ~isfield(stimfile.task{taskNum}{phaseNum},'taskID')
  if ~isfield(stimfile.task{taskNum}{phaseNum},'collectEyeData')
    disp(sprintf('(getTaskEyeTraces) **** No taskID field found in task. This stimfile was probably generated with an older version of mgl/task. You need to update your mgl code. ****'));
    return
  else
    % for mglEyelink V1 messages only one task could collect data
    if (stimfile.task{taskNum}{phaseNum}.collectEyeData == 1)
      taskID = 0;
      phaseNum = 1; % start with phase 1 in the recorded mgl messages 
    else
      taskID = NaN;
    end
  end
else
  taskID = stimfile.task{taskNum}{phaseNum}.taskID;
end

% check eye tracker type
eyeTrackerType = stimfile.myscreen.eyeTrackerType;
if ~isequal(eyeTrackerType,'Eyelink')
  disp(sprintf('(getTaskEyeTraces) Loading of eye tracker type %s not implemented',eyeTrackerType));
  return
else
  eyeTrackerFileExt = 'edf';
end

% get the filename
if isfield(stimfile.myscreen.eyetracker,'datafilename')
  eyeTrackerFilename = stimfile.myscreen.eyetracker.datafilename;
else
  disp(sprintf('(getTaskEyeTraces) No associated eyetracker datafilename found in myscreen'));
  return
end

% check for file, should be in myscreen directory
eyeTrackerFilename = fullfile(stimfile.stimfilePath,sprintf('%s.%s',eyeTrackerFilename,eyeTrackerFileExt));
if ~isfile(eyeTrackerFilename)
  disp(sprintf('(getTaskEyeTraces) Could not find eye tracker file %s',eyeTrackerFilename));
  return
end

% replace tilde
if exist('mlrReplaceTilde') == 2
  eyeTrackerFilename = mlrReplaceTilde(eyeTrackerFilename);
else
  if ~isempty(findstr('~',eyeTrackerFilename))
    disp(sprintf('(getTaskEyeTraces) The ~ in filename %s may not be parsed correctly',eyeTrackerFilename));
  end
end

% load the file
disppercent(-inf,sprintf('(getTaskEyeTraces) Opening edf file %s',eyeTrackerFilename));
sprintf('\n');
edf = mglEyelinkEDFRead(eyeTrackerFilename,0);
disppercent(inf);
if isempty(edf),return,end

% get all the messages that match our taskID
% get the number of trials
edf.nTrials = max(edf.mgl.trialNum(find(edf.mgl.taskID == taskID)));

if removeBlink
    % blink window, extra padding (50ms)
    if removeBlink==1 % == true because we can't pass logical and test with getArgs, use 1-eps for 1s
        blinkWindow = ceil(0.005*edf.samplerate); % default
    elseif removeBlink > 1 % assume in samples
        blinkWindow = ceil(removeBlink);
    elseif removeBlink < 1 % assume seconds
        blinkWindow = ceil(removeBlink*edf.samplerate);        
    end
    for Bn = 1:numel(edf.blinks.startTime)
        blinks = (edf.gaze.time >= edf.blinks.startTime(Bn)-blinkWindow & ...
                         edf.gaze.time <= edf.blinks.endTime(Bn)+blinkWindow);
        edf.gaze.x(blinks) = NaN;
        edf.gaze.y(blinks) = NaN;
        edf.gaze.pupil(blinks) = NaN;
    end
end

% now process each trial
disppercent(-inf,sprintf('(getTaskEyeTraces) Extracting trial by trial data for %i trials',edf.nTrials));
% get the start times

for i = 1:max(edf.mgl.trialNum(edf.mgl.taskID == taskID))
  % get start time
  % find the segment 0 message
  %%% I think this should be segment 1--at least in my code segment 1 == seg1
  %%% and that seems to be what updateTask writes out. The seg==0 often includes
  %%% deadtime related to waiting for backtics, user start, etc
  segmentOneTime = edf.mgl.time((edf.mgl.taskID == taskID) &  ...
                                   (edf.mgl.trialNum==i) &  ...
                                   (edf.mgl.segmentNum==1));
  segNumTime = edf.mgl.time((edf.mgl.taskID == taskID) &  ...
			    (edf.mgl.trialNum==i) &  ...
			    (edf.mgl.segmentNum==segNum));
  
  % call this the startTime
  if ~isempty(segNumTime)
    startTime(i) = segNumTime;
  else
    startTime(i) = nan;
  end
  % end time is the start of the next trial
  if i > 1
    if ~isempty(segmentOneTime)
      endTime(i-1) = segmentOneTime;
    else
      endTime(i-1) = nan;
    end
  end
end

% make the end time of the last trial, at most as long as the longest trial so far
if ~isempty(endTime)
  maxTrialLen = max(endTime-startTime(1:end-1))+1;
  endTime(end+1) = min(max(edf.gaze.time),startTime(end)+maxTrialLen-1);
else
  % if we have only one trial, then use till the end of the data
  endTime(end+1) = max(edf.gaze.time);
  maxTrialLen = endTime-startTime+1;
end
maxTrialLen = maxTrialLen+dataPad;

% now get time between samples (isn't the sample rate availible?)
% timeBetweenSamples = median(diff(edf.gaze.time));
% get the time between samples in milliseconds
timeBetweenSamples = (1/edf.samplerate)*1000;

% figure out how large to make data array
e.eye.xPos = nan(edf.nTrials,ceil(maxTrialLen/timeBetweenSamples));
e.eye.yPos = nan(edf.nTrials,ceil(maxTrialLen/timeBetweenSamples));
e.eye.pupil = nan(edf.nTrials,ceil(maxTrialLen/timeBetweenSamples));

% put in time in seconds
e.eye.time = (0:(size(e.eye.xPos,2)-1))/edf.samplerate;

% go through each trial and populate traces
warning('off','MATLAB:interp1:NaNinY');
for iTrial = 1:edf.nTrials
  disppercent(iTrial/edf.nTrials);
  % the times for this trial
  thisTrialTimes = startTime(iTrial):timeBetweenSamples:endTime(iTrial);
  % get data for xPos, yPos and pupil traces form edf data 
  e.eye.xPos(iTrial,1:length(thisTrialTimes)) = interp1(edf.gaze.time,edf.gaze.x,thisTrialTimes,'linear',nan);
  e.eye.yPos(iTrial,1:length(thisTrialTimes)) = interp1(edf.gaze.time,edf.gaze.y,thisTrialTimes,'linear',nan);
  e.eye.pupil(iTrial,1:length(thisTrialTimes)) = interp1(edf.gaze.time,edf.gaze.pupil,thisTrialTimes,'linear',nan);
end
warning('on','MATLAB:interp1:NaNinY');
disppercent(inf);

% convert to device coordinates
w = stimfile.myscreen.screenWidth;
h = stimfile.myscreen.screenHeight;
xPix2Deg = stimfile.myscreen.imageWidth/w;
yPix2Deg = stimfile.myscreen.imageHeight/h;

e.eye.xPos = (e.eye.xPos-(w/2))*xPix2Deg;
e.eye.yPos = ((h/2)-e.eye.yPos)*yPix2Deg;
%% I dot want to extrapolate. Thus, I created this new gazetime

keyboard
gaze_time=edf.gaze.time(1):edf.gaze.time(end);
keyboard

% Get fixations
fixstart=ismember(gaze_time,edf.fixations.startTime); 
fixend=ismember(gaze_time,edf.fixations.endTime);
fixend=fixend*2;
fixind=fixstart+fixend;
fixs=fixind;

% Get blinks
blinkstart=ismember(gaze_time,edf.blinks.startTime);
blinkend=ismember(gaze_time,edf.blinks.endTime);
blinkend=blinkend*2;
blinkind=blinkend+blinkstart;
blinks=blinkind;

% Get saccades
saccstart=ismember(gaze_time,edf.saccades.startTime);
saccend=ismember(gaze_time,edf.saccades.endTime);
saccend=saccend*2;
saccind=saccend+saccstart;
saccs=saccind;

% Start of fixation is 9, fixation duration is 1, and end of the fixation is 2
for i=2:(length(gaze_time)-1);
	if fixs(i-1)==1 & fixs(i+1)<2;
		fixs(i)=1;
	end
end
for i=1:length(gaze_time);
	if fixend(i)==2;
		fixs(i-1)=1;
	elseif fixstart(i)==1;
		fixs(i)=9;
	end
end

% Start of blinks is 9, blink duration is 1, and end of the blink is 2
for i=2:(length(gaze_time)-1);
	if blinks(i-1)==1 & blinks(i+1)<2;
		blinks(i)=1;
	end
end
for i=1:length(gaze_time);
	if blinkend(i)==2;
		blinks(i-1)=1;
	elseif blinkstart(i)==1;
		blinks(i)=9;
	end
end

% Start of saccade is 9, saccade duration is 1, and end of the saccade is 2
for i=2:(length(gaze_time)-1);
	if saccs(i-1)==1 & saccs(i+1)<2;
		saccs(i)=1;
	end
end
for i=1:length(gaze_time);
	if saccs(i)==2;
		saccs(i-1)=1;
	elseif saccstart(i)==1;
		saccs(i)=9;
	end
end

e.eye.fixations = nan(edf.nTrials,ceil(maxTrialLen));
e.eye.blinks = nan(edf.nTrials,ceil(maxTrialLen));
e.eye.saccades = nan(edf.nTrials,ceil(maxTrialLen));

e.eye.event_time = (0:(size(e.eye.saccades,2)-1))/1000;

% go through each trial and populate traces
for iTrial = 1:edf.nTrials
  disppercent(iTrial/edf.nTrials);
  % the times for this trial
  thisTrialTimes = startTime(iTrial):1:endTime(iTrial);
  % get data for xPos, yPos and pupil traces form edf data 
  [c1]=find(gaze_time==(thisTrialTimes(1)));
	[c2]=find(gaze_time==(thisTrialTimes(end)));
	e.eye.fixations(iTrial,1:(c2-c1+1)) = fixs(c1:c2);
  e.eye.saccades(iTrial,1:(c2-c1+1)) = saccs(c1:c2);
  e.eye.blinks(iTrial,1:(c2-c1+1)) = blinks(c1:c2);
end
disppercent(inf);

% display figure
if dispFig
figure;
distance=sqrt((e.eye.xPos.^2)+(e.eye.yPos.^2));
meandist=nanmean(distance,1);
  yPos = meandist;
  xPos = e.eye.time;
  plot(xPos,yPos,'LineWidth',2);
  hold on
  cueappX=repmat(2,1,length(e.eye.time));
  plot(cueappX,yPos,'--r','LineWidth',1);
xlabel('Times (Secs)');
ylabel('Eye position (deg)');
title('Mean eye position by trial type');
ylim([0,10]);
xlim([0,5]);
x=0;end;
