% Purpose:  Wrapper function for power analysis to inform Simran's experiment.
%           Input: # of trials to include (ideally a subset sampled with replacement)
%                  # of subjects to include (sampled with replacement)
%                  # of bootstrap iterations to run


function p = power4simran(n_trials,n_subj,n_boot)
addpath(genpath('~/apps'));

subjList = {'AS' 'DT' 'KT' 'MJ' 'MM' 'RF' 'SO' 'SP' 'SX' 'YS'};

parfor b = 1:n_boot
   % randomly choose (with replacement) a number of subjects from the list above
   subj = subjList(randi(numel(subjList),1,n_subj));

   % collect performance for each subject with the pre-defined # of (randomly sampled) trials
   boot_perf = nan(numel(subj),2,9,2);
   for s = 1:numel(subj)
      boot_perf(s,:,:,:) = compute_performance(subj{s},'endo',n_trials);
   end

   % perform repeated measures ANOVA and look at the main effect of eccentricity only
   rm_both_ecc = simple_mixed_anova(boot_perf(:,:,1:9,:));
   rm_single_ecc = simple_mixed_anova(boot_perf(:,:,1:9,1));

   % store p-value for the main effect of attention
   both_ecc(b) = rm_both_ecc{3,5};
   single_ecc(b) = rm_single_ecc{3,5};
   fprintf('Done %i/%i\n',b,n_boot);
end
p.both_ecc = both_ecc;
p.single_ecc = single_ecc;

filename = sprintf('./power_for_simran/n_subj%i.mat',n_subj);
save(filename,'p','n_boot','n_trials');


%% Compute subject performance
function perf = compute_performance(subj,whichAttn,n_trials,varargin)

%% Parse optional inputs
optionalIn = {'dataDir'};
optionalVal = {sprintf('../%s/data/%s/homo/',whichAttn,subj)};
opt = parseOptionalInputs(optionalIn,optionalVal,varargin);


%% Load subject's data
% parse stim files
data = parseFiles(opt.dataDir,1,{'tilt' 'cues' 'ecc' 'loc' 'sfs' 'accuracy' 'brokenTrial.trialIdx' 'response' 'reactionTime'});

% remove trials in which fixation was broken
data = structfun(@(x) x(~data.trialIdx(~isnan(data.trialIdx))),data,'UniformOutput',0);

% transform tilts and responses to target present (CW; 1) and target absent (CCW; 0)
data.tilt(data.tilt==-1) = 0;
data.response(data.response==1) = 0;
data.response(data.response==2) = 1;

% get tested SFs and ecc
sfVal = unique(data.sfs);
eccVal = unique(data.ecc);


%% Compute performance within sets
% create factor structures for:
% cue
cues.val = data.cues;
cues.label = {'neutral' 'valid'};

% sfs
sfs.val = data.sfs;
sfs.label = cellfun(@num2str,num2cell(unique(sfs.val)),'UniformOutput',0);

% ecc
ecc.val = data.ecc;
ecc.label = cellfun(@num2str,num2cell(unique(ecc.val)),'UniformOutput',0);

% create condition structures for target
target.val = data.tilt;
target.label = {'absent' 'present'};

% compute d-prime
perf = condParser(data.response,target,cues,sfs,ecc);

% sample, with replacement, the desired # of trials
perf.raw = cellfun(@(x) x(randi(numel(x),1,n_trials)),perf.raw,'uniformoutput',0);

% Hautas, 1995 approach
% add 0.5 to the # of hits and # of false alarms; add 1 to the # of signal and noise trials
% take proportion and then compute d-prime
fa_hits = cellfun(@sum, perf.raw)+0.5;
noise_sig = cellfun(@numel, perf.raw)+1;
perf.perf = fa_hits./noise_sig;

% compute Hautas d-prime
perf.perf = squeeze(diff(norminv(perf.perf),[],1));
perf = perf.perf;
perf(perf<0) = 0;
