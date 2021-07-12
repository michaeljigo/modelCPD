% Purpose:  This function will run simulations to determine the necessary sample size to obtain a power of 80% for the main effect
%           of attention across various spatial frequencies. Power will be determined using a bootstrap approah:
%           1. Sample, with replacement, n subjects from the available pool of subjects that ran the experiment.
%           2. Perform ANOVA on resampled dataset.
%           3. Repeat Steps 1 & 2 a large number of times.
%           4. Count the number of times the main effects and interactions were significant.
%           
%           The desired # of subjects should produce significant effects 80% of the time (i.e., power=0.8).
%
% By:       Michael Jigo
% Edited:   06.22.21
%
% Input:    dataPath    string for the path to the saved data for each experiment
%           nSubj       # of subjects to simulate
%           nSim        # of repetitions to simulate
%
% Output:   pval        structure with fields for each statsitical comparison in the ANOVA

function pval = simulatePower_bootstrapSubj(dataPath,nSubj,nSim)

%% Load data
load(dataPath);
cs = log(cs);

%% Resample, with replacement, nSubj from the pool of subjects available
tic
for s = 1:nSim
   rng('shuffle');
   sampleIdx = randi(size(cs,1),1,nSubj);
   if numel(size(cs))==2
      simCS = cs(sampleIdx,:);

      % perform one-way repeated measures anova
      anovaOut = simple_mixed_anova(simCS);
      pval.attn(s) = anovaOut{3,6};
   elseif numel(size(cs))==3
      simCS = cs(sampleIdx,:,:);
   
      % perform repeated-measures anova and store p-values for main effects and interactions
      anovaOut = simple_mixed_anova(simCS);
      pval.attn(s) = anovaOut{3,6};
      pval.sfs(s) = anovaOut{5,6};
      pval.attn_sfs(s) = anovaOut{7,6};
   elseif numel(size(cs))==4
      simCS = cs(sampleIdx,:,:,:);

      anovaOut = simple_mixed_anova(simCS);
      pval.attn(s) = anovaOut{3,6};
      pval.sfs(s) = anovaOut{5,6};
      pval.ecc(s) = anovaOut{7,6};
      pval.attn_sfs(s) = anovaOut{9,6};
      pval.attn_ecc(s) = anovaOut{11,6};
      pval.sfs_ecc(s) = anovaOut{13,6};
      pval.attn_sfs_ecc(s) = anovaOut{15,6};
   end
end
toc
