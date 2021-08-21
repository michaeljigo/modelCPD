% Purpose:  Initialize parameters for model optimization.
%           Model parameters include:
%           stimdrive      stimulus drive parameters
%           supdrive       suppressive drive parameters
%           attn           attentional gain parameters
%
% Input:    varargin             specify desired parameter name and its value in matching pairs (e.g., init_model_parameters('cg_max',1))
%              attn_neg_spread   only useful if difference-of-gaussians (dog) desired
%              attn_neg_amp      only useful if difference-of-gaussians desired
%              attn_amp_baseline unused for manuscript
%              spatial_profile   'cosine' or 'dog'
%              fix_with_input    fixes input parameters during optimization
%
% Output:   stimdrive            structure of stimulus drive parameters
%           supdrive             structure of suppressive drive parameters
%           attn                 structure of attentional gain parameters
%
% By:       Michael Jigo


function [stimdrive supdrive attn] = init_parameters(varargin) 

%% Initialize list of possible parameters
in = {'cg_max'
      'cg_slope'
      'freq_max' 
      'freq_slope'
      'freq_min'
      'bw_max' 
      'bw_slope'
      'sup_space'
      'sup_ori' 
      'sup_freq' 
      'attn_freq_max'
      'attn_freq_slope'
      'attn_bw' 
      'attn_amp_max'
      'attn_amp_slope'
      'attn_spread' 
      'attn_baseline'
      'attn_sup_amp'
      'attn_sup_spread'
      'sf_profile'
      'stimdrive'          % stimulus drive struct; start points, bounds
      'supdrive'           % suppressive drive struct; start points, bounds
      'attn'               % attention struct; start points, bounds
      'fix_with_input'};
val = [repmat({[]},1,numel(in)-1) 0];
start_params = parseOptionalInputs(in,val,varargin);



%% Stimulus drive and contrast gain parameters
% cg_max       maximum contrast gain (at fovea)
% cg_slope     change in contrast gain with eccentricity
% freq_max     maximum SF preference (at fovea)
% freq_slope   shift in SF preference with eccentricity
% bw_max       SF bandwidth
% bw_slope     change in SF bandwidth with eccentricity (fixed at 0)
names      = { 'cg_max'    'cg_slope'       'freq_max'   'freq_slope'   'freq_min'     'bw_max'    'bw_slope'};
bnds       = { [1 3]        [-1 -0.05]       [0 3.5]      [-1 0]        [0.5 0.5]      [1 2]       [0 0]};
plaus_bnds = { [1 2]        [-0.3 -0.1]      [2 3]        [-0.7 -0.6]   [0.5 0.5]      [1 2]       [0 0]};


if isempty(start_params.stimdrive)
   % put parameters into stimdrive structure
   for n = 1:numel(names)
      stimdrive.bnd.(names{n}) = bnds{n};
      stimdrive.plaus_bnd.(names{n}) = plaus_bnds{n};
   
      if isempty(start_params.(names{n}))
         % randomly choose starting parameters
         stimdrive.(names{n}) = rand(1)*diff(bnds{n})+min(bnds{n});
      else
         % use the inputted start parameters
         stimdrive.(names{n}) = start_params.(names{n});
   
         if start_params.fix_with_input
            % this option will adjust the bounds so that the parameter value will remain fixed at what's inputted
            stimdrive.bnd.(names{n}) = repmat(start_params.(names{n}),1,2);
            stimdrive.plaus_bnd.(names{n}) = repmat(start_params.(names{n}),1,2);
         end
      end
   end
else
   stimdrive = start_params.stimdrive;
end
stimdrive.param_names = names;



%% Suppressive drive parameters
% sup_space    surround suppression pool width; value is a scalar defining how width relates to channel SF preference
% sup_ori      cross-orientation suppression pool; value is in degrees
% sup_freq     cross-frequency suppression pool
names =        {'sup_space'   'sup_ori'   'sup_freq'};
bnds =         {[4 4]         [180 180]   [2 2]};
plaus_bnds =   {[4 4]         [180 180]   [2 2]};

if isempty(start_params.supdrive)
   % put parameters into supdrive structure
   for n = 1:numel(names)
      supdrive.bnd.(names{n}) = bnds{n};
      supdrive.plaus_bnd.(names{n}) = plaus_bnds{n};
   
      if isempty(start_params.(names{n}))
         % randomly choose starting parameters
         supdrive.(names{n}) = rand(1)*diff(bnds{n})+min(bnds{n});
      else
         % use inputted start parameters
         supdrive.(names{n}) = start_params.(names{n});
         
         if start_params.fix_with_input
            % this option will adjust the bounds so that the parameter value will remain fixed at what's inputted
            supdrive.bnd.(names{n}) = repmat(start_params.(names{n}),1,2);
            supdrive.plaus_bnd.(names{n}) = repmat(start_params.(names{n}),1,2);
         end
      end
   end
else
   supdrive = start_params.supdrive;
end
supdrive.param_names = names;



%% Attention parameters
% attn_freq_max      maximum attentional SF preference (at fovea)
% attn_freq_slope    change in attentional SF preference with eccentricity (in octaves)
% attn_bw            attention SF bandwidth (in octaves)
% attn_amp_max       maximum attentional amplitude
% attn_amp_slope     change in attentional amplitude with eccentricity (fixed at 0)
% attn_spread        spatial spread of attention (in degrees)
% attn_baseline      baseline attention modulation
% attn_sup_amp       scalar of suppressive surround amplitude (difference-of-Gaussians only)
% attn_sup_spread    scalar of suppressive surround spread (difference-of-Gaussians only)

names =      {'attn_freq_max'    'attn_freq_slope'     'attn_bw'     'attn_amp_max'   'attn_amp_slope'   'attn_spread'     'attn_baseline'      'attn_sup_amp'   'attn_sup_spread'};
bnds =       {[1 4]              [-1 -0.1]            [1 3]         [1 20]            [0 0]             [1 5]             [0.1 1]              [0.1 0.95]          [1 5]};
plaus_bnds = {[2 3]              [-0.2 -0.1]           [1 2]         [1 3]             [0 0]             [3 4]             [0.5 1]              [0.5 0.8]        [2 4]};


% put parameters into attn structure
if isempty(start_params.attn)
   for n = 1:numel(names)
      attn.bnd.(names{n}) = bnds{n};
      attn.plaus_bnd.(names{n}) = plaus_bnds{n};

      if isempty(start_params.(names{n}))
         % randomly choose starting parameters
         attn.(names{n}) = rand(1)*diff(bnds{n})+min(bnds{n});
      else
         % use inputted start parameters
         attn.(names{n}) = start_params.(names{n});
         
         if start_params.fix_with_input
            % this option will adjust the bounds so that the parameter value will remain fixed at what's inputted
            attn.bnd.(names{n}) = repmat(start_params.(names{n}),1,2);
            attn.plaus_bnd.(names{n}) = repmat(start_params.(names{n}),1,2);
         end
      end
   end
else
   attn = start_params.attn;
end
attn.param_names = names;
