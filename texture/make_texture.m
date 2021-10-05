% ----------------------------------------------------------------------
% function [targ_tex notarg_tex texture_params]=make_texture(varargin)
% ----------------------------------------------------------------------
% Purpose:  Generate texture of oriented lines. One texture will contain
%           a subpatch of oriented lines differing from the background. 
%
%           If run without any input arguments, the defaults will be used.
% 
% ----------------------------------------------------------------------
% Input(s)
% px_per_deg         : pixels per visual degree in image 
% im_size            : [height,width] of final texture image
% line_size          : [width,height] dimension of individual lines (in degrees of visual angle, dva)
% line_jitter        : spatial jitter of individual line positions (in dva)
% line_spacing_row   : vertical spacing between lines (in dva)
% line_spacing_col   : horizontal spacing between lines (in dva)
% targ_array_size    : [rows,columns] in the target patch
% ori_targ           : orientation of lines in target patch
% ori_bg         : orientation of lines in background
% ori_bw             : bandwidth of oriented lines (larger bandwidth = more variable line tilts)
% targ_ecc           : horizontal eccentricity of target patch
% init_lines         : number of lines to initialize texture grid with (smaller numbers speed up texture generation, but must be large enough for texture size )
%
%
% e.g. generate_line_texture('px_per_deg',32,'im_size',[5 5],'line_size',[0.1 0.4],'line_jitter',0,'line_spacing_row',0.68,'line_spacing_col',0.71,'targ_array_size',[3 3],'ori_targ',-45 ,'ori_bg',45,'ori_bw',0,'targ_ecc',0)
%
% ----------------------------------------------------------------------
% Output(s)
% targ_tex           : texture image with target patch
% notarg_tex         : texture image without target patch
% texture_params     : parameters that define texture
% ----------------------------------------------------------------------
% Function created by Michael Jigo
% Last update : 2020-11-05
% ----------------------------------------------------------------------

function [targ_tex notarg_tex texture_params] = make_texture(varargin)

%% Set default texture parameters
in = {'px_per_deg' ...           % pixels per degree
      'im_size' ...              % image size in degrees
      'line_size' ...            % [width height] of line elements in degrees
      'line_jitter' ...          % spatial jitter between lines in degrees
      'line_spacing_row' ...     % spacing between rows in degrees
      'line_spacing_col' ...     % spacing between columns in degrees
      'targ_array_size' ...      % [rows columns] of target array
      'ori_targ' ...             % orientation of target array
      'ori_bg' ...           % orientation of surrounding array
      'ori_bw' ...               % orientaiton bandwidth of all lines
      'targ_ecc' ...             % eccentricity, in degrees, of target patch
      'init_lines' ...           % # of lines to initiliaze the target array
      'bg_color'};           % background color - 'white' or 'gray'

% default values
val = {32 ...                    % px_per_deg
      [20 20] ...                % im_size
      [0.1 0.4] ...              % line_size
      0 ...                      % line_jitter
      0.71 ...                   % line_spacing_row
      0.68 ...                   % line_spacing_col
      [3 3] ...                  % targ_array_size
      -45 ...                    % ori_targ
      45 ...                     % ori_bg
      0 ...                      % ori_bw
      0 ...                      % targ_ecc
      100 ...                    % init_lines
      'white'};                  % bg_color                   
      
texture_params = parseOptionalInputs(in,val,varargin);

%% Give warning about overlap
   % inputted line spacing
   texture_params.line_spacing_row = round(1e3*texture_params.line_spacing_row)./1e3;
   texture_params.line_spacing_col = round(1e3*texture_params.line_spacing_col)./1e3;

   % compute minimum spacing for absolutely no overlap
   no_overlap = sqrt(2*max(texture_params.line_size).^2);
   if any([texture_params.line_spacing_row texture_params.line_spacing_col] < no_overlap)
      warning('Line spacing is below %.2f. Lines will overlap, especially when ori_bw is high',no_overlap);
   end



%%%%%%% Create an initial large texture then cut the desired # of rows and columns out of the large texture
%% Generate line indices
n_lines = texture_params.init_lines;
n_lines = round(n_lines/2);
[line_row_idx line_col_idx targ_row targ_col jitter ori ori_targ] = get_line_indices(n_lines,texture_params);


%% Make vertical line with the specified dimensions
% line template (tmp) will be a square that can contain a 45-oriented line with a few pixels of extra blank space
tmp_size = ceil(sqrt(sum(texture_params.line_size.^2))*1e1)./1e1;
tmp_pix = ceil(texture_params.px_per_deg*tmp_size)+5;
tmp_line = zeros(tmp_pix,tmp_pix);

% get coordinates for line
tmp_cntr = ceil(median(1:tmp_pix));
line_dim = floor(texture_params.line_size*texture_params.px_per_deg);
col_idx = ceil(tmp_cntr-line_dim(1)/2):floor(tmp_cntr+line_dim(1)/2);
row_idx = ceil(tmp_cntr-line_dim(2)/2):floor(tmp_cntr+line_dim(2)/2);

% draw line
tmp_line(row_idx,col_idx) = 1;


%% Initialize texture grid
% pre-allocate matrix
pix_width = ceil(((max(line_col_idx)-min(line_col_idx))+tmp_size)*texture_params.px_per_deg);
pix_height = ceil(((max(line_row_idx)-min(line_row_idx))+tmp_size)*texture_params.px_per_deg);
notarg_tex = zeros(pix_height,pix_width);
targ_tex = notarg_tex;
cntr_px = round(size(notarg_tex)/2);


%% Generate texture grid
% convert center locations and sizes to pixels
pix_col_centers = line_col_idx*texture_params.px_per_deg;
pix_row_centers = line_row_idx*texture_params.px_per_deg;
pix_line = unique(size(tmp_line));
pix_jitter = jitter*texture_params.px_per_deg;

% re-center row and column centers in terms of pixels
pix_row_centers = pix_row_centers+cntr_px(1);
pix_col_centers = pix_col_centers+cntr_px(2);

x_size = nan(1,2);
y_size = nan(1,2);

% initialize target
for r = 1:numel(line_row_idx)
   for c = 1:numel(line_col_idx)
      tmp = [];
      
      % set row indices
      start_row = (pix_row_centers(r)+pix_jitter(r,c))-pix_line/2;
      end_row = (pix_row_centers(r)+pix_jitter(r,c))+pix_line/2-1;
      row_idx = round(start_row:end_row);

      % keep track of negative indices and remove corresponding values from line
      out_row_idx = row_idx<=0 | row_idx>=size(targ_tex,1);
      row_idx = row_idx(~out_row_idx);
   
      % start of removed line index
      if sum(out_row_idx)>0
         line_row_keep = find(out_row_idx==0);
         line_row_keep = max(1,min(line_row_keep)):max(line_row_keep);
      else
         line_row_keep = 1:size(tmp_line,1);
      end

      % set column indices
      start_col = (pix_col_centers(c)+pix_jitter(r,c))-pix_line/2;
      end_col = (pix_col_centers(c)+pix_jitter(r,c))+pix_line/2-1;
      col_idx = round(start_col:end_col);

      % keep track of negative indices and remove corresponding values from line
      out_col_idx = col_idx<=0 | col_idx>=size(targ_tex,2);
      col_idx = col_idx(~out_col_idx);

      % start of removed line index
      if sum(out_col_idx)>0
         line_col_keep = find(out_col_idx==0);
         line_col_keep = max(1,min(line_col_keep)):max(line_col_keep);
      else
         line_col_keep = 1:size(tmp_line,2);
      end

      % now that sizes are set, rotate line before adjusting line size
      bg_line = imrotate(tmp_line,ori(r,c),'nearest','crop');
      targ_line = imrotate(tmp_line,ori_targ(r,c),'nearest','crop');

      % adjust line size
      targ_line = targ_line(line_row_keep,line_col_keep);
      bg_line = bg_line(line_row_keep,line_col_keep);

      % insert target and background lines
      if targ_row(r) && targ_col(c)
         tmp(:,:,1) = targ_tex(row_idx,col_idx);
         tmp(:,:,2) = targ_line;
         tmp = sum(tmp,3,'omitnan');
         targ_tex(row_idx,col_idx) = targ_tex(row_idx,col_idx)+targ_line;
         
         % keep track of the dimensions of the target
         x_size(:,1) = min([x_size(:,1); col_idx']);
         x_size(:,2) = max([x_size(:,2); col_idx']);
         y_size(:,1) = min([y_size(:,1); row_idx']);
         y_size(:,2) = max([y_size(:,2); row_idx']);
      else
         tmp(:,:,1) = targ_tex(row_idx,col_idx);
         tmp(:,:,2) = bg_line;
         tmp = sum(tmp,3,'omitnan');
         targ_tex(row_idx,col_idx) = targ_tex(row_idx,col_idx)+bg_line;
      end

      % generate no-target grid
      notarg_tex(row_idx,col_idx) = notarg_tex(row_idx,col_idx)+bg_line; 
   end
end
texture_params.targ_patch_size = [diff(x_size) diff(y_size)]./texture_params.px_per_deg;


%% Now crop to desired image size
% cut out indices
im_halfwidth_px = texture_params.im_size/2*texture_params.px_per_deg;
row_px = [-im_halfwidth_px(1) im_halfwidth_px(1)]+cntr_px(1);
row_px = row_px(1):row_px(2)-1;
col_px = [-im_halfwidth_px(2) im_halfwidth_px(2)]+cntr_px(2);
col_px = col_px(1):col_px(2)-1;

% crop image
row_px = floor(row_px); col_px = floor(col_px);
targ_tex = targ_tex(row_px,col_px);
notarg_tex = notarg_tex(row_px,col_px);
return











%%%% Helper functions
function [line_row  line_col targ_cntr_row targ_cntr_col jitter ori ori_targ] = get_line_indices(n_lines,texture_params);

targ_rows = texture_params.targ_array_size(1);
if floor(targ_rows/2)*2 == targ_rows;
   % even # of target lines: center pixel will be centered between two lines
   line_row = sort([texture_params.line_spacing_row/2:texture_params.line_spacing_row:(texture_params.line_spacing_row*n_lines) -texture_params.line_spacing_row/2:-texture_params.line_spacing_row:(-texture_params.line_spacing_row*n_lines)]);
  
   % set target centers
   n_targ_per_side = targ_rows/2;
   targ_cntr_row = n_targ_per_side*[-texture_params.line_spacing_row texture_params.line_spacing_row];
   targ_cntr_row = line_row>=min(targ_cntr_row) & line_row<=max(targ_cntr_row);
else
   % odd # of target lines: center pixel will be on a line
   line_row = -(texture_params.line_spacing_row*n_lines):texture_params.line_spacing_row:(texture_params.line_spacing_row*n_lines);
   line_row = round(line_row*1e3)./1e3;
   
   % set target centers
   n_targ_per_side = (targ_rows-1)/2;
   targ_cntr_row = n_targ_per_side*[-texture_params.line_spacing_row texture_params.line_spacing_row];
   targ_cntr_row = line_row>=min(targ_cntr_row) & line_row<=max(targ_cntr_row);
end

% COLUMNS 
targ_cols = texture_params.targ_array_size(2);
if floor(targ_cols/2)*2 == targ_cols;
   % even # of target lines: center pixel will be centered between two lines
   line_col = sort([texture_params.line_spacing_col/2:texture_params.line_spacing_col:(texture_params.line_spacing_col*n_lines) -texture_params.line_spacing_col/2:-texture_params.line_spacing_col:(-texture_params.line_spacing_col*n_lines)]);
  
   % set target centers
   n_targ_per_side = targ_cols/2; % # of target elements on either side of central line
   targ_cntr_col = n_targ_per_side*[-texture_params.line_spacing_col texture_params.line_spacing_col];
   targ_cntr_col = round(targ_cntr_col.*1e6)./1e6;
   targ_cntr_col = line_col>=min(targ_cntr_col) & line_col<=max(targ_cntr_col);

   % shift target eccentricity based on column spacing
   numShift    = round(texture_params.targ_ecc./texture_params.line_spacing_col);
   newCntrCol  = zeros(1,numel(targ_cntr_col));   
   newCntrCol(find(targ_cntr_col)+numShift) = 1;
   targ_cntr_col = newCntrCol;
else
   % odd # of target lines: center pixel will be on a line
   line_col = -(texture_params.line_spacing_col*n_lines):texture_params.line_spacing_col:(texture_params.line_spacing_col*n_lines);
   line_col = round(line_col*1e3)./1e3;
   
   % set target centers
   [~,new_cntr] = min(abs(line_col-texture_params.targ_ecc)); % center target on desired (horizontal) eccentricity
   n_targ_per_side = (targ_cols-1)/2; % compute # of target elements on either side of central line
   targ_cntr_col = n_targ_per_side*[-texture_params.line_spacing_col texture_params.line_spacing_col]+line_col(new_cntr);
   targ_cntr_col = round(targ_cntr_col.*1e6)./1e6;
   targ_cntr_col = line_col>=min(targ_cntr_col) & line_col<=max(targ_cntr_col);
end


%% Add jitter 
% spatial jitter
jitter = unifrnd(-texture_params.line_jitter/2,texture_params.line_jitter/2,[numel(line_row) numel(line_col)]);

% orientation jitter
ori = texture_params.ori_bg+unifrnd(-texture_params.ori_bw/2,texture_params.ori_bw/2,[numel(line_row) numel(line_col)]);
ori_targ = texture_params.ori_targ+unifrnd(-texture_params.ori_bw/2,texture_params.ori_bw/2,[numel(line_row) numel(line_col)]);


   


function opt = parseOptionalInputs(names,val,updateVals)

% if updateVals is a cell array, that means that varargin was used as input in a nested function. expand cell array
if numel(updateVals)==1 && iscell(updateVals)
   updateVals = updateVals{1};
end

% check for length between default names and values
if ~isequal(numel(names),numel(val))
   error('WHOOPS! Mismatch in number of names annd number of default values.')
end

% align actual inputs with the possible inputs
[~,optIdx] = ismember(names,updateVals(1:2:end));
[~,userIdx] = ismember(updateVals(1:2:end),names);

% remove names (in updateVals) that are not part of the defaults specified in "names"
invalidIn = find(userIdx==0)*2;
updateVals([invalidIn invalidIn-1]) = [];
userIdx(userIdx==0) = [];

if isempty(userIdx)
   % if none of the to-be-updated values are not included, just set the defaults
   for i = 1:length(names)
      opt.(names{i}) = val{i};
   end
else
   % if some relevant parameters are inputted, update the corresponding values
   userIn = cell(1,length(names)); userIn(userIdx) = updateVals(2:2:end);
   userIn(optIdx==0) = val(optIdx==0);

   % create the input variables
   for i = 1:length(names)
      opt.(names{i}) = userIn{i};
   end
end
