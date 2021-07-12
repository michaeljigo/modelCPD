% Purpose:  Create a landolt square with specified width and side length.
%           If function is run without any inputs, a default Landolt square will be generated.
%
% By:       Michael Jigo


function [im params] = make_landolt_square(varargin)

%% Set default texture parameters
in = {'px_per_deg' ...           % pixels per degree
      'line_length' ...          % image size in degrees
      'line_width' ...           % width of line elements in degrees
      'gap_size' ...             % length of gap (in arc minutes)
      'im_size' ...              % final size of image
      'display' ...              % display square
      'background_color'};       % background color - 'white' or 'gray'

% default values
val = {32 ...                    % px_per_deg
      1 ...                      % line_length
      0.05 ...                   % line_width
      20 ...                     % gap_size
      1 ...                      % im_size
      0 ...                      % display
      'black'};                  % backgrnd_color                   
      
params = parseOptionalInputs(in,val,varargin);


%% initialize image matrix
   px_size  = params.px_per_deg*params.line_length;
   im       = zeros(px_size);


%% create a square with the proper line width
   px_line_width = round(params.px_per_deg*params.line_width); % this is the # of pixels to leave untouched on each side of the square
   empty_idx = (px_line_width+1):(px_size-px_line_width);
   im(empty_idx,empty_idx) = 1;


%% add gap to the square
   % convert from arc minutes to degrees
   deg_gap = params.gap_size/60;
   px_gap = round(deg_gap*params.px_per_deg);

   % add gap to top of the square
   gap_row_idx = 1:px_line_width;
   col_cntr = round(px_size/2);
   gap_halflen = round(px_gap/2);
   startidx = col_cntr-gap_halflen+1;
   endidx = startidx+px_gap-1;
   im(gap_row_idx,startidx:endidx) = 1;


%% pad
   pad_deg = params.im_size-params.line_length; % # of degrees to add to image
   pad_px = round((params.px_per_deg*pad_deg)/2);
   im = padarray(im,pad_px,1);
   im = padarray(im',pad_px,1);
   im = im';


%% change background color
   if ismember(params.background_color,'black')
      % the square will be white and the background will be black
      newim = im;
      newim(im==0) = 1;
      newim(im==1) = 0;
      im = newim;
   end


%% display landolt square
   if params.display
      figure; imagesc(im); colormap gray; axis square; set(gca,'visible','off');
   end
