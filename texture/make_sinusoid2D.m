% Purpose:  Generate a 2-dimensional grating. Grating image will always be square (i.e., width=height)
%
% By:       Michael Jigo
%           05.25.21

function [cosinusoid, sinusoid, imap] = make_sinusoid2D(varargin)

%% Default parameters for grating images 
in = {'freq'            % center SF of sinusoid
      'amp'             % amplitude
      'phase'           % phase
      'ori'             % orientation (0=vertical)
      'im_width'        % width of image
      'px_per_deg'      % pixels per degree
      'window_width'    % width of window function
      'window_type'};   % window type ('cosine' or 'gaussian')
      

val = {2                % center SF of sinusoid
       1                % amplitude
       0                % phase
       0                % orientation (0=vertical)
       4                % width of image
       32               % pixels per degree
       0                % width of window function
       'cosine'};       % window type ('cosine' or 'gaussian')

params = parseOptionalInputs(in,val,varargin);


%% Make grating
% calculate orientation
params.ori = deg2rad(params.ori);
a=cos(params.ori)*params.freq;
b=sin(params.ori)*params.freq;

% get a grid of x and y coordinates that has the correct number of pixels
pix_width = round(params.im_width*params.px_per_deg);
im_halfwidth = params.im_width/2;
[x,y] = meshgrid(linspace(-im_halfwidth,im_halfwidth,pix_width));
[ph,r] = cart2pol(x,y);
imap.x = x; 
imap.y = y;
imap.r = r;
imap.ph = ph;

% create grating
cosinusoid = params.amp*cos(2*pi*(a*x+b*y)+params.phase);
sinusoid = params.amp*sin(2*pi*(a*x+b*y)+params.phase);


%% Make window
if params.window_width>0
   space = linspace(-im_halfwidth,im_halfwidth,pix_width);
   switch params.window_type
      case 'cosine'
         % create cosine window function (i.e., Hann window)
         window = make_cosine_fun(r,0,params.window_width,1);
      case 'gaussian'
         window = make_gaussian(space,0,params.window_width,1);
         window = window'*window;
   end
else
   % no-windowing function
   window = ones(size(x));
end

%% Make windowed-grating
cosinusoid = cosinusoid.*window;
sinusoid = sinusoid.*window;
