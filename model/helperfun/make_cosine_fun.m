% Purpose:  This function creates a cosine function centered at "center", over the values defined by the 
%           vector "space". Also creates an annulus of raised cosines that tile space around the center.
% By:       Michael Jigo
%
% Inputs:   space    --    vector specifying spatial location
%           center   --    center position of the raised cosine
%           width    --    specifies the full width at half-maximum (i.e., bandwidth)
%           height   --    if specified, this is the amplitude of the cosine. Otherwise, it is scaled to unit
%                          volume.
% Outputs:  rc       --    vector containing cosine function centered at the inputted "center" location
%           annulus  --    vector containing cosine functions that are shifted from the center location but
%                          tile space such that the sum of rc and annulus will equal a plateau.

function [rc annulus] = make_cosine_fun(space,center,width,height,exponent)

width = width/2;
rc = 0.5+0.5*cos((space-center)./width*pi);
annulus = 0.5+0.5*cos((space-(center+width))./width*pi);
rc(space>=center+width | space<=center-width) = 0;
annulus(space>=center+2*width | space<=center-2*width) = 0;

if exist('exponent','var') && ~isempty(exponent)
   rc = rc.^exponent;
   annulus = annulus.^exponent;
end

if ~exist('height','var')
   rc = rc./sum(rc);
   annulus = annulus./sum(annulus);
else
   rc = rc*height;
   annulus = annulus*height;
end

