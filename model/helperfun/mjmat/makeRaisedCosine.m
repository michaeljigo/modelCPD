function rc = makeRaisedCosine(space,center,width,height)
% Usage:    rc = makeRaisedCosine(space,center,width,height)
% Purpose:  This function creates a raised cosine centered at "center", over the values defined by the 
%           vector "space".
% Inputs:   space    --    vector specifying spatial location
%           center   --    center position of the raised cosine
%           width    --    specifies how wide the cosine lobe from minima to minima
%           height   --    if specified, this is the amplitude of the cosine. Otherwise, it is scaled to unit
%                          volume.


rc = 0.5+0.5*cos((space-center)/width*pi);
rc(space>=center+width | space<=center-width) = 0;

if ~exist('height','var')
   rc = rc./sum(rc);
else
   rc = rc*height;
end
