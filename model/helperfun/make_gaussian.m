% Purpose:  This function creates a Gaussian centered at "center", over the values defined by the 
%           vector "space".
% Inputs:   space    --    vector specifying spatial location
%           center   --    center position of the raised cosine
%           width    --    specifies the full-width-at-half-maximum (FWHM) 
%           height   --    if specified, this is the amplitude of the Gaussian. Otherwise, it is scaled to unit
%                          volume.
% Outputs:  gaus     --    vector containing Gaussian centered at the inputted "center" location

function gaus = make_gaussian(space,center,width,height)

gaus = exp((-4*log(2)*(space-center).^2)./width.^2);

if ~exist('height','var')
   gaus = gaus./sum(gaus);
else
   gaus = gaus*height;
end

