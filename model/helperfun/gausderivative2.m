% Purpose:  Generate 2nd-derivative of Gaussian function with specified parameters.

%function y = gausderivative2(position,amplitude,width,center,baseline)
function y = gausderivative2(position,center,width,amplitude,baseline)

%y = (2*amplitude)./(sqrt(3*width).*pi^0.25).*exp(-((position-center).^2./(2*width.^2))).*(1-((position-center).^2./width.^2))+baseline;

y = amplitude*exp(-((position-center).^2./(2*width.^2))).*(1-((position-center).^2./width.^2))+baseline;
