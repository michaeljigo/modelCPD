% Purpose:  This function will create an asymmetric Gaussian with a plateau. There will be two sigma terms and a single exponent that will control
%           the curviness of the function. If the exponent is greater than 2, it will become more like a plateau. An exponent equal to 2 
%           corresponds to a normal distribution.

function [ag neg_ag pos_ag] = asymGaussian_twoPeak(x,center,sigmaLeft, sigmaRight, exponent, amplitude, dist, baseline)

% input can be either 6 scalars or 1 vector
if nargin==2
   p = center;
   center = p(1);
   sigmaLeft = p(2);
   sigmaRight = p(3);
   exponent = p(4);
   amplitude = p(5);
   dist = p(6);
   if length(p)==7
      baseline = p(7);
   else
      baseline = 0;
   end
elseif nargin<7
   error('Need to input a vector of 7 values or 7 separate inputs');
end
if ~exist('baseline','var')
   baseline = 0;
end

% evaluate function
negCenter = center-dist;
posCenter = center+dist;

pos_ag = nan(1,numel(x));
neg_ag = nan(1,numel(x));

neg_ag(x<=negCenter) = amplitude*exp(-(abs(x(x<=negCenter)-negCenter)./sigmaLeft.^2).^exponent)+baseline;
neg_ag(x>negCenter) = amplitude*exp(-(abs(x(x>negCenter)-negCenter)./sigmaRight.^2).^exponent)+baseline;

pos_ag(x<=posCenter) = amplitude*exp(-(abs(x(x<=posCenter)-posCenter)./sigmaLeft.^2).^exponent)+baseline;
pos_ag(x>posCenter) = amplitude*exp(-(abs(x(x>posCenter)-posCenter)./sigmaRight.^2).^exponent)+baseline;

ag = neg_ag+pos_ag;
