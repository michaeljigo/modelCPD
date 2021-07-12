% Naka-Rushton formula (Michaelis-Menten equation):

% response = response_max*((intensity^n)/(intensity^n+half_sat^n))

% Where:
% out = response
% x = intensity
% maxResp = response_max
% halfResp = half_sat
% n = slope

function out = nakaRushton(x,params)

maxResp = params(1);
xAtHalfResp = params(2);
n = params(3);
if length(params)==4
   baseline = params(4);
elseif length(params)==3
   baseline  = 0;
else
   error('Params must be a 1x4 vector');
end

out = maxResp*((x.^n)./(x.^n+xAtHalfResp^n))+baseline;
