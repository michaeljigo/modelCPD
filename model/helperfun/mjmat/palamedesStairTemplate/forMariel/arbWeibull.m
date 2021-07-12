% arbWeibull   Evaluation of an arbritrary function; based on Weibull
% funciton. only accepts positive values of alpha.
%
%   syntax: y = arbWeibull(params, x, threshPerformance)
%
%   threshPerformance = the performance at which you wish to hold performance
%
%   y = arbWeibull(params, x, threshPerformance), where 'params' contains the four
%   parameters of a Psychometric Funtion (i.e., [alpha beta gamma lambda]),
%   returns the Psychometric Function evaluated at values in 'x'. 'x' is
%   array of any size.
%
%   x = arbWeibull(params, y, threshPerformance, 'Inverse') returns the x-value at
%   which the Psychometric Function evaluates to y.
%
%   dydx = arbWeibull(params, x, threshPerformance, 'Derivative') returns the
%   derivative (slope of tangent line) of the Psychometric Function
%   evaluated at x. WILL DO DERIVATION LATER
%
%   'params' need not have four entries. A two element vector will be
%   interpreted as [alpha beta], a three element vector as [alpha beta
%   gamma]. Missing elements in 'params' will be assigned a value of 0.

function y = arbWeibull(params,x,threshPerformance,varargin)

if min(x) < 0
    message = 'This function (based off the Weibull) is not defined for negative stimulus intensity values. ';
    message = [message 'In case your stimulus intensity values are log transformed values, you '];
    message = [message 'should use the Gumbel (or ''log-Weibull'') function (PAL_Gumbel).'];
    message = [message 'Visit www.palamedestoolbox.org/weibullandfriends.html for more '];
    message = [message 'information'];
    error('PALAMEDES:negativeXweibull',message);
end

[alpha, beta, gamma, lambda] = PAL_unpackParamsPF(params);

% make base of function
base = 1/(1-((threshPerformance-gamma)/(1-gamma-lambda)));

if ~isempty(varargin)
    if strncmpi(varargin{1}, 'Inverse',3)
        term1 = (log(1-lambda-x)-log(1-gamma-lambda))/log(base);
        term2 = ((1/alpha)^beta)^-1;
        y = (abs(term1*term2))^(1/beta);
    end
else
    y = gamma+(1-gamma-lambda).*(1-base.^(-1*(x./alpha).^beta));
end