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
%
%
% ----------------------------------------------------------------------
% Function created by Michael Jigo
% Last update : 2020-12-03
% ----------------------------------------------------------------------
% ----------------------------------------------------------------------
% RF=arbWeibull(params,x,threshPerformance,varargin)
% ----------------------------------------------------------------------
% Purpose:  Function defines a psychometric function adapted from a Weibull.
%           As a result, it only accepts positive values of the alpha parameter.
%           More information regarding alpha described in the "Inputs" section.
%           
%           The PF is arbitrary to the extent that its "base" depends on the 
%           inputted threshPerformance variable (see more info on 
%           threshPerformance below). Such a flexible "base" allows the 
%           adaptive procedure to target any performance level an experimenter
%           desires.
%
% ----------------------------------------------------------------------
% Input(s)
% params             : 1x4 vector defining parameters of the psychometric function
%                      [alpha beta gamma lambda]
%                      alpha  = location parameter
%                      beta   = slope
%                      gamma  = lower asymptote
%                      lambda = lapse rate (1-upper asymptote)
% x                  : 1xm vector of candidate alpha parameters that would achieve desired performance
% threshPerformance  : desired convergence performance level 
% varargin           : 'Inverse' (optional input), computes inverse of psychometric function
%
% ----------------------------------------------------------------------
% Output(s)
% y                  : probability of a correct response at each candidate alpha parameter
% ----------------------------------------------------------------------
% Function adapted by Michael Jigo from the Palamedes toolbox
% Last update: 02.19.21
% ----------------------------------------------------------------------


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
