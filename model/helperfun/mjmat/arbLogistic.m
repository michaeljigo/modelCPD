%PAL_arbLogistic   Evaluation of an arbitrary Logistic Psychometric
%   Function. The "arbitrary" refers to the fact that the location (aka threshold) 
%   parameter is flexible. It does not have to correspond to the midpoint between
%   the upper and lower asymptotes. Instead, it can be any arbitrary proportion. 
%   This is particularly useful when trying to titrate performance in a task to a 
%   specific value of proportion correct.
%
%   syntax: y = arbLogistic(params, x, perfLevel)
%
%   y = arbLogistic(params, x, perfLevel), where 'params' contains the four
%   parameters of a Psychometric Funtion (i.e., [alpha beta gamma lambda]),
%   returns the Psychometric Function evaluated at values in 'x'. 'x' is
%   array of any size. Note that beta is the inverse of the normal
%   distribution's standard deviation (or 'sigma'). 'perfLevel' is a scalar that 
%   determines the proportion correct when evaluated at the location parameter.
%
%   'params' need not have four entries. A two element vector will be
%   interpreted as [alpha beta], a three element vector as [alpha beta
%   gamma]. Missing elements in 'params' will be assigned a value of 0.
%
%   This example returns the function value at threshold when gamma 
%   ('guess-rate') and lambda ('lapse-rate') both equal 0:
%       
%   y = PAL_CumulativeNormal([1 2 0 0], 1) returns:
%
%   y = 0.5000
%
%Introduced: Palamedes version 1.0.0 (NP)
%Modified: Palamedes version 1.0.2, 1.1.1, 1.2.0, 1.4.0, 1.4.4, 1.6.3 
%   (see History.m)

function y = arbLogistic(params, x, perfLevel, varargin)

if ~exist('perfLevel','var')
   perfLevel = 0.75;
end

[alpha, beta, gamma, lambda] = PAL_unpackParamsPF(params);

%% make function evaluate to perfLevel when x=alpha
perfLevel = (perfLevel-gamma)./(1-gamma-lambda);
% adjust alpha parameter to reflect desired performance level
alpha = alpha-(1./beta).*log(perfLevel./(1-perfLevel));

%% evaluate logistic function
y = gamma + (1 - gamma - lambda).*(1./(1+exp(-1*(beta).*(x-alpha))));
