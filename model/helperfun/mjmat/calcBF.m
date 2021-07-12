% Purpose:  This function will compute the Bayes factor and probability/evidence in favor of null hypothesis (H0) and 
%           the alternative hypothesis (H1). Based on Raftery 1995, described in Masson 2011.
%

function [bf, pH0, pH1] = calcBF(nObs,sse_effect,sse_error,effect_df)

% calculate partial eta squared
partialEtaSqr = sse_effect./sse_error;

% calculate delta BIC
deltaBIC = nObs*log(partialEtaSqr)+effect_df*log(nObs);

% calculate bayes factor
bf = exp(deltaBIC/2);

% calculate probabilities
pH0 = bf/(bf+1);
pH1 = 1-pH0;
