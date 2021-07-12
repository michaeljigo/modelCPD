% Purpose:        Compute negative log-likelihood from predicted model response and observed data.
% By:             Michael Jigo
%                 04.07.21
% Adapted from:   Luigi Acerbi 
%                 (https://github.com/lacerbi/vbmc/wiki#i-do-not-have-a-likelihood-function-but-another-loss-function-can-i-still-use-vbmc)

function ll = log_likelihood(pred_y,data_x,data_y)
   %  compute sigma (std(pred_y-data_y))
   sigma = log(std(pred_y-data_y));
   sigma2 = exp(2*sigma);

   % # of independently tsested levels
   N = numel(data_x);

   % compute likelihood
   ll = -0.5*sum((pred_y - data_y).^2)/sigma2 - 0.5*N*log(2*pi*sigma2);
