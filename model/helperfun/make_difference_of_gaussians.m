% Purpose:  Generate difference-of-Gaussians (DoG) function.
%           Two Gaussians define the DoG, one with a positive amplitude and another with a negative amplitude.
%           Each Gaussian is defined independently from each other, but share the same center.
%
% By:       Michael Jigo
% Edited:   07.01.21

function [dog pos_gaus neg_gaus] = make_difference_of_gaussians(position,center,pos_amp,neg_amp,pos_width,neg_width)

   % make positive gaussian
      pos_gaus = make_gaussian(position,center,pos_width,pos_amp);
   
   % negative gaussian
      neg_gaus = make_gaussian(position,center,neg_width,neg_amp);

   % create difference-of-Gaussians
      dog      = pos_gaus-neg_gaus;
