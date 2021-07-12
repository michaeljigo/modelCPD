function y = fisherTrans(x)
% x is fisher transformed

y = 0.5*(log((1+x)./(1-x)));