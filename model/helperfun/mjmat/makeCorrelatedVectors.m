% Make two correlated vectors with mean, mu, and standard deviation, sd.
% Covariance is specified by covarMat
% n = number of rows of vectors
% r = expected correlation coefficient

% U is a nx2 matrix where each column is an output vector

% Adapted from Hammad Munawar (www.hammadmunawar.weebly.com)
% IAA, Islamabad, Pakistan

function U = makeCorrelatedVectors(n,mu,sd,r,corrLimit,printOut)

% define defaults
if ieNotDefined('printOut')
    printOut = 1;
end

% correlations should always be within the expected range set by corrLimit
if ieNotDefined('corrLimit')
    corrLimit = 0.02;
end

% make covariance matrix
covarMat = eye(2); % make covariance matrix
covarMat(covarMat==0) = r;

verifyR = -inf;
while abs(abs(verifyR)-abs(r))>corrLimit
    L = chol(covarMat); % cholesky factorization
    U = mu+sd.*rand(n,2); % generate random vectors with specified mean and sd
    U=U*L; % set correlations
    verifyR = corr(U(:,1),U(:,2));
end

if printOut
    % Print output
    fprintf('\n mu = [%.3f %.3f] \n',mean(U));
    fprintf('sigma = [%.4f %.4f] \n', std(U));
    fprintf('r = %.3f \n',verifyR); % verify correlation
end
