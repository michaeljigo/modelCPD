% performs leave-one-out cross validation on a polynomial fit of order p
% and computes the mean squared error (cv) across the left-out samples

% p can be a vector and cv will be computed for each order of the
% polynomial

function cv = LOOCV_polynomial(x,y,p,dispRes)
if ieNotDefined('dispRes')
    dispRes = 1;
end

% initialize matrix that will hold result
cv = nan(length(p),length(x));
for pp = 1:length(p)
    for i = 1:length(x)
        % compute index for the non-left out samples
        idx = setdiff(1:length(x),i);
        % do fit with subset of data
        fitCoeff = polyfit(x(idx),y(idx),p(pp));
        % compute error of test
        cv(pp,i) = sum((y(i)-polyval(fitCoeff,x(i))).^2);
    end
end
% compute the mean squared error of the test samples
cv = mean(cv,2);

% if a vector of polynomial orders were inputted, plot the cv for each
% polynomial order
if length(p)>1 && dispRes
    figure('Name','Cross validation result');
    scatter(p,cv,50,'k','filled'); hold on
    plot(p,cv,'k-');
    xlabel('polynomial order');
    ylabel('CV error');
end

