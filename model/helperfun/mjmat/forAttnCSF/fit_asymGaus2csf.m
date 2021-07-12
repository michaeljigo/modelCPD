% Purpose:  This function will fit the asymmetric Gaussian function to contrast sensivitity at various spatial frequencies.

function csfParams = fit_asymGaus2csf(sfVal,cs)

% set up solver
options = optimoptions(@fmincon,'Algorithm','active-set','Display','none','MaxIter',1e5,...
   'FunctionTolerance',1e-6,'StepTolerance',1e-50,'MaxFunctionEvaluations',1e50);
lb = [log2(0.5) 1 1 1];
ub = [log2(8) 5 5 500];
init = rand(1,numel(lb)).*(ub-lb)+lb;
fitFun = @(params)fit_asymGaus(sfVal,cs,params);
problem = createOptimProblem('fmincon','objective',fitFun,'x0',init,'lb',lb,'ub',ub,'options',options);
ms = MultiStart('FunctionTolerance',1e-6,'XTolerance',1e-3,'Display','off','StartPointsToRun','bounds','UseParallel',1);
csfParams = run(ms,problem,200);

% objective function for fitting double-exponential function
function cost = fit_asymGaus(sfVal,thresh,p)

if size(thresh,1)>1
   thresh = thresh';
end

% get dexp model output
newP = nan(1,6);
newP(1:3) = p(1:3);
newP(4) = 2;
newP(5) = p(4);
newP(6) = 1;
p = newP;

model = evalCSF('asymGaus',log2(sfVal),p);

% compute cost
cost = sqrt(nanmean((log(model)-log(thresh)).^2));
