% compute log-likelihood

function ll = computeLL(p,nCorr,nTrial,computeL)
if ieNotDefined('computeL')
    computeL = 0;
end
nCorr = double(nCorr);
nTrial = double(nTrial);
p = double(p);

if any(p==1)
    p(p==1) = 1-1e-6;
end
p(p==0) = 1e-6;

if computeL
    % compute the likelihood
    ll = prod((p.^nCorr).*(1-p).^(nTrial-nCorr));
else
    ll = real(sum(nCorr.*log(p) + (nTrial-nCorr).*log(1-p)));
end
