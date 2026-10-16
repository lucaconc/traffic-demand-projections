function [mu,info] = fitPriorGLS(A,Y,prior,priorCV)
%FITPRIORGLS  Generalized least squares with a prior demand (Cascetta, 1984; Bell, 1991).
%
%   mu = fitPriorGLS(A,Y,prior,priorCV)
%
%   Minimizes, over mu >= 0,
%       (ybar - A*mu)' * Wy^-1 * (ybar - A*mu) + (mu - prior)' * Wp^-1 * (mu - prior),
%   where ybar are the mean sensor counts, Wy = diag(s_i^2/N) their estimated
%   variances (sample variances, floored at the mean counts' Poisson variance and
%   at 1/N), and Wp = diag((priorCV*prior_j)^2) the prior variances, i.e. a
%   multiplicative prior error of relative size priorCV. Stacking both terms gives
%   a nonnegative least-squares problem solved by lsqnonneg. Unlike the other
%   estimators, it needs a prior demand; with an exact prior it returns it.
[N,m] = size(Y); n = size(A,2);
assert(size(A,1) == m,'fitPriorGLS:Size','A must have one row per sensor.');
prior = double(prior(:)); assert(numel(prior) == n,'fitPriorGLS:Prior','prior must have one entry per class.');
ybar = mean(Y,1).';
vy = max(var(Y,0,1).', max(ybar,1)) / N;               % variance of the mean counts
vp = max((priorCV*prior).^2, 1e-8*max(1,max(prior))^2);  % prior variances, kept positive
M = [A ./ sqrt(vy); diag(1./sqrt(vp))];
d = [ybar ./ sqrt(vy); prior ./ sqrt(vp)];
mu = lsqnonneg(M,d);
info = struct('ybar',ybar,'countVariance',vy,'priorVariance',vp, ...
    'countMisfit',sum((ybar - A*mu).^2 ./ vy),'priorMisfit',sum((mu - prior).^2 ./ vp));
end
