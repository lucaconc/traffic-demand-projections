function [mu,info] = fitBootstrapGLS(A,Y,order,B,rs)
%FITBOOTSTRAPGLS Bootstrap-weighted nonnegative GLS on joint cumulants.
%
%   [mu,info] = fitBootstrapGLS(A,Y,order,B,rs)
%
%   Minimizes (k - C mu)' W (k - C mu) over mu >= 0, where k are the sample
%   joint cumulants through ORDER, C their linear system, and W the inverse of
%   the bootstrap covariance of k (B resamples of the intervals, ridge 1e-8 of
%   its mean variance). Solved as NNLS after the Cholesky factor of W.
if nargin < 5 || isempty(rs), rs = RandStream('mt19937ar','Seed',1); end
[C,subsets] = cumulantSystem(A,order);
k = sampleJointCumulants(Y,subsets);
N = size(Y,1); K = numel(k); boot = zeros(B,K);
for b = 1:B
    idx = randi(rs,N,N,1);
    boot(b,:) = sampleJointCumulants(Y(idx,:),subsets).';
end
S = cov(boot);
S = S + 1e-8*max(mean(diag(S)),eps)*eye(K);
[R,p] = chol(inv(S));
if p > 0, R = diag(1./sqrt(max(diag(S),eps))); end   % fallback: diagonal weights
mu = lsqnonneg(R*C,R*k,optimset('Display','off'));
info = struct('rank',rank(C),'equations',K);
end
