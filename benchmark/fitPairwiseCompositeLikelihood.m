function [mu,info] = fitPairwiseCompositeLikelihood(A,Y,mu0,varargin)
%FITPAIRWISECOMPOSITELIKELIHOOD Pairwise composite likelihood of sensor pairs.
%
%   Maximizes the sum over sensor pairs (i,k) of the exact bivariate log
%   likelihood of (Y_i,Y_k) = (U_i+W, U_k+W), with independent Poisson
%   U_i (paths crossing i but not k), U_k and W (paths crossing both), using
%   the analytic gradient and fmincon with bounds mu >= 1e-8. Its identifying
%   system is the cumulant system through order two.
%
%   Name/value: MaxIterations (200), MaxSeconds (Inf): the optimizer stops when
%   the time is exceeded and info.timeLimited is true (the estimate is then not
%   converged).
ip = inputParser; addParameter(ip,'MaxIterations',200); addParameter(ip,'MaxSeconds',Inf);
parse(ip,varargin{:}); oo = ip.Results; t0 = tic; timeLimited = false;
A = double(A~=0); [m,n] = size(A);
pairs = nchoosek(1:m,2);
% Pairs that no path crosses jointly have independent counts, so their bivariate
% log-likelihood is the sum of two Poisson marginals: they are evaluated in closed
% form, once per sensor with multiplicity deg(i) (same objective, much faster).
P = struct('both',{},'onlyI',{},'onlyK',{},'u',{},'w',{}); deg = zeros(m,1);
for t = 1:size(pairs,1)
    i = pairs(t,1); k = pairs(t,2);
    if ~any(A(i,:)&A(k,:)), deg(i) = deg(i)+1; deg(k) = deg(k)+1; continue; end
    [u,~,g] = unique(double(Y(:,[i k])),'rows');
    q = numel(P)+1;
    P(q).both = A(i,:)&A(k,:); P(q).onlyI = A(i,:)&~A(k,:); P(q).onlyK = A(k,:)&~A(i,:);
    P(q).u = u; P(q).w = accumarray(g,1);
end
Mg = struct('A',A,'deg',deg,'S',sum(double(Y),1).','G',sum(gammaln(double(Y)+1),1).','N',size(Y,1));
obj = @(x) negPairwise(x,P,Mg);
lb = 1e-8*ones(n,1);
x0 = max(double(mu0(:)),0.05*max(mean(double(mu0)),1));
opts = optimoptions('fmincon','Algorithm','interior-point','SpecifyObjectiveGradient',true, ...
    'Display','off','MaxIterations',oo.MaxIterations,'OptimalityTolerance',1e-8,'StepTolerance',1e-10, ...
    'OutputFcn',@stopOnTime);
[mu,~,flag] = fmincon(obj,x0,[],[],[],[],lb,[],[],opts);
info = struct('rank',rank(cumulantSystem(A,2)),'exitflag',flag,'timeLimited',timeLimited,'seconds',toc(t0), ...
    'coveredPairs',numel(P),'pairs',size(pairs,1));
    function stop = stopOnTime(~,~,~)
        stop = toc(t0) > oo.MaxSeconds;
        if stop, timeLimited = true; end
    end
end

function [f,g] = negPairwise(x,P,Mg)
f = 0; g = zeros(numel(x),1);
% uncovered pairs: deg(i) copies of sensor i's Poisson marginal log-likelihood
lam = Mg.A*x(:); use = Mg.deg > 0;
for i = find(use).'
    if lam(i) <= 0
        if Mg.S(i) > 0, f = Inf; end
        continue
    end
    f = f - Mg.deg(i)*(Mg.S(i)*log(lam(i)) - Mg.N*lam(i) - Mg.G(i));
    if nargout > 1, g = g - Mg.deg(i)*(Mg.S(i)/lam(i) - Mg.N)*Mg.A(i,:).'; end
end
for t = 1:numel(P)
    lb_ = sum(x(P(t).both)); li = sum(x(P(t).onlyI)); lk = sum(x(P(t).onlyK));
    y1 = P(t).u(:,1); y2 = P(t).u(:,2); w = P(t).w;
    l0 = bivLogPmf(y1,y2,lb_,li,lk);
    f = f - sum(w.*l0);
    if nargout > 1
        dB = exp(bivLogPmf(y1-1,y2-1,lb_,li,lk)-l0) - 1;
        dI = exp(bivLogPmf(y1-1,y2,lb_,li,lk)-l0) - 1;
        dK = exp(bivLogPmf(y1,y2-1,lb_,li,lk)-l0) - 1;
        g = g - (sum(w.*dB)*P(t).both(:) + sum(w.*dI)*P(t).onlyI(:) + sum(w.*dK)*P(t).onlyK(:));
    end
end
end

function lp = bivLogPmf(y1,y2,lb,li,lk)
% log P(U_i+W=y1, U_k+W=y2), W~Pois(lb), U_i~Pois(li), U_k~Pois(lk); -Inf if y<0.
% Vectorized over the rows; the sum over w uses a running log-sum-exp.
y1 = y1(:); y2 = y2(:); lp = -inf(size(y1));
ok = y1 >= 0 & y2 >= 0;
if ~any(ok), return; end
a = y1(ok); b = y2(ok); wmax = max(min(a,b));
acc = -inf(size(a));
for w = 0:wmax
    t = poisLog(w,lb) + poisLog(a-w,li) + poisLog(b-w,lk);
    t(w > min(a,b)) = -inf;
    mx = max(acc,t); fin = isfinite(mx);
    acc(fin) = mx(fin) + log(exp(acc(fin)-mx(fin)) + exp(t(fin)-mx(fin)));
end
lp(ok) = acc;
end

function v = poisLog(k,lam)
% log Poisson probability; -Inf for negative k. gammaln is evaluated only on
% nonnegative counts (MATLAB's gammaln rejects negative arguments).
v = -inf(size(k)); ok = k >= 0;
if lam <= 0, v(ok & k==0) = 0; return; end
v(ok) = k(ok)*log(lam) - lam - gammaln(k(ok)+1);
end
