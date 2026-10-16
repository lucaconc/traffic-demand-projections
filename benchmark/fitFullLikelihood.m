function [mu,info] = fitFullLikelihood(A,Y,mu0,varargin)
%FITFULLLIKELIHOOD Exact maximum-likelihood estimate for small designs.
%
%   [mu,info] = fitFullLikelihood(A,Y,mu0)
%
%   For Y = A X with independent X_j ~ Poisson(mu_j), the joint probability
%   table of Y on the box [0, max Y_i] is built exactly by successive
%   convolution: starting from a point mass at 0, path j shifts the table by
%   x*a_j with Poisson(mu_j) weights. Shifts only increase the counts, so the
%   truncation to the box is exact for every observed y. The derivatives use
%   the same table:
%     d p(y)/d mu_j = p(y-a_j) - p(y)
%     d2 p(y)/d mu_j d mu_k = p(y-a_j-a_k) - p(y-a_j) - p(y-a_k) + p(y)
%   and the log-likelihood is maximized over mu >= 0 by projected Newton steps
%   with backtracking (Levenberg damping when the Hessian is not negative
%   definite). Distinct observations are weighted by their counts.
%
%   Name/value: MaxCells (4e6: largest table), MaxIter (100), TolStep (1e-10),
%   EvaluateOnly (false: true returns the log-likelihood at mu0, without fitting).
%   info: loglik, iterations, cells, exitflag (1 converged, 0 limit).
p = inputParser; addParameter(p,'MaxCells',4e6); addParameter(p,'MaxIter',100);
addParameter(p,'TolStep',1e-10); addParameter(p,'EvaluateOnly',false); parse(p,varargin{:}); o = p.Results;
A = double(A~=0); [m,n] = size(A); Y = double(Y);
[U,~,g] = unique(Y,'rows'); w = accumarray(g,1);
dims = max(U,[],1) + 1;
if prod(dims) > o.MaxCells
    error('fitFullLikelihood:TooLarge','The probability table would have %d cells.',prod(dims));
end
if o.EvaluateOnly                                   % log-likelihood at mu0 itself (zeros allowed)
    mu = max(double(mu0(:)),0); L = loglik(mu);
    info = struct('loglik',L,'iterations',0,'cells',prod(dims),'exitflag',1); return
end
mu = max(double(mu0(:)), 1e-3*max(mean(mu0),1));
L = loglik(mu); flag = 0; lam = 0; it = 0;
for it = 1:o.MaxIter
    [~,gr,He] = evaluate(mu);
    step = newtonStep(gr,He,lam,mu);
    t = 1; improved = false;
    while t > 1e-8
        cand = max(mu + t*step, 0);
        Lc = loglik(cand);
        if Lc > L
            dmu = norm(cand-mu); mu = cand; L = Lc; improved = true; break
        end
        t = t/2;
    end
    if ~improved, lam = max(10*lam,1e-6*max(abs(diag(He)))); if lam > 1e12, break; end, continue; end
    lam = lam/10;
    if dmu <= o.TolStep*(1+norm(mu)), flag = 1; break; end
end
info = struct('loglik',L,'iterations',it,'cells',prod(dims),'exitflag',flag);

    function Lv = loglik(x)
        P = table_(x); pv = P(idx(U));
        if any(pv <= 0), Lv = -Inf; else, Lv = w.'*log(pv); end
    end
    function [Lv,gr,He] = evaluate(x)
        P = table_(x); pv = P(idx(U)); Lv = w.'*log(max(pv,realmin));
        D = zeros(size(U,1),n);
        for j = 1:n, D(:,j) = shifted(P,U,A(:,j).') - pv; end
        R = D./pv;
        gr = (w.'*R).';
        He = zeros(n);
        for j = 1:n
            for k = j:n
                d2 = shifted(P,U,(A(:,j)+A(:,k)).') - shifted(P,U,A(:,j).') - shifted(P,U,A(:,k).') + pv;
                He(j,k) = w.'*(d2./pv - R(:,j).*R(:,k)); He(k,j) = He(j,k);
            end
        end
    end
    function P = table_(x)
        % exact joint table on the box by successive convolution
        P = zeros([dims 1]); P(1) = 1;
        for j = 1:n
            a = A(:,j).'; xmax = min(floor((dims-1)./max(a,eps)));
            if ~isfinite(xmax), xmax = 0; end
            pk = exp((0:xmax)*log(max(x(j),realmin)) - x(j) - gammaln((0:xmax)+1));
            Q = zeros(size(P));
            for s = 0:xmax
                Q = Q + pk(s+1)*shiftTable(P,s*a,dims);
            end
            P = Q;
        end
    end
    function li = idx(V)
        c = num2cell(V+1,1); li = sub2ind([dims 1],c{:});
        if m == 1, li = V(:,1)+1; end
    end
    function v = shifted(P,V,d)
        % p(V - d), zero where a coordinate becomes negative
        Vs = V - d; v = zeros(size(V,1),1); ok = all(Vs >= 0,2);
        if any(ok), v(ok) = P(idx(Vs(ok,:))); end
    end
end

function S = shiftTable(P,d,dims)
% S(y) = P(y - d) on the box (zero where y - d is negative)
S = zeros(size(P));
src = arrayfun(@(k) 1:dims(k)-d(k), 1:numel(dims), 'UniformOutput', false);
dst = arrayfun(@(k) 1+d(k):dims(k), 1:numel(dims), 'UniformOutput', false);
if any(cellfun(@isempty,src)), return; end
S(dst{:}) = P(src{:});
end

function step = newtonStep(gr,He,lam,mu)
% Newton ascent step; Levenberg damping if -He is not positive definite.
n = numel(gr); M = -He + lam*eye(n);
[R,pd] = chol((M+M.')/2);
if pd ~= 0
    e = eig((M+M.')/2); M = M + (abs(min(e)) + 1e-8*max(1,max(abs(e))))*eye(n);
    R = chol((M+M.')/2);
end
step = R\(R.'\gr);
% do not push variables at zero further negative
step(mu <= 0 & step < 0) = 0;
end
