function [mu,info] = estimateMuEfficientGMM(A,Y,varargin)
%ESTIMATEMUEFFICIENTGMM Efficiently weighted projected scores (two-step GMM).
%
%   [mu,info] = estimateMuEfficientGMM(A,Y)
%   [mu,info] = estimateMuEfficientGMM(A,Y,'ProjectionMatrix',C,'Mu0',mu1)
%
%   For each projection Z_h = c_h'Y (compound Poisson with jump intensities
%   lambda_hk = sum_{j: b_hj=k} mu_j, b_h = A'c_h) the moment conditions are
%     rho_hk(Z) - 1,   rho_hk(z) = p_h(z-k)/p_h(z)   (one per jump size k; their
%                      equal-weight combination is the projected-likelihood score)
%     central-moment conditions of orders 1..r, r the identifying order (>= 2):
%       (Z-m) ; (Z-m)^2 - k2 ; (Z-m)^3 - k3 ; (Z-m)^4 - k4 - 3 k2^2 ;
%       (Z-m)^5 - k5 - 10 k2 k3 ,  with m = k1 and k_q = sum_k k^q lambda_hk
%       (the cumulants of a compound Poisson law). With the subset-sum
%       projections they carry the information of the joint cumulants through
%       order r, so the estimator is asymptotically at least as efficient as
%       cumulant matching (GLS) on the same data.
%   All have expectation zero at the true mu. Step one fits the projected
%   composite likelihood with the same projections (estimateMuOptimizedProjections);
%   step two minimizes gbar(mu)' W gbar(mu) over mu >= 0, W the inverse of the
%   covariance of the moment conditions (eigenvalues floored at Ridge times the
%   largest), by projected Levenberg-Marquardt with the analytic Jacobian
%   (d rho_k / d lambda_l = rho_{k+l} - rho_k rho_l), started from the step-one
%   estimate and from NNLS on the projections' cumulant equations (linear in mu);
%   the lower objective is kept, the weights are recomputed there and the fit is
%   repeated (WeightIterations).
%
%   Default projections: the nonzero rows of the cumulant system, i.e. the sums
%   of the sensors of every covered subset of size <= r, with r the smallest
%   order whose cumulant system identifies mu (at least 2), capped at
%   MaxProjections subsets in order of size.
%
%   Robust version ('Robust',true): every condition is redescending (a count the
%   model cannot explain contributes almost nothing, as in the contaminated
%   likelihood) and centred by its exact expectation under the model (computed
%   from the projection's probability table), so it has mean zero on clean data:
%     h_k(Z) - E_mu h_k(Z),  h_k(z) = R_k(z) - R_0(z),
%       R_s(z) = (1-eps) p(z-s) / ((1-eps) p(z) + eps q(z)), q the geometric
%       background law of the contaminated likelihood (mean: the projection's
%       sample mean, floored at one); h_k are the components of the contaminated
%       score, and both terms vanish for a spike;
%     R_0(Z) psi_c(D) - E_mu[.]  and  R_0(Z) min(D^2,c^2) - E_mu[.],
%       D = (Z - mean)/sd under the model, psi_c the Huber function (clip at +-c);
%       the model weight R_0 (close to 0 for a spike) makes them redescending.
%   The pilot is the projected CL with the contaminated likelihood (Epsilon), and
%   the weight matrix uses only intervals whose model weight R_0 is at least
%   WeightThreshold in every projection (all intervals if too few remain).
%   Derivatives are analytic: d R_s/d lambda_l = R_{s+l} - R_s - R_s (R_l - R_0),
%   the chain rule through D for the clipped conditions, and the shift identity
%   sum_z p(z-l) h(z) = sum_z p(z) h(z+l) for the centring terms.
%
%   info: C, nMoments, firstStep (mu), firstStepSeconds, gmmSeconds, seconds,
%   iterations, objective, exitflag (1 converged, 0 iteration limit), robust,
%   cleanShare (share of intervals used for the weights, robust version).
%
%   Name/value: ProjectionMatrix ([]), Mu0 ([]: step one), MaxProjections (300),
%   Shrinkage ('auto'): weight matrix. 'none': inverse sample covariance of the
%   conditions (efficient weights). 'lw': covariance shrunk toward its diagonal
%   (Ledoit-Wolf), robust when it is near-singular (many conditions, many zero
%   rates). 'auto': efficient weights, and a refit with 'lw' if Hansen's J =
%   N*objective exceeds df+10*sqrt(2*df)+50 (df = conditions - parameters);
%   info.weighting and info.Jstat report the outcome. Steps are projected
%   Levenberg-Marquardt with an active set: rates at zero whose gradient points
%   outward stay fixed.
%   Ridge (1e-10), MaxIter (200), TolStep (1e-9), Display ('off'), Robust (false),
%   Epsilon (0.01), HuberC (5), WeightThreshold (0.5), WeightIterations (2: the
%   robust weights are recomputed once at the first robust estimate and the fit
%   repeated), NumericJacobian (false: analytic derivatives), MomentOrder ([]:
%   the identifying order, at least 2; plain version), CumulantStart (true: also
%   start from NNLS on the cumulant equations; false: only from the pilot or
%   Mu0, about half the time), PlainPilot ('pcl': the
%   projected CL on its own small design, fast and close to the optimum;
%   'pclFull': on all GMM projections, slow; 'cumulants': NNLS on the
%   projections' cumulant equations only), RobustPilotDesign ('full': the
%   contaminated projected CL on all GMM projections, most accurate; 'small':
%   on its own small design, faster), Pilot ('cleanGMM':
%   robust PCL, then the plain GMM on the intervals it considers clean;
%   'robustPCL': the robust PCL alone).

ip = inputParser;
addParameter(ip,'ProjectionMatrix',[]);
addParameter(ip,'Mu0',[]);
addParameter(ip,'MaxProjections',300);
addParameter(ip,'Ridge',[]);   % truncation tolerance; [] = 1e-5/N (vanishes as N grows)
addParameter(ip,'MaxIter',200);
addParameter(ip,'Shrinkage','auto',@(x)any(strcmpi(x,{'auto','lw','none'})));
addParameter(ip,'RemoveRedundant',true,@(x)islogical(x)&&isscalar(x));   % drop exactly dependent conditions
addParameter(ip,'Covariance','auto',@(x)any(strcmpi(x,{'auto','sample','model'})));   % source of S
addParameter(ip,'SimulationSize',2000); addParameter(ip,'SimulationSeed',1);
addParameter(ip,'TolStep',1e-9);
addParameter(ip,'Display','off');
addParameter(ip,'Robust',false);
addParameter(ip,'Epsilon',0.01);
addParameter(ip,'HuberC',5);
addParameter(ip,'WeightThreshold',0.5);
addParameter(ip,'WeightIterations',2);
addParameter(ip,'DualPilot',true);   % also run the two-step fit from the cumulant pilot; keep the lower CU criterion
addParameter(ip,'MomentOrder',[]);
addParameter(ip,'CumulantStart',true);
addParameter(ip,'RobustPilotDesign','full',@(x)any(strcmpi(x,{'full','small'})));
addParameter(ip,'PlainPilot','pcl',@(x)any(strcmpi(x,{'pcl','pclFull','cumulants'})));
addParameter(ip,'NumericJacobian',false);
addParameter(ip,'Pilot','cleanGMM',@(x)any(strcmpi(x,{'cleanGMM','robustPCL'})));
parse(ip,varargin{:}); o = ip.Results;
t0 = tic;
A = double(A~=0); [m,n] = size(A); Y = double(Y); N = size(Y,1);
% Sensors that never counted a vehicle force every class crossing them to rate
% zero (their mean count is zero, and rates are nonnegative); every consistent
% estimator gives zero there. These classes are fixed at zero and the GMM is
% fitted on the other sensors and classes, which stay nonzero and distinct.
% Without this, silent sensors duplicate projections and make the condition
% covariance deficient, so that the identification check would refuse the fit.
silent = all(Y == 0,1);
if any(silent) && ~all(silent)
    zc = any(A(silent,:),1);
    if all(zc)
        mu = zeros(n,1); info = struct('status','every class crosses a silent sensor','fixedZeroClasses',n,'silentSensors',nnz(silent)); return
    end
    keepS = ~silent; keepC = ~zc; args = varargin;
    args = restrictArg(args,'ProjectionMatrix',@(Cm) uniqueNonzeroRows(double(Cm(:,keepS))));
    args = restrictArg(args,'Mu0',@(v) v(keepC));
    [muR,info] = estimateMuEfficientGMM(A(keepS,keepC),Y(:,keepS),args{:});
    mu = zeros(n,1); mu(keepC) = muR(:);
    info.fixedZeroClasses = nnz(zc); info.silentSensors = nnz(silent);
    return
end
[Cdef,rIdent] = defaultProjections(A,o.MaxProjections);
if isempty(o.ProjectionMatrix), C = Cdef; else, C = double(o.ProjectionMatrix); end
% Default: the identifying order (at least 2); an explicit MomentOrder may be 0
% (score conditions only) to 5.
if isempty(o.MomentOrder), r = min(max(rIdent,2),5); else, r = min(max(round(o.MomentOrder),0),5); end
H = size(C,1); B = C*A;                                % H x n jump sizes
% ---- step one -------------------------------------------------------------
tf = tic;
if isempty(o.Mu0)
    if o.Robust
        if strcmpi(o.RobustPilotDesign,'small')     % faster, slightly less accurate
            mu1 = estimateMuOptimizedProjections(A,Y,'Display','off','Contamination',o.Epsilon);
        else
            mu1 = estimateMuOptimizedProjections(A,Y,'ProjectionMatrix',C,'Display','off','Contamination',o.Epsilon);
        end
        if strcmpi(o.Pilot,'cleanGMM')
            ok = cleanIntervals(A,Y,C,mu1,o);
            if nnz(ok) >= max(50,0.5*N)
                mu1 = estimateMuEfficientGMM(A,Y(ok,:),'ProjectionMatrix',C,'Mu0',mu1, ...
                    'MomentOrder',r,'CumulantStart',o.CumulantStart);
            end
        end
    elseif strcmpi(o.PlainPilot,'pcl')
        mu1 = estimateMuOptimizedProjections(A,Y,'Display','off');     % own small design
    elseif strcmpi(o.PlainPilot,'pclFull')
        mu1 = estimateMuOptimizedProjections(A,Y,'ProjectionMatrix',C,'Display','off');
    else
        mu1 = [];                                       % set below: cumulant-equation start
    end
else
    mu1 = double(o.Mu0(:));
end
firstSeconds = toc(tf);
% ---- per-projection data ----------------------------------------------------
P = buildProjections(Y,C,B);
nMom = sum(arrayfun(@(p) numel(p.jumps)+2, P));
if ~o.Robust, nMom = sum(arrayfun(@(p) numel(p.jumps)+r, P)); end
if o.Robust
    [mu,info] = robustStepTwo(P,mu1,nMom,n,N,o);
    info.C = C; info.firstStep = mu1; info.firstStepSeconds = firstSeconds;
    info.seconds = toc(t0); info.gmmSeconds = info.seconds - firstSeconds; info.robust = true;
    return
end
% ---- starts: NNLS on the projections' cumulant equations (and step one) -------
if r == 0 && isempty(mu1), mu1 = estimateMuOptimizedProjections(A,Y,'Display','off'); end   % no cumulant start at order 0
tc = tic; if r >= 1, mu2 = cumulantStart(P,Y,C,r,n); else, mu2 = max(mu1(:),0); end
if isempty(mu1), mu1 = mu2; firstSeconds = firstSeconds + toc(tc); end
% Efficient weights first; if Hansen's J statistic N*Q lies far beyond its
% chi-square range (df = conditions - parameters), the weight matrix has failed
% (typically near-singular, e.g. many zero rates): refit with weights shrunk
% toward the diagonal (Ledoit-Wolf). 'none' / 'lw' force one weighting.
% Exactly dependent conditions are removed: by Panjer's identity
% sum_k k*lambda_k*(rho_k(z)-1) = z - m_c, the first-moment condition of every
% projection is a combination of its ratio conditions, and for a projection whose
% sensors all have singleton projections, z - m_c = sum_i c_i*m_i*(rho_i,1 - 1)
% makes its largest-jump ratio condition a combination of the others.
if o.RemoveRedundant, o.Keep = nonRedundantConditions(P,C,r,max(mu1(:),0)); else, o.Keep = true(1,nMom); end
% Source of the covariance S of the retained conditions: the sample if it can
% estimate it (N-1 >= M), otherwise the model, from max(SimulationSize,20M)
% intervals simulated at the current estimate (rank independent of N; valid
% under the model).
M = nnz(o.Keep); covSource = lower(o.Covariance);
if strcmp(covSource,'auto'), if N-1 >= M, covSource = 'sample'; else, covSource = 'model'; end, end
o.CovSource = covSource;
o.SimCtx = struct('A',A,'C',C,'B',B,'Ms',max(o.SimulationSize,20*M),'seed',o.SimulationSeed);
if strcmpi(o.Shrinkage,'lw'), mode = 'lw'; else, mode = 'none'; end
[best,it,flag,dEff,identified] = weightedFit(P,mu1,mu2,N,nMom,n,r,o,mode);
% The first weights are computed at the pilot, and a poor pilot can lead the
% fit into a poor basin. With DualPilot, the two-step fit is repeated from the
% cumulant pilot, and the two results are compared by the continuously updated
% criterion gbar(mu)'*W(mu)*gbar(mu), with the weights at each candidate itself,
% since two-step objectives computed with different weights are not comparable.
pilotUsed = 'pcl';
if o.DualPilot && identified && norm(mu2(:)-mu1(:)) > 1e-8*(1+norm(mu1(:)))
    [best2,it2,flag2,dEff2,id2] = weightedFit(P,mu2,mu1,N,nMom,n,r,o,mode);
    it = it + it2;
    if id2 && cuCriterion(P,best2.mu,N,nMom,n,r,o) < cuCriterion(P,best.mu,N,nMom,n,r,o)
        best = best2; flag = flag2; dEff = dEff2; pilotUsed = 'cumulants';
    end
end
J = NaN; df = dEff - n; weighting = mode; status = 'ok';
if ~identified
    % The retained conditions do not identify mu (rank of G'WG below n):
    % the GMM is not defined, and the pilot (muC by default) is returned.
    best = struct('mu',max(mu1(:),0),'Q',NaN); status = 'not identified by the retained conditions: pilot returned';
elseif df > 0
    % Hansen's J, approx. chi-square with d - n degrees of freedom (correct model, interior mu)
    Jmax = df + 10*sqrt(2*df) + 50; J = N*best.Q;
    if strcmpi(o.Shrinkage,'auto') && ~(J <= Jmax)
        [bestL,itL,flagL,~,idL] = weightedFit(P,mu1,mu2,N,nMom,n,r,o,'lw');
        if idL, it = it + itL; best = bestL; flag = flagL; J = N*best.Q; weighting = 'lw (fallback)'; end
    end
end
mu = best.mu; Q = best.Q;
info = struct('C',C,'nMoments',nMom,'momentOrder',r,'firstStep',mu1,'cumulantStart',mu2, ...
    'firstStepSeconds',firstSeconds,'gmmSeconds',toc(t0)-firstSeconds,'seconds',toc(t0), ...
    'iterations',it,'objective',Q,'exitflag',flag,'robust',false,'cleanShare',1, ...
    'weighting',weighting,'Jstat',J,'Jdf',df,'conditionsKept',nnz(o.Keep),'effectiveDim',dEff, ...
    'covariance',covSource,'status',status,'fixedZeroClasses',0,'silentSensors',0,'pilotUsed',pilotUsed);
end

function [best,it,flag,dEff,identified] = weightedFit(P,mu1,mu2,N,nMom,n,r,o,mode)
% Two-step fit with the given weighting ('none': efficient weights; 'lw': shrunk).
% identified is false if, at the start, G'WG has rank below n (the retained
% conditions, weighted by W, do not identify mu); the fit is then not run.
oW = o; oW.Shrinkage = mode;
best = struct('mu',max(mu1(:),0),'Q',Inf); it = 0; flag = 0; dEff = NaN; identified = true;
for wit = 1:max(1,o.WeightIterations)
    if isfield(o,'CovSource') && strcmp(o.CovSource,'model')
        G1 = simulatedConditions(o.SimCtx,best.mu,nMom,r);
    else
        G1 = momentsPerObservation(P,best.mu,N,nMom,r);
    end
    if isfield(o,'Keep'), G1 = G1(:,o.Keep); end
    [W,dEff] = conditionWeights(G1,oW);
    if wit == 1
        [~,Jg0] = meanMoments(P,best.mu,nMom,n,r); [~,Jg0] = keepRows(o,zeros(size(Jg0,1),1),Jg0);
        Hs0 = Jg0.'*W*Jg0; Hs0 = (Hs0+Hs0.')/2; ev = eig(Hs0);
        if dEff < n || ~(min(ev) > 1e-10*max(abs(ev))), identified = false; return; end
    end
    starts = {best.mu};
    if wit == 1 && o.CumulantStart && norm(mu2-best.mu) > 1e-12*(1+norm(mu2)), starts{end+1} = mu2; end
    best.Q = Inf;
    for st = 1:numel(starts)
        [muS,QS,itS,flS] = lmFit(P,starts{st},W,nMom,n,r,o);
        it = it + itS;
        if QS < best.Q, best = struct('mu',muS,'Q',QS); flag = flS; end
    end
end
end

% -------------------------------------------------------------------------
function [C,r] = defaultProjections(A,maxProj)
% Sums over every covered sensor subset of size <= r (r: identifying order, >= 2).
[m,n] = size(A); C = zeros(0,m); r = 2;
for ord = 1:m
    rows = coveredSubsets(A,ord);
    if rank(subsetRows(A,rows)) == n, r = max(2,ord); break; end
end
rows = coveredSubsets(A,r);
for k = 1:numel(rows)
    c = zeros(1,m); c(rows{k}) = 1; C(end+1,:) = c; %#ok<AGROW>
    if size(C,1) >= maxProj, break; end
end
end

function rows = coveredSubsets(A,order)
m = size(A,1); rows = {};
for s = 1:min(order,m)
    cb = nchoosek(1:m,s);
    for k = 1:size(cb,1)
        if any(all(A(cb(k,:),:),1)), rows{end+1} = cb(k,:); end %#ok<AGROW>
    end
end
end

function R = subsetRows(A,rows)
R = zeros(numel(rows),size(A,2));
for k = 1:numel(rows), R(k,:) = all(A(rows{k},:),1); end
end

function lp = compoundLogPmf(lambda,jumps,zmax)
% log P(Z=z), z=0..zmax, Z compound Poisson with intensities lambda at jumps
% (Panjer recursion in the log domain). If the compiled compoundLogPmfMex is on the
% path (buildPaperMex), it is used; it computes the same terms in the same order.
% Otherwise, under MATLAB, a scalar loop is used for
% up to four jump sizes: the JIT compiles it, which is much faster than the
% small-vector version; both compute the same terms in the same order.
persistent useScalar useMex
if isempty(useScalar), useScalar = exist('OCTAVE_VERSION','builtin') == 0; end
if isempty(useMex), useMex = exist('compoundLogPmfMex','file') == 3; end
if useMex                                   % compiled C version (buildPaperMex), same terms and order
    lp = compoundLogPmfMex(double(lambda(:)),double(jumps(:)),double(zmax));
    return
end
lp = -inf(zmax+1,1); tot = sum(lambda); lp(1) = -tot;
K = max(jumps); coef = zeros(1,K); coef(jumps) = jumps(:).'.*lambda(:).';
if useScalar && K <= 4
    lc = -inf(1,K); pos = coef > 0; lc(pos) = log(coef(pos));
    for z = 1:zmax
        m = -Inf; any_ = false;
        for k = 1:min(K,z)
            if pos(k)
                t = lc(k) + lp(z-k+1);
                if t > -Inf
                    if ~any_ || t > m, m = t; end
                    any_ = true;
                end
            end
        end
        if ~any_, continue; end
        acc = 0;
        for k = 1:min(K,z)
            if pos(k)
                t = lc(k) + lp(z-k+1);
                if t > -Inf, acc = acc + exp(t-m); end
            end
        end
        lp(z+1) = m + log(acc) - log(z);
    end
    return
end
for z = 1:zmax
    k = 1:min(z,numel(coef)); k = k(coef(k)>0);
    if isempty(k), continue; end
    t = log(coef(k)) + lp(z-k+1).';
    mx = max(t); if ~isfinite(mx), continue; end
    lp(z+1) = mx + log(sum(exp(t-mx))) - log(z);
end
end

function [rho,lp] = ratios(p,lambda)
% rho(:,s) = p(u - s)/p(u) for s = 1..2*maxJump (needed for the Jacobian)
K = max(p.jumps); lp = compoundLogPmf(lambda,p.jumps,p.zmax);
u = p.u; rho = zeros(numel(u),2*K);
for s = 1:2*K
    v = -inf(numel(u),1); ok = u-s >= 0; v(ok) = lp(u(ok)-s+1);
    rho(:,s) = exp(v - lp(u+1));
end
end

function G = momentsPerObservation(P,mu,N,nMom,r)
% N x nMom matrix of moment conditions at mu (for the weight matrix).
if nargin < 5, r = 2; end
G = zeros(N,nMom); col = 0;
for h = 1:numel(P)
    p = P(h); lambda = (p.J.'*mu).'; rho = ratios(p,lambda);
    for k = 1:numel(p.jumps)
        col = col+1; G(:,col) = rho(p.idx,p.jumps(k)) - 1;
    end
    kap = lambda(:).' * (p.jumps(:).^(1:5)); cm = centralModel(kap);     % cumulants of orders 1..5
    dev = p.z - kap(1);
    for q = 1:r, col = col+1; G(:,col) = dev.^q - cm(q); end
end
end

function [g,Jg] = meanMoments(P,mu,nMom,n,r)
% Mean moment conditions and their Jacobian with respect to mu.
if nargin < 5, r = 2; end
g = zeros(nMom,1); Jg = zeros(nMom,n); row = 0;
for h = 1:numel(P)
    p = P(h); lambda = (p.J.'*mu).'; w = p.w/sum(p.w);
    rho = ratios(p,lambda); K = numel(p.jumps);
    for k = 1:K
        row = row+1; jk = p.jumps(k);
        g(row) = w.'*rho(:,jk) - 1;
        if nargout > 1
            dl = zeros(1,K);                            % d/d lambda_l
            for l = 1:K
                jl = p.jumps(l);
                dl(l) = w.'*(rho(:,jk+jl) - rho(:,jk).*rho(:,jl));
            end
            Jg(row,:) = (p.J*dl.').';
        end
    end
    kap = lambda(:).' * (p.jumps(:).^(1:5)); cm = centralModel(kap);     % cumulants of orders 1..5
    dev = p.u - kap(1); avg = w.' * (dev(:).^(0:5));   % avg(q+1) = mean of dev^q
    l = p.jumps(:);
    for q = 1:r
        g(row+q) = avg(q+1) - cm(q);
        switch q                                   % d cm(q) / d lambda_l
            case 1, dcm = zeros(size(l));
            case 2, dcm = l.^2;
            case 3, dcm = l.^3;
            case 4, dcm = l.^4 + 6*kap(2)*l.^2;
            case 5, dcm = l.^5 + 10*(l.^2*kap(3) + l.^3*kap(2));
        end
        dl = -q*l*avg(q) - dcm;                    % d/dlambda_l of mean dev^q: -q l mean dev^(q-1)
        Jg(row+q,:) = (p.J*dl).';
    end
    row = row+r;
end
end

% ===== robust version =========================================================
function [mu,info] = robustStepTwo(P,mu1,nMom,n,N,o)
% Weights from the intervals the pilot considers clean, then projected
% Levenberg-Marquardt on the bounded, model-centred conditions.
for h = 1:numel(P)
    zbar = max(mean(P(h).z),1); th = zbar/(1+zbar);          % geometric background law
    P(h).logq1 = log(1-th); P(h).logth = log(th);
end
mu = max(mu1(:),0); itTotal = 0;
for wit = 1:max(1,o.WeightIterations)
[G1,clean] = robustPerObservation(P,mu,N,nMom,o);
keep = clean;
if nnz(keep) < max(50,2*nMom), keep = true(N,1); end
oW = o; if strcmpi(oW.Shrinkage,'auto'), oW.Shrinkage = 'none'; end
if isempty(oW.Ridge), oW.Ridge = 1e-10; end          % robust variant (not in the paper): previous tolerance
W = conditionWeights(G1(keep,:),oW);
lam = 1e-3; flag = 0;
[g,Jg] = robustMeanMoments(P,mu,nMom,n,o);
Q = g.'*W*g; it = 0;
for it = 1:o.MaxIter
    grad = Jg.'*W*g; Hgn = Jg.'*W*Jg; improved = false;
    while lam < 1e12
        step = boundedStep(Hgn,grad,lam,mu);
        cand = max(mu+step,0);
        gc = robustMeanMoments(P,cand,nMom,n,o);        % conditions only
        Qc = gc.'*W*gc;
        if Qc < Q
            [g,Jg] = robustMeanMoments(P,cand,nMom,n,o); % Jacobian only when accepted
            dmu = norm(cand-mu); mu = cand; Q = Qc;
            lam = max(lam/10,1e-12); improved = true; break
        end
        lam = lam*10;
    end
    if ~improved || dmu <= o.TolStep*(1+norm(mu)), flag = 1; break; end
end
itTotal = itTotal + it;
end
info = struct('nMoments',nMom,'iterations',itTotal,'objective',Q,'exitflag',flag,'cleanShare',mean(keep));
end

function T = robustTable(p,lambda,o)
% Model probability table and the bounded conditions on the grid 0..Zt.
mean_ = sum(p.jumps.*lambda); sd = sqrt(max(sum(p.jumps.^2.*lambda),eps));
Zt = max(p.zmax, ceil(mean_ + 12*sd + 10));
K = max(p.jumps);
lp = compoundLogPmf(lambda,p.jumps,Zt); z = (0:Zt).';
lq = p.logq1 + z*p.logth;
a = log(1-o.Epsilon) + lp; b = log(o.Epsilon) + lq;
mx = max(a,b); logD = mx + log(exp(a-mx)+exp(b-mx));
R = zeros(Zt+1,2*K+1);                                  % R(:,s+1) = (1-eps) p(z-s)/D(z)
for sft = 0:2*K
    v = -inf(Zt+1,1); ok = z-sft >= 0; v(ok) = lp(z(ok)-sft+1);
    R(:,sft+1) = exp(log(1-o.Epsilon) + v - logD);
end
D_ = (z - mean_)/sd; c = o.HuberC;
R0 = R(:,1);
psi = max(-c,min(c,D_)); inside = abs(D_) < c; sq0 = min(D_.^2,c^2);
hub = R0.*psi; sq = R0.*sq0;                                % redescending through R_0
H = R(:,2:K+1) - R0;                                        % contaminated score components
pz = exp(lp);
T = struct('R',R,'H',H,'hub',hub,'sq',sq,'EH',pz.'*H,'Ehub',pz.'*hub,'Esq',pz.'*sq, ...
    'pz',pz,'psi',psi,'inside',inside,'sq0',sq0,'D',D_,'sd',sd,'Zt',Zt);
end

function gh = robustConditionsMean(p,lambda,o)
% Mean over the observed counts of the model-centred bounded conditions.
T = robustTable(p,lambda,o); w = p.w/sum(p.w); u = p.u+1;
K = numel(p.jumps); gh = zeros(K+2,1);
for k = 1:K, jk = p.jumps(k); gh(k) = w.'*T.H(u,jk) - T.EH(jk); end
gh(K+1) = w.'*T.hub(u) - T.Ehub;
gh(K+2) = w.'*T.sq(u) - T.Esq;
end

function [g,Jg] = robustMeanMoments(P,mu,nMom,n,o)
% Conditions and their Jacobian with respect to mu (analytic, or central
% differences in each projection's jump intensities if NumericJacobian).
g = zeros(nMom,1); Jg = zeros(nMom,n); row = 0;
for h = 1:numel(P)
    p = P(h); lambda = (p.J.'*mu).'; K = numel(p.jumps);
    if nargout < 2
        gh = robustConditionsMean(p,lambda,o); dl = zeros(K+2,K);
    elseif o.NumericJacobian
        gh = robustConditionsMean(p,lambda,o);
        dl = zeros(K+2,K);
        for l = 1:K
            st = 1e-5*max(lambda(l),1); lp_ = lambda; lm = lambda; lp_(l) = lp_(l)+st;
            lm(l) = max(lm(l)-st,0); dl(:,l) = (robustConditionsMean(p,lp_,o) - robustConditionsMean(p,lm,o))/(lp_(l)-lm(l));
        end
    else
        [gh,dl] = robustConditionsAnalytic(p,lambda,o);
    end
    g(row+(1:K+2)) = gh; Jg(row+(1:K+2),:) = dl*p.J.'; row = row+K+2;
end
end

function [gh,dl] = robustConditionsAnalytic(p,lambda,o)
% Mean conditions of one projection and their derivatives in its jump
% intensities. Grid functions f(z) on 0..Zt; observed mean uses weights w on p.u.
T = robustTable(p,lambda,o); w = p.w/sum(p.w); u = p.u+1;
K = numel(p.jumps); jumps = p.jumps; R = T.R; pz = T.pz; Zt = T.Zt;
gh = zeros(K+2,1); dl = zeros(K+2,K);
% dR_s/dlambda for jump size jl: R_{s+jl} - R_s - R_s (R_jl - R_0)
dR = @(sft,jl) R(:,sft+jl+1) - R(:,sft+1) - R(:,sft+1).*(R(:,jl+1) - R(:,1));
shiftUp = @(f,jl) [f(jl+1:end); repmat(f(end),jl,1)];     % f(z+jl) on the grid (edge: last value)
for l = 1:K
    jl = jumps(l);
    % d/dlambda of the expectation E[f] = sum_z pz f: sum_z pz (f(z+jl) - f(z)) + sum_z pz df
    dR0 = dR(0,jl);
    dD = -jl/T.sd - T.D*jl^2/(2*T.sd^2);
    dhub = dR0.*T.psi + R(:,1).*T.inside.*dD;
    dsq = dR0.*T.sq0 + R(:,1).*(T.inside.*(2*T.D.*dD));
    for k = 1:K
        f = T.H(:,jumps(k)); df = dR(jumps(k),jl) - dR0;     % columns of H are indexed by jump size
        dl(k,l) = w.'*df(u) - (pz.'*(shiftUp(f,jl) - f) + pz.'*df);
    end
    dl(K+1,l) = w.'*dhub(u) - (pz.'*(shiftUp(T.hub,jl) - T.hub) + pz.'*dhub);
    dl(K+2,l) = w.'*dsq(u) - (pz.'*(shiftUp(T.sq,jl) - T.sq) + pz.'*dsq);
end
for k = 1:K, gh(k) = w.'*T.H(u,jumps(k)) - T.EH(jumps(k)); end
gh(K+1) = w.'*T.hub(u) - T.Ehub;
gh(K+2) = w.'*T.sq(u) - T.Esq;
end

function [G,clean] = robustPerObservation(P,mu,N,nMom,o)
% Per-interval conditions (for the weights) and the clean-interval flags.
G = zeros(N,nMom); clean = true(N,1); col = 0;
for h = 1:numel(P)
    p = P(h); lambda = (p.J.'*mu).'; T = robustTable(p,lambda,o); zi = p.z+1; K = numel(p.jumps);
    for k = 1:K, col = col+1; G(:,col) = T.H(zi,p.jumps(k)) - T.EH(p.jumps(k)); end
    G(:,col+1) = T.hub(zi) - T.Ehub; G(:,col+2) = T.sq(zi) - T.Esq; col = col+2;
    clean = clean & T.R(zi,1) >= o.WeightThreshold;
end
end

function ok = cleanIntervals(A,Y,C,mu,o)
% Intervals whose projected counts all have model weight R_0 >= WeightThreshold
% under the contaminated law at mu.
ok = true(size(Y,1),1); B = C*A;
for h = 1:size(C,1)
    z = Y*C(h,:).'; jumps = unique(B(h,B(h,:)>0));
    lambda = arrayfun(@(k) sum(mu(B(h,:)==k)),jumps);
    lp = compoundLogPmf(lambda,jumps,max(z));
    zb = max(mean(z),1); th = zb/(1+zb);
    a = log(1-o.Epsilon) + lp(z+1); b = log(o.Epsilon) + log(1-th) + z*log(th);
    mx = max(a,b); R0 = exp(a - (mx + log(exp(a-mx)+exp(b-mx))));
    ok = ok & R0 >= o.WeightThreshold;
end
end

% ===== plain version: higher-order conditions, starts and fit ===================
function [mu,Q,it,flag] = lmFit(P,mu0,W,nMom,n,r,o)
% Projected Levenberg-Marquardt on gbar' W gbar from mu0.
mu = max(mu0(:),0); lam = 1e-3; flag = 0;
[g,Jg] = meanMoments(P,mu,nMom,n,r); [g,Jg] = keepRows(o,g,Jg);
Q = g.'*W*g; it = 0;
for it = 1:o.MaxIter
    grad = Jg.'*W*g; Hgn = Jg.'*W*Jg; improved = false;
    while lam < 1e12
        step = boundedStep(Hgn,grad,lam,mu);
        cand = max(mu+step,0);
        gc = meanMoments(P,cand,nMom,n,r); gc = keepRows(o,gc);   % conditions only
        Qc = gc.'*W*gc;
        if Qc < Q
            [g,Jg] = meanMoments(P,cand,nMom,n,r); [g,Jg] = keepRows(o,g,Jg);   % Jacobian only when accepted
            dmu = norm(cand-mu); mu = cand; Q = Qc;
            lam = max(lam/10,1e-12); improved = true; break
        end
        lam = lam*10;
    end
    if ~improved || dmu <= o.TolStep*(1+norm(mu)), flag = 1; break; end
    pg = mu - max(mu - grad,0);                      % projected gradient
    if norm(pg,inf) <= 1e-10*max(1,norm(grad,inf)) && dmu <= 1e-6*(1+norm(mu)), flag = 1; break; end
end
end

function mu = cumulantStart(P,Y,C,r,n) %#ok<INUSL>
% NNLS on the projections' cumulant equations k_q(Z_h) = sum_j b_hj^q mu_j.
rows = zeros(0,n); rhs = zeros(0,1);
for h = 1:numel(P)
    z = P(h).z; m = mean(z); c = z - m;
    m2 = mean(c.^2); m3 = mean(c.^3); m4 = mean(c.^4); m5 = mean(c.^5);
    kq = [m, m2, m3, m4-3*m2^2, m5-10*m3*m2];
    b = (P(h).J*P(h).jumps(:)).';                     % jump size of each class in projection h
    for q = 1:r
        rows(end+1,:) = b.^q; rhs(end+1,1) = kq(q); %#ok<AGROW>
    end
end
sc = 1./max(sqrt(sum(rows.^2,2)),eps);                 % equilibrate the orders
mu = lsqnonneg(rows.*sc,rhs.*sc,optimset('Display','off'));
% Strictly positive start: a zero rate can make observed counts impossible under
% the pilot, which would give infinite conditions and no weight matrix.
mu = max(mu, 0.02*max(mean(mu),eps));
end

function cm = centralModel(kap)
% Model central moments of orders 1..5 (targets of mean (Z-m)^q) from cumulants.
cm = [0, kap(2), kap(3), kap(4)+3*kap(2)^2, kap(5)+10*kap(2)*kap(3)];
end

function step = dampedStep(Hgn,grad,lam)
% Levenberg-Marquardt step -(H + lam*D)\grad, with the damping scale D floored
% relative to the largest diagonal entry (weakly determined directions get real
% regularization), solved by Cholesky; if the damped matrix is still not
% positive definite, a minimum-norm least-squares step is used instead. The
% accept/reject test of the caller guarantees descent in any case.
d = diag(Hgn); d = max(d, 1e-8*max([d; eps]));
M = Hgn + lam*diag(d); M = (M+M.')/2;
[R,p] = chol(M);
if p == 0
    step = -(R\(R.'\grad));
else
    w = warning('off','MATLAB:nearlySingularMatrix'); w2 = warning('off','MATLAB:singularMatrix');
    c = onCleanup(@() warning([w w2]));
    step = -(pinv(M)*grad);
end
end

function step = boundedStep(Hgn,grad,lam,mu)
% Projected Levenberg-Marquardt step for mu >= 0 (Bertsekas-type active set):
% rates at zero whose gradient points outward (descent would make them
% negative) are held fixed, and the damped step is computed for the others.
% Clipping a full step instead wastes it on the bound and stalls when many
% rates are zero.
act = (mu <= 0) & (grad > 0);
step = zeros(size(mu));
free = ~act;
if any(free), step(free) = dampedStep(Hgn(free,free),grad(free),lam); end
end

function [W,dEff] = conditionWeights(G,o)
% Weight matrix from the conditions G (N x M), on the correlation scale.
% 'none': inverse of the sample covariance on its numerically nonzero spectrum
%   (eigenvalues of the correlation matrix below tau times the largest are
%   dropped, not floored, so near-dependent directions get no weight), with
%   tau = o.Ridge, or tau_N = 1e-5/N by default, which vanishes as N grows, so that
%   asymptotically no direction of a nonsingular covariance is dropped.
% 'lw': the same after Ledoit-Wolf shrinkage toward the diagonal, with intensity
%   delta = sum of estimated variances of the off-diagonal sample covariances /
%   sum of their squares; for a fixed number of conditions delta = O(1/N), so the
%   shrinkage vanishes and the weights remain asymptotically efficient.
% Conditions without variation get no weight. dEff: number of retained
% directions (effective dimension of the conditions).
N = size(G,1); M = size(G,2); W = zeros(M); dEff = 0;
Gc = G - mean(G,1);
S = (Gc.'*Gc)/max(N-1,1); S = (S+S.')/2;
v = diag(S); live = isfinite(v) & v > 1e-12*max([v(isfinite(v)); realmin]);
if ~any(live), return; end
Sl = S(live,live);
if strcmpi(o.Shrinkage,'lw')
    Gl = Gc(:,live); Q2 = Gl.^2;
    varS = max((Q2.'*Q2)/N - ((Gl.'*Gl)/N).^2, 0)/N;      % estimated variance of each s_ij
    off = ~eye(size(Sl)); den = sum(Sl(off).^2);
    if den > 0, delta = min(1, max(0, sum(varS(off))/den)); else, delta = 1; end
    Sl = (1-delta)*Sl + delta*diag(diag(Sl));
end
sd = sqrt(diag(Sl)); Rm = Sl./(sd*sd.');
tol = o.Ridge; if isempty(tol), tol = max(1e-12, 1e-5/N); end   % tau_N = 1e-5/N -> 0
[V,D] = eig((Rm+Rm.')/2); d = diag(D); keep = d > tol*max(d);
dEff = nnz(keep);
W(live,live) = (V(:,keep)*diag(1./d(keep))*V(:,keep).')./(sd*sd.');
W = (W+W.')/2;
end

function [g,Jg] = keepRows(o,g,Jg)
% Apply the mask of non-redundant conditions.
if ~isfield(o,'Keep'), return; end
g = g(o.Keep);
if nargin > 2, Jg = Jg(o.Keep,:); end
end

function keep = nonRedundantConditions(P,C,r,muPilot)
% Mask over the stacked conditions [rho_k - 1 (k = jumps), central moments q = 1..r]
% of every projection. Dropped, because exact linear combinations of others:
%  - every first-moment condition (Panjer: sum_k k*lambda_k*(rho_k - 1) = z - m_c);
%  - for each non-singleton projection whose sensors all have singleton
%    projections (z - m_c = sum_i c_i*m_i*(rho_i,1 - 1)), the ratio condition
%    with the largest coefficient k*lambda_k at the pilot; the identity makes it
%    redundant only where lambda_k > 0, so nothing is dropped if every jump
%    class of the projection is empty at the pilot;
%  - the order-q moment condition (q >= 2) of a binary projection S with |S| > q
%    when every subset T of S with |T| <= q is a binary projection of the design
%    (polarization: (sum_S u_i)^q is a fixed combination of (sum_T u_i)^q).
C = double(C); isBinary = all(C == 0 | C == 1, 2);
isSingle = isBinary & sum(C,2) == 1;
singleSensors = find(any(C(isSingle,:),1));
key = @(S) sprintf('%d,',sort(S));
have = containers.Map('KeyType','char','ValueType','logical');
for h = find(isBinary).', have(key(find(C(h,:)))) = true; end
keep = false(1,0);
for h = 1:numel(P)
    kr = true(1,numel(P(h).jumps));
    nz = find(C(h,:) ~= 0);
    if ~isSingle(h) && ~isempty(kr) && all(ismember(nz,singleSensors))
        coef = P(h).jumps(:).' .* (muPilot(:).'*P(h).J);   % k * lambda_k at the pilot
        [cmax,kStar] = max(coef);
        if cmax > 0, kr(kStar) = false; end
    end
    km = true(1,r); if r >= 1, km(1) = false; end
    if isBinary(h)
        for q = 2:r
            if numel(nz) > q && allSubsetsPresent(nz,q,have,key), km(q) = false; end
        end
    end
    keep = [keep kr km]; %#ok<AGROW>
end
end

function ok = allSubsetsPresent(S,q,have,key)
% True if every nonempty subset of S with at most q elements is in the design.
ok = true;
for t = 1:q
    T = nchoosek(S,t);
    for k = 1:size(T,1)
        if ~isKey(have,key(T(k,:))), ok = false; return; end
    end
end
end

function P = buildProjections(Y,C,B)
% Per-projection data: projected counts, their unique values and the jump classes.
P = struct('z',{},'u',{},'w',{},'idx',{},'jumps',{},'J',{},'zmax',{});
for h = 1:size(C,1)
    z = Y*C(h,:).';
    [u,~,idx] = unique(z);
    jumps = unique(B(h,B(h,:)>0));
    J = double(B(h,:).' == jumps);                     % n x K: class j has jump jumps(k)
    P(h) = struct('z',z,'u',u,'w',accumarray(idx,1),'idx',idx,'jumps',jumps(:).','J',J,'zmax',max(z));
end
end

function G = simulatedConditions(ctx,mu,nMom,r)
% Conditions at mu on ctx.Ms intervals simulated from the model at mu: their
% covariance estimates S(mu) with a rank that does not depend on the sample size.
st = rng; rng(ctx.seed,'twister');
Ysim = poissrnd(repmat(max(mu(:),0).',ctx.Ms,1)) * double(ctx.A~=0).';
rng(st);
Psim = buildProjections(Ysim,ctx.C,ctx.B);
G = momentsPerObservation(Psim,mu,ctx.Ms,nMom,r);
end

function args = restrictArg(args,name,f)
% Apply f to the value of a name/value argument, if present.
for k = 1:2:numel(args)-1
    if (ischar(args{k}) || isstring(args{k})) && strcmpi(args{k},name) && ~isempty(args{k+1})
        args{k+1} = f(args{k+1});
    end
end
end

function C = uniqueNonzeroRows(C)
% Drop zero and duplicate rows of a projection matrix (after removing sensors).
C = C(any(C ~= 0,2),:);
[~,ia] = unique(C,'rows','stable'); C = C(ia,:);
end

function Q = cuCriterion(P,mu,N,nMom,n,r,o)
% Continuously updated GMM criterion: weights evaluated at mu itself.
if isfield(o,'CovSource') && strcmp(o.CovSource,'model')
    G1 = simulatedConditions(o.SimCtx,mu,nMom,r);
else
    G1 = momentsPerObservation(P,mu,N,nMom,r);
end
if isfield(o,'Keep'), G1 = G1(:,o.Keep); end
W = conditionWeights(G1,o);
[g,Jg] = meanMoments(P,mu,nMom,n,r); g = keepRows(o,g,Jg);
Q = g(:).'*W*g(:);
if ~isfinite(Q), Q = Inf; end
end
