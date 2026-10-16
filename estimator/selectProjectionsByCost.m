function [C,info] = selectProjectionsByCost(A,Y,muPilot,varargin)
%SELECTPROJECTIONSBYCOST Choose a cheap and informative set of projections.
%
%   [C,info] = selectProjectionsByCost(A,Y,muPilot)
%   [mu,~]   = estimateMuOptimizedProjections(A,Y,'ProjectionMatrix',C)
%
%   Candidates: binary sums over covered sensor subsets of size <= MaxBinary,
%   and integer projections on covered sensor groups S of size 2..MaxGroup with
%   weights L^0, L^1, ... on the sensors of S, for L = 2 (cheap, approximate)
%   and L = L_S, a few standard deviations above the group's counts (close to
%   the exact joint likelihood of Y_S, Section "A single integer projection").
%
%   Cost of a candidate c: (largest observed c'y + 1) x (number of jump classes),
%   the work of one evaluation of its recursion.
%
%   Two phases.
%   1. Identification: from the empty set, add the candidate with the largest
%      increase of log det(DEps*I + sum of Fisher informations) per
%      cost^CostExponent (the criterion of the binary selection, finite at every
%      step), until rank(R_C) = n. No Godambe criterion is used here.
%   2. Efficiency: from that identifying set, add the candidate with the largest
%      increase of log det G per cost^CostExponent, where G = H K^{-1} H is the
%      Godambe information of the equal-weight composite likelihood of the set:
%      H = sum of the candidates' Fisher informations, K = second moment of their
%      summed scores (redundant projections inflate K). Stops when the relative
%      gain per cost falls below TolGain, the total cost exceeds Budget, or after
%      MaxProjections.
%   H and K are estimated from scores at muPilot on either
%     'simulate' (default) Ms = max(SimulationSize, 20n) intervals drawn from the
%                model at muPilot: a model-based plug-in estimate of the design's
%                Godambe information at the pilot, of full rank whatever N;
%     'sample'   the observed intervals (CriterionSamples of them), with K shrunk
%                toward its diagonal (Ledoit-Wolf), since it has rank <= N-1.
%   A relative ridge (Ridge, 1e-12 of the mean eigenvalue) only guards against
%   rounding. G estimates a design criterion at the pilot; G <= I holds for the
%   population quantities at the pilot, not for the true rates.
%
%   Name/value (defaults): MaxBinary 3, MaxGroup 3, IncludeInteger true,
%   Budget Inf, TolGain 1e-3, MaxProjections 4*size(A,2), Ridge 1e-12,
%   Criterion 'simulate', SimulationSize 2000, Seed 1,
%   CostExponent 1 (gain per cost^CostExponent: 1 cost-aware, 0 cost ignored),
%   CriterionSamples 500 (random intervals used for H and K; costs use all data).
%   info: cost (per selected projection), totalCost, logdetPath, candidates
%   (number), isInteger (per selected projection), L (per selected projection).
p = inputParser;
addParameter(p,'MaxBinary',3); addParameter(p,'MaxGroup',3); addParameter(p,'IncludeInteger',true);
addParameter(p,'Budget',Inf); addParameter(p,'TolGain',1e-3); addParameter(p,'MaxProjections',[]);
addParameter(p,'Ridge',1e-12);
addParameter(p,'Criterion','simulate',@(x)any(strcmpi(x,{'simulate','sample'})));
addParameter(p,'SimulationSize',2000); addParameter(p,'Seed',1);
addParameter(p,'DEps',1e-6);   % epsilon of the summed-information criterion in phase 1
addParameter(p,'CostExponent',1);
addParameter(p,'CriterionSamples',500);   % intervals used to estimate the criterion (all if fewer)
parse(p,varargin{:}); o = p.Results;
A = double(A~=0); Y = double(Y); [m,n] = size(A); N = size(Y,1);
mu = max(double(muPilot(:)),1e-3*max(mean(muPilot),1e-3));
if isempty(o.MaxProjections), o.MaxProjections = 4*n; end

% ---- candidate family ----------------------------------------------------
cand = zeros(0,m); isInt = false(0,1); Lval = zeros(0,1);
for s = 1:min(o.MaxBinary,m)
    S = nchoosek(1:m,s);
    for k = 1:size(S,1)
        if any(all(A(S(k,:),:),1)) || s == 1
            c = zeros(1,m); c(S(k,:)) = 1; cand(end+1,:) = c; isInt(end+1,1) = false; Lval(end+1,1) = 1; %#ok<AGROW>
        end
    end
end
if o.IncludeInteger
    for s = 2:min(o.MaxGroup,m)
        S = nchoosek(1:m,s);
        for k = 1:size(S,1)
            g = S(k,:);
            if ~any(any(A(g,:),1) & sum(A(g,:),1) >= 2), continue; end   % some path crosses >= 2 sensors of g
            ymax = max(max(Y(:,g))); LS = ymax + 1 + ceil(3*sqrt(max(mean(Y(:,g),1))));
            for L = unique([2 LS])
                c = zeros(1,m); c(g) = L.^(0:s-1); cand(end+1,:) = c; isInt(end+1,1) = true; Lval(end+1,1) = L; %#ok<AGROW>
            end
        end
    end
end
% remove candidates with identical jump vectors (identical information)
B = cand*A; [~,ia] = unique(B,'rows','stable'); cand = cand(ia,:); isInt = isInt(ia); Lval = Lval(ia);
nc = size(cand,1);

% ---- per-candidate scores, information and cost ----------------------------
% Costs use the largest observed counts; the criterion uses scores on simulated
% intervals at the pilot ('simulate') or on observed ones ('sample').
if strcmpi(o.Criterion,'simulate')
    Ms = max(o.SimulationSize, 20*n);
    st = rng; rng(o.Seed,'twister');
    Ycrit = poissrnd(repmat(mu(:).',Ms,1)) * A.';
    rng(st);
else
    if N > o.CriterionSamples
        Ycrit = Y(unique(round(linspace(1,N,o.CriterionSamples))),:);   % evenly spaced (iid intervals)
    else
        Ycrit = Y;
    end
end
Nc = size(Ycrit,1);
U = cell(nc,1); H = cell(nc,1); cost = zeros(nc,1); JR = cell(nc,1); trH = zeros(nc,1);
for k = 1:nc
    c = cand(k,:); b = round(c*A);
    [U{k}, nJ] = candidateScores(b, Ycrit*c.', mu);
    H{k} = (U{k}.'*U{k})/Nc;                 % Fisher information (outer product of scores)
    trH(k) = trace(H{k}); JR{k} = jumpRows(c,A);
    cost(k) = (max(Y*c.')+1)*max(nJ,1);
end
cost = cost/min(cost);                       % relative to the cheapest candidate
w = cost.^o.CostExponent;                    % 1: gain per cost; 0: cost ignored

% ---- phase 1: identification by summed information per cost ---------------------
% log det(eps*I + sum of Fisher informations), as in the binary selection
% (Section 3.4), is finite at every step; gains per cost^CostExponent; continues
% until the jump-class matrix has rank n.
sel = []; avail = true(nc,1); total = 0; Rsel = zeros(0,n); rk = 0;
Hsum = zeros(n); ldD = sum(log(o.DEps*ones(n,1)));
while rk < n
    best = -Inf; kb = 0;
    for k = find(avail).'
        if rank([Rsel; JR{k}]) <= rk && rk > 0 && isempty(find(avail,1)), continue; end
        [R1,p1] = chol((Hsum+H{k}+(Hsum+H{k}).')/2 + o.DEps*eye(n));
        if p1, continue; end
        r = (2*sum(log(diag(R1))) - ldD)/w(k);
        if r > best, best = r; kb = k; end
    end
    if kb == 0, break; end
    sel(end+1) = kb; avail(kb) = false; total = total + cost(kb); %#ok<AGROW>
    Hsum = Hsum + H{kb}; [R1,~] = chol((Hsum+Hsum.')/2 + o.DEps*eye(n)); ldD = 2*sum(log(diag(R1)));
    Rsel = [Rsel; JR{kb}]; rk = rank(Rsel); %#ok<AGROW>
    if numel(sel) >= o.MaxProjections, break; end
end
nIdent = numel(sel); identified = rk == n;

% ---- phase 2: Godambe information gain per cost ------------------------------
Hs = zeros(n); Vs = zeros(Nc,n);
for k = sel, Hs = Hs + H{k}; Vs = Vs + U{k}; end
ld = criterion(Hs,Vs,o); path = ld;
while identified && numel(sel) < o.MaxProjections && any(avail)
    best = -Inf; kb = 0;
    for k = find(avail).'
        if total + cost(k) > o.Budget, continue; end
        r = (criterion(Hs + H{k}, Vs + U{k}, o) - ld)/w(k);
        if r > best, best = r; kb = k; end
    end
    if kb == 0, break; end
    gainKb = best*w(kb);
    if gainKb < o.TolGain*max(abs(ld),1)*w(kb), break; end
    sel(end+1) = kb; avail(kb) = false; %#ok<AGROW>
    Hs = Hs + H{kb}; Vs = Vs + U{kb}; ld = ld + gainKb; total = total + cost(kb); path(end+1) = ld; %#ok<AGROW>
end
C = cand(sel,:);
info = struct('cost',cost(sel),'totalCost',total,'logdetPath',path,'candidates',nc, ...
    'isInteger',isInt(sel),'L',Lval(sel),'rank',rk,'identified',identified, ...
    'identificationProjections',nIdent,'criterionIntervals',Nc,'criterion',lower(o.Criterion));
end

function ld = criterion(Hs,Vs,o)
% log det of the Godambe information H K^{-1} H, with K the second moment of the
% summed scores (their mean is zero under the model at the pilot). For observed
% intervals ('sample'), K is shrunk toward its diagonal (Ledoit-Wolf), since it
% has rank <= N-1. The relative ridge only guards against rounding.
[Nc,n] = size(Vs); K = (Vs.'*Vs)/Nc; K = (K+K.')/2;
if strcmpi(o.Criterion,'sample')
    Q2 = Vs.^2; varK = max((Q2.'*Q2)/Nc - K.^2, 0)/Nc;
    off = ~eye(n); den = sum(K(off).^2);
    if den > 0, delta = min(1,max(0,sum(varK(off))/den)); else, delta = 1; end
    K = (1-delta)*K + delta*diag(diag(K));
end
sH = max(trace(Hs)/n,realmin); sK = max(trace(K)/n,realmin);
[RH,pH] = chol((Hs+Hs.')/2 + o.Ridge*sH*eye(n)); [RK,pK] = chol(K + o.Ridge*sK*eye(n));
if pH || pK, ld = -Inf; return; end
ld = 4*sum(log(diag(RH))) - 2*sum(log(diag(RK)));
end

function [U,nJ] = candidateScores(b,z,mu)
% Per-interval score of the projection's log-likelihood at mu (N x n).
n = numel(mu); N = numel(z); K = max(b); lam = accumarray(b(b>0).',mu(b>0),[K 1]).';
jumps = find(lam > 0); nJ = numel(jumps);
zmax = max(z); lp = -inf(zmax+1,1); lp(1) = -sum(lam);
coef = (1:K).*lam;
for zz = 1:zmax
    k = jumps(jumps <= zz); if isempty(k), continue; end
    t = log(coef(k)) + lp(zz-k+1).'; mx = max(t); if ~isfinite(mx), continue; end
    lp(zz+1) = mx + log(sum(exp(t-mx))) - log(zz);
end
U = zeros(N,n);
for j = find(b > 0)
    ok = z - b(j) >= 0; r = zeros(N,1);
    r(ok) = exp(min(lp(z(ok)-b(j)+1) - lp(z(ok)+1),700));
    U(:,j) = r - 1;
end
end

function R = jumpRows(C,A)
% Jump-class matrix R_C of a projection family.
B = round(C*A); R = zeros(0,size(A,2));
for h = 1:size(B,1)
    for k = unique(B(h,B(h,:)>0)), R(end+1,:) = double(B(h,:)==k); end %#ok<AGROW>
end
end
