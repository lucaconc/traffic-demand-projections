function [mu,info] = estimateMuWeightedProjections(A,Y,varargin)
%ESTIMATEMUWEIGHTEDPROJECTIONS Projected composite likelihood with weights chosen
% from the covariance of the projection scores.
%
%   [mu,info] = estimateMuWeightedProjections(A,Y)                       % 'score' (default)
%   [mu,info] = estimateMuWeightedProjections(A,Y,'WeightUpdates','iterated')
%   [mu,info] = estimateMuWeightedProjections(A,Y,'Weights','scalar')
%
%   Both start from the equal-weight estimate muC of estimateMuOptimizedProjections
%   (options in 'PCLOptions') and keep its projection design C.
%
%   'scalar': one nonnegative weight per projection. With U_h the score of
%   projection h at muC for each interval, K_hk = mean of U_h'*U_k, and the
%   weights minimize 0.5*w'*K*w - w'*diag(K) over w >= 0 (the composite score
%   closest to the full-likelihood score, Huang and Ferrari 2024, with the sign
%   constraint so that the criterion stays a weighted likelihood). The composite
%   likelihood is then refitted with these weights from muC.
%
%   'onestep': one optimally weighted Newton step on the likelihood-ratio
%   conditions rho_{h,k}-1 of the projections (the components of the projection
%   scores), from muC: mu = muC - (G'WG)\(G'W gbar), with rates at zero whose step
%   is negative held at zero. W is the inverse covariance of the conditions; if
%   Hansen's J = N*gbar'*W*gbar exceeds df+10*sqrt(2*df)+50, W is recomputed with
%   the covariance shrunk toward its diagonal (Ledoit-Wolf).
%
%   'score' (default; alias 'fullscore'): minimizes mean(U)'*pinv(S)*mean(U),
%   the stacked score vector weighted by its covariance. WeightUpdates:
%   'twostep' estimates S at muC and once more at the first solution;
%   'iterated' re-estimates S at each new solution until mu stops changing
%   (at most 20 times), so that S is the covariance at the estimate itself.
%   Details: the composite score weighted by the full covariance of the
%   stacked projection scores U = (U_1,...,U_H), U_h = grad log p_{c_h}(Z_h):
%   solves G'*pinv(S)*mean(U) = 0 (Gauss-Newton with Levenberg-Marquardt
%   damping and an active set at zero), S and G re-estimated once at the
%   first solution. S is singular (each U_h repeats its likelihood-ratio
%   conditions once per class); its pseudo-inverse handles this exactly.
%
%   'equal': the equal-weight composite likelihood (muC itself).
%
%   info: muC, C, weights (scalar), weighting and Jstat (onestep), seconds.
ip = inputParser;
addParameter(ip,'Weights','score',@(x)any(strcmpi(x,{'score','fullscore','scalar','onestep','equal'})));
addParameter(ip,'WeightUpdates','twostep',@(x)any(strcmpi(x,{'twostep','iterated'})));
addParameter(ip,'MaxIter',100,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
addParameter(ip,'PCLOptions',{},@iscell);
addParameter(ip,'Ridge',1e-10,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
parse(ip,varargin{:}); o = ip.Results;
t0 = tic;
A = double(A); Y = double(Y); [N,~] = size(Y); n = size(A,2);
pclOpts = [{'Display','off'}, o.PCLOptions];
[muC,infoC] = estimateMuOptimizedProjections(A,Y,pclOpts{:});
muC = max(double(muC(:)),0); C = double(infoC.C); H = size(C,1);
info = struct('weightUpdates',0,'muC',muC,'C',C,'weights',ones(H,1)/H,'weighting','','Jstat',NaN, ...
    'pclInfo',infoC,'seconds',NaN,'method',lower(o.Weights));
if strcmpi(o.Weights,'equal'), mu = muC; info.seconds = toc(t0); return; end

% Per-projection probability tables and ratios at muC
T = projectionTables(A,Y,C,muC);

switch lower(o.Weights)
    case 'scalar'
        % K(h,k) = mean_r U_h(r,:)*U_k(r,:)', with U_h(r,j) = rho_{h,b_hj}(z_hr) - 1
        K = zeros(H);
        Uall = zeros(N,H,n);
        for h = 1:H, Uall(:,h,:) = reshape(T(h).score,N,1,n); end
        for j = 1:n
            Uj = Uall(:,:,j); K = K + Uj.'*Uj;
        end
        K = K/N; K = (K+K.')/2;
        d = diag(K);
        Kr = K + 1e-10*max(trace(K)/H,realmin)*eye(H);
        R = chol(Kr);
        w = lsqnonneg(R, R.'\d, optimset('Display','off'));   % min 0.5 w'Kw - w'd, w >= 0
        if ~any(w > 0), w = ones(H,1); end
        w = w/sum(w);
        [mu,infoW] = estimateMuOptimizedProjections(A,Y,pclOpts{:},'ProjectionMatrix',C, ...
            'ProjectionWeights',w,'Mu0',muC);
        mu = max(double(mu(:)),0);
        info.weights = w; info.refitInfo = infoW;
    case 'onestep'
        % Conditions rho_{h,k}-1 for every projection and positive jump size
        G = zeros(0,n); g = zeros(N,0);
        for h = 1:H
            ks = T(h).jumps;
            for a = 1:numel(ks)
                k = ks(a);
                g(:,end+1) = T(h).rho(:,k) - 1; %#ok<AGROW>
                % d rho_k / d mu_j = rho_{k+b_j} - rho_k rho_{b_j}  (classes with b_j > 0)
                row = zeros(1,n);
                for j = find(T(h).b(:).' > 0)
                    l = T(h).b(j);
                    row(j) = mean(T(h).rho(:,k+l) - T(h).rho(:,k).*T(h).rho(:,l));
                end
                G(end+1,:) = row; %#ok<AGROW>
            end
        end
        gbar = mean(g,1).'; M = numel(gbar);
        df = max(M-n,1); Jmax = df + 10*sqrt(2*df) + 50;
        W = weightsFrom(g,'none',o.Ridge); info.weighting = 'efficient';
        J = N*(gbar.'*W*gbar);
        if ~(J <= Jmax)
            W = weightsFrom(g,'lw',o.Ridge); info.weighting = 'lw (fallback)';
            J = N*(gbar.'*W*gbar);
        end
        info.Jstat = J; info.Jdf = df;
        grad = G.'*W*gbar;                                 % half-gradient of gbar'W gbar
        Hgn = G.'*W*G;
        free = muC > 0 | grad < 0;                        % rates at zero move only inward
        step = zeros(n,1);
        if any(free)
            Hf = Hgn(free,free); Hf = (Hf+Hf.')/2;
            Hf = Hf + 1e-10*max(trace(Hf)/nnz(free),realmin)*eye(nnz(free));
            step(free) = -(Hf\grad(free));
        end
        mu = max(muC + step,0);
    case {'score','fullscore'}
        mu = muC; info.weighting = 'pinv'; itTot = 0;
        if strcmpi(o.WeightUpdates,'iterated'), nW = 20; else, nW = 2; end
        for wit = 1:nW                                  % weights at the current solution, then refit
            muPrev = mu;
            T = projectionTables(A,Y,C,mu);
            [U,~] = stackedScores(T,n,N);
            W = pinvWeights(U);
            [mu,it,Q] = weightedGN(A,Y,C,mu,W,n,N,o.MaxIter);
            itTot = itTot + it;
            if wit >= 2 && norm(mu-muPrev) <= 1e-6*(1+norm(mu)), break; end
        end
        info.iterations = itTot; info.objective = Q; info.weightUpdates = wit;
end
info.seconds = toc(t0);
end

function T = projectionTables(A,Y,C,mu)
% For every projection: jump sizes b, observed counts z, and the ratios
% rho(:,k) = p(z-k)/p(z) for k = 1..2*max(b) (needed up to k+l), and the
% per-interval score U(:,j) = rho(:,b_j) - 1 (0 for classes with b_j = 0).
H = size(C,1); n = size(A,2);
T = struct('b',{},'jumps',{},'rho',{},'score',{});
for h = 1:H
    b = round(C(h,:)*A); z = Y*C(h,:).';
    Kmax = max(b); zmax = max(z);
    lp = logPmf(b,mu,zmax);
    rho = zeros(numel(z),2*Kmax);
    for k = 1:2*Kmax
        ok = z-k >= 0;
        lr = -inf(numel(z),1); lr(ok) = lp(z(ok)-k+1) - lp(z(ok)+1);
        rho(:,k) = exp(min(lr,700));
    end
    U = zeros(numel(z),n);
    for j = find(b > 0), U(:,j) = rho(:,b(j)) - 1; end
    T(h) = struct('b',b,'jumps',unique(b(b>0)),'rho',rho,'score',U);
end
end

function lp = logPmf(b,mu,zmax)
% log P(Z = z), z = 0..zmax, for Z = sum_j b_j X_j, X_j ~ Poisson(mu_j)
% (Panjer recursion in the log domain; compoundLogPmfMex when available).
K = max(b); lam = zeros(1,K);
for k = 1:K, lam(k) = sum(mu(b == k)); end
jumps = find(lam > 0); lamJ = lam(jumps);
if exist('compoundLogPmfMex','file') == 3
    lp = compoundLogPmfMex(lamJ(:),jumps(:),double(zmax)); return
end
lp = -inf(zmax+1,1); lp(1) = -sum(lamJ);
coef = zeros(1,K); coef(jumps) = jumps.*lamJ;
for z = 1:zmax
    k = 1:min(z,K); k = k(coef(k) > 0);
    if isempty(k), continue; end
    t = log(coef(k)) + lp(z-k+1).';
    m = max(t); if ~isfinite(m), continue; end
    lp(z+1) = m + log(sum(exp(t-m))) - log(z);
end
end

function W = weightsFrom(G,mode,ridge)
% Inverse covariance of the conditions G (N x M): 'none' efficient (eigenvalue
% floor ridge), 'lw' shrunk toward the diagonal (Ledoit-Wolf intensity);
% conditions without variation get no weight.
N = size(G,1); M = size(G,2);
Gc = G - mean(G,1); S = (Gc.'*Gc)/max(N-1,1); S = (S+S.')/2;
v = diag(S); live = isfinite(v) & v > 1e-12*max([v(isfinite(v)); realmin]);
W = zeros(M); if ~any(live), return; end
Sl = S(live,live); Gl = Gc(:,live);
if strcmpi(mode,'lw')
    Q2 = Gl.^2; varS = max((Q2.'*Q2)/N - ((Gl.'*Gl)/N).^2,0)/N;
    off = ~eye(size(Sl)); den = sum(Sl(off).^2);
    if den > 0, delta = min(1,max(0,sum(varS(off))/den)); else, delta = 1; end
    Sl = (1-delta)*Sl + delta*diag(diag(Sl));
end
sd = sqrt(diag(Sl)); Rm = Sl./(sd*sd.');
[V,D] = eig((Rm+Rm.')/2); d = diag(D); d = max(d,ridge*max(d));
W(live,live) = (V*diag(1./d)*V.')./(sd*sd.');
W = (W+W.')/2;
end

function [U,G] = stackedScores(T,n,N)
% Stacked projection scores U (N x sum_h n_h), U_hj = rho_{h,b_hj} - 1 for the
% classes with b_hj > 0, and the mean derivative G (columns: mu).
U = zeros(N,0); G = zeros(0,n);
for h = 1:numel(T)
    js = find(T(h).b(:).' > 0);
    for j = js
        k = T(h).b(j);
        U(:,end+1) = T(h).rho(:,k) - 1; %#ok<AGROW>
        if nargout > 1
            row = zeros(1,n);
            for i = js
                l = T(h).b(i);
                row(i) = mean(T(h).rho(:,k+l) - T(h).rho(:,k).*T(h).rho(:,l));
            end
            G(end+1,:) = row; %#ok<AGROW>
        end
    end
end
end

function W = pinvWeights(U)
% Pseudo-inverse of the covariance of the conditions, on the correlation scale.
N = size(U,1); M = size(U,2);
Uc = U - mean(U,1); S = (Uc.'*Uc)/max(N-1,1); S = (S+S.')/2;
v = diag(S); live = isfinite(v) & v > 1e-12*max([v(isfinite(v)); realmin]);
W = zeros(M); if ~any(live), return; end
sd = sqrt(v(live)); Rm = S(live,live)./(sd*sd.');
[V,D] = eig((Rm+Rm.')/2); d = diag(D); keep = d > 1e-10*max(d);
W(live,live) = (V(:,keep)*diag(1./d(keep))*V(:,keep).')./(sd*sd.');
W = (W+W.')/2;
end

function [mu,it,Q] = weightedGN(A,Y,C,mu,W,n,N,maxIter)
% Minimize mean(U)'*W*mean(U) over mu >= 0 by projected Levenberg-Marquardt.
T = projectionTables(A,Y,C,mu); [U,G] = stackedScores(T,n,N);
g = mean(U,1).'; Q = g.'*W*g; lam = 1e-3; it = 0;
for it = 1:maxIter
    grad = G.'*W*g; Hgn = G.'*W*G;
    act = (mu <= 0) & (grad > 0); free = ~act; improved = false;
    while lam < 1e12
        step = zeros(n,1);
        Hf = Hgn(free,free); Hf = (Hf+Hf.')/2; dd = max(diag(Hf),1e-8*max([diag(Hf);eps]));
        step(free) = -((Hf + lam*diag(dd))\grad(free));
        cand = max(mu+step,0);
        Tc = projectionTables(A,Y,C,cand); Uc = stackedScores(Tc,n,N); gc = mean(Uc,1).';
        Qc = gc.'*W*gc;
        if isfinite(Qc) && Qc < Q
            dmu = norm(cand-mu); mu = cand; Q = Qc; T = Tc;
            [U,G] = stackedScores(T,n,N); g = mean(U,1).';
            lam = max(lam/10,1e-12); improved = true; break
        end
        lam = lam*10;
    end
    if ~improved || dmu <= 1e-9*(1+norm(mu)), break; end
end
end
