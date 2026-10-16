function R = runIdentificationBenchmark(varargin)
%RUNIDENTIFICATIONBENCHMARK Methods compared under controlled identification difficulty.
%
%   R = runIdentificationBenchmark()
%   R = runIdentificationBenchmark('Replicates',50,'SampleSizes',[100 1000],'UseParallel',true)
%
%   Routing matrices (benchmarkRoutingCases): easy (running example), intermediate
%   (all 7 nonzero columns of {0,1}^3), difficult (all 15 of {0,1}^4), complex
%   (all 31 of {0,1}^5). Rates are
%   strictly positive, with fixed relative intensities times each RateScale.
%   Every replicate draws Nmax intervals once; each sample size uses the first N
%   (nested), and every method sees exactly the same data (paired comparison).
%
%   Methods (R.methods):
%     NNLS on the sensor means           (cumulant system through order 1)
%     cumulant matching, orders 2 to 5  (unweighted NNLS on sample joint cumulants;
%                                         only orders <= number of sensors)
%     bootstrap-weighted GLS             (cumulants through the identifying order)
%     projected composite likelihood     (estimateMuOptimizedProjections, defaults)
%     projected GMM                      (estimateMuEfficientGMM: efficiently weighted
%                                         projected scores, mean and variance conditions)
%     pairwise composite likelihood      (exact bivariate likelihoods of all sensor
%                                         pairs; skipped, as not feasible, when a pair
%                                         has more than PairwiseMaxUnique distinct values)
%   R.rank(c,k) is the rank of method k's identifying system for case c; a method
%   with rank below the number of paths cannot identify the rates (its error floor
%   reflects missing information, not poor optimization).
%
%   Output: relative error ||mu_hat-mu||/||mu||, seconds, failure (error or
%   non-finite estimate) and feasibility per case, scale, N, replicate and method;
%   R.raw (table), R.summary (medians), saved to OutputDir as MAT and CSV.
%
%   Rate profiles (RateProfiles): 'balanced' (relative rates between 0.5 and 1.5)
%   and 'heavyPairs' (rates of the paths crossing exactly two sensors five times
%   larger, renormalized to mean one). A rank-deficient method is flat along its
%   unidentified directions, where its estimate is set by the optimizer (an
%   interior-point method ends near the analytic centre of the flat set); with
%   balanced rates that centre can lie close to the truth, so the second profile
%   exposes the error floors of such methods.
%
%   Name/value (defaults): Cases {'easy','intermediate','difficult'} (add 'complex' for
%   the 5-sensor, 31-path matrix), Replicates 100, SampleSizes [100 500 2000 10000],
%   RateScales [2 10], BootstrapResamples 200, PairwiseMaxSensors 4 (the pairwise CL
%   is skipped on larger designs), PairwiseMaxUnique 1500 (distinct
%   count pairs per sensor pair above which the pairwise CL is skipped),
%   PairwiseMaxSeconds 30 (time cap per pairwise fit; a capped fit counts as failed),
%   Seed 1, UseParallel false, OutputDir 'benchmark_results'.

p = inputParser;
addParameter(p,'Replicates',100,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
addParameter(p,'SampleSizes',[100 500 2000 10000],@(x)isnumeric(x)&&all(x>=10));
addParameter(p,'RateScales',[2 10],@(x)isnumeric(x)&&all(x>0));
addParameter(p,'RateProfiles',{'balanced','heavyPairs'},@(x)iscell(x)||isstring(x));
addParameter(p,'BootstrapResamples',200,@(x)isnumeric(x)&&isscalar(x)&&x>=20);
addParameter(p,'PairwiseMaxUnique',1500,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
addParameter(p,'PairwiseMaxSeconds',30,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'PairwiseMaxSensors',4,@(x)isnumeric(x)&&isscalar(x)&&x>=2);
addParameter(p,'CostExponent',0.5,@(x)isnumeric(x)&&isscalar(x)&&x>=0);   % cost-aware selection (method 10)
addParameter(p,'CostMaxGroup',2,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
addParameter(p,'Cases',{'easy','intermediate','difficult'},@(x)iscell(x)||isstring(x));
addParameter(p,'Seed',1,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'UseParallel',false,@(x)islogical(x)||isnumeric(x));
addParameter(p,'OutputDir','benchmark_results',@(x)ischar(x)||isstring(x));
parse(p,varargin{:}); o = p.Results;
% All output paths lie outside the library (see safeOutputPath).
o = guardOutputFields(o,{'OutputDir'});
Ns = sort(o.SampleSizes(:).'); Nmax = max(Ns);
cases = benchmarkRoutingCases();
cases = cases(ismember(string({cases.name}),string(o.Cases)));     % default: easy, intermediate, difficult
assert(~isempty(cases),'runIdentificationBenchmark:Cases','No routing matrix selected (Cases).');
methods = {'NNLS (means)','Cumulants, order 2','Cumulants, order 3','Cumulants, order 4', ...
           'Cumulants, order 5','Bootstrap GLS','Projected CL','Projected GMM','Pairwise CL','Cost-aware CL'};
combos = struct('profile',{},'scale',{});
for pr = string(o.RateProfiles(:)).'
    for sc_ = o.RateScales(:).', combos(end+1) = struct('profile',char(pr),'scale',sc_); end %#ok<AGROW>
end
nM = numel(methods); nC = numel(cases); nS = numel(combos); nN = numel(Ns); nR = o.Replicates;
% Ranks of the identifying systems and applicability
rnk = nan(nC,nM); applicable = true(nC,nM);
for c = 1:nC
    A = cases(c).A; m = size(A,1);
    rnk(c,1) = rank(cumulantSystem(A,1));
    for ord = 2:5
        applicable(c,ord) = ord <= m;
        if applicable(c,ord), rnk(c,ord) = rank(cumulantSystem(A,ord)); end
    end
    rnk(c,6) = rank(cumulantSystem(A,cases(c).identifyingOrder));
    rnk(c,7) = NaN;                                  % filled from the fits below
    rnk(c,8) = rank(cumulantSystem(A,max(2,cases(c).identifyingOrder)));
    rnk(c,9) = rank(cumulantSystem(A,2));
    rnk(c,10) = size(A,2);                                    % the selection completes the rank
end
workers = 0; if o.UseParallel, workers = Inf; end
firstErr = containers.Map();
ERR = nan(nC,nS,nN,nR,nM); SEC = ERR; FAIL = false(size(ERR)); FEAS = true(size(ERR)); PCLRANK = nan(nC,nS,nN,nR);
for c = 1:nC
    for s = 1:nS
        fprintf('Benchmark: %s, %s rates, scale %g\n',cases(c).name,combos(s).profile,combos(s).scale);
        e = nan(nR,nN,nM); t = e; f = false(nR,nN,nM); fe = true(nR,nN,nM); pr = nan(nR,nN); ms = cell(nR,nM);
        cs = cases(c); sc = combos(s).scale; cs.relRates = profileRates(cases(c),combos(s).profile);
        parfor (r = 1:nR, workers)
            [e(r,:,:),t(r,:,:),f(r,:,:),fe(r,:,:),pr(r,:),ms(r,:)] = benchmarkReplicate(cs,sc,c,s,r,Ns,Nmax,o,applicable(c,:)); %#ok<PFBNS>
        end
        ERR(c,s,:,:,:) = permute(e,[2 1 3]); SEC(c,s,:,:,:) = permute(t,[2 1 3]);
        FAIL(c,s,:,:,:) = permute(f,[2 1 3]); FEAS(c,s,:,:,:) = permute(fe,[2 1 3]);
        PCLRANK(c,s,:,:) = permute(pr,[2 1]);
        for m = 1:nM
            mk = ms(:,m); mk = mk(~cellfun(@isempty,mk));
            if ~isempty(mk) && ~isKey(firstErr,methods{m}), firstErr(methods{m}) = mk{1}; end
        end
    end
    v = PCLRANK(c,:,:,:); rnk(c,7) = max(v(:));
end
% Long table
rows = cell(nC*nS*nN*nR*nM,9); k = 0;
for c = 1:nC, for s = 1:nS, for n = 1:nN, for r = 1:nR, for m = 1:nM %#ok<ALIGN>
    if ~applicable(c,m), continue; end
    k = k+1;
    rows(k,:) = {cases(c).name,combos(s).profile,combos(s).scale,Ns(n),r,methods{m},ERR(c,s,n,r,m),SEC(c,s,n,r,m), ...
        FAIL(c,s,n,r,m) || ~FEAS(c,s,n,r,m)};
end, end, end, end, end
rows = rows(1:k,:);
raw = cell2table(rows,'VariableNames',{'Case','Profile','RateScale','N','Replicate','Method','RelError','Seconds','Failed'});
raw.Case = string(raw.Case); raw.Method = string(raw.Method); raw.Profile = string(raw.Profile);
% Summary: medians and quartiles over replicates (successful fits only)
G = findgroups(raw.Case,raw.Profile,raw.RateScale,raw.N,raw.Method);
okFit = ~raw.Failed & isfinite(raw.RelError);
nG = max(G); S = cell(nG,10);
for g = 1:nG
    sel = G==g; ok = sel & okFit; first = find(sel,1);
    q = quantileLocal(raw.RelError(ok),[0.25 0.5 0.75]);
    S(g,:) = {raw.Case(first),raw.Profile(first),raw.RateScale(first),raw.N(first),raw.Method(first), ...
        q(2),q(1),q(3),median(raw.Seconds(ok)),mean(raw.Failed(sel))};
end
summary = cell2table(S,'VariableNames',{'Case','Profile','RateScale','N','Method','MedianError','Q1Error', ...
    'Q3Error','MedianSeconds','FailureRate'});
R = struct('cases',cases,'methods',{methods},'rank',rnk,'applicable',applicable, ...
    'sampleSizes',Ns,'rateScales',o.RateScales,'rateProfiles',{cellstr(string(o.RateProfiles))},'replicates',nR,'bootstrap',o.BootstrapResamples, ...
    'raw',raw,'summary',summary,'settings',o);
R.errorMessages = firstErr;
k_ = keys(firstErr);
if ~isempty(k_), fprintf('First error of each method that failed:\n'); end
for i_ = 1:numel(k_), fprintf('  %-24s %s\n',k_{i_},firstErr(k_{i_})); end
folder = char(o.OutputDir); if ~isfolder(folder), mkdir(folder); end
save(fullfile(folder,'benchmark_results.mat'),'R','-v7.3');
writetable(raw,fullfile(folder,'benchmark_raw.csv'));
writetable(summary,fullfile(folder,'benchmark_summary.csv'));
fprintf('Benchmark saved to %s\n',folder);
end

% -------------------------------------------------------------------------
function [e,t,f,fe,pr,msg] = benchmarkReplicate(cs,scale,c,s,r,Ns,Nmax,o,applicable)
% One replicate: nested samples, all methods on the same data.
nN = numel(Ns); nM = 10;                        % keep equal to numel(methods)
e = nan(nN,nM); t = nan(nN,nM); f = false(nN,nM); fe = true(nN,nM); pr = nan(1,nN); msg = cell(1,nM);
A = cs.A; mu = scale*cs.relRates(:); n = numel(mu);
rs = RandStream('mt19937ar','Seed',o.Seed + 1e6*c + 1e4*s + r);
X = zeros(Nmax,n);
for j = 1:n, X(:,j) = poissonDraw(rs,mu(j),Nmax); end
Yall = X*A.';
for a = 1:nN
    Y = Yall(1:Ns(a),:);
    rsB = RandStream('mt19937ar','Seed',o.Seed + 7e6 + 1e6*c + 1e4*s + 10*r + a);
    try, fitCumulantNNLS(A,Y,1); catch, end   % untimed warm-up: first-call costs are not charged to method 1
    for m = 1:nM
        if ~applicable(m), continue; end
        t0 = tic;
        try
            switch m
                case 1, est = fitCumulantNNLS(A,Y,1);
                case {2,3,4,5}, est = fitCumulantNNLS(A,Y,m);
                case 6, est = fitBootstrapGLS(A,Y,cs.identifyingOrder,o.BootstrapResamples,rsB);
                case 7
                    [est,info] = estimateMuOptimizedProjections(A,Y,'Display','off');
                    pr(a) = info.projectionDesignRank;
                case 8, est = estimateMuEfficientGMM(A,Y);
                case 9
                    if size(A,1) > o.PairwiseMaxSensors
                        fe(a,m) = false; if isempty(msg{m}), msg{m} = 'skipped: more than PairwiseMaxSensors sensors'; end; continue
                    end
                    if maxUniquePairs(Y) > o.PairwiseMaxUnique
                        fe(a,m) = false; if isempty(msg{m}), msg{m} = 'skipped: too many distinct count pairs (PairwiseMaxUnique)'; end; continue
                    end
                    [est,ipw] = fitPairwiseCompositeLikelihood(A,Y,fitCumulantNNLS(A,Y,2),'MaxSeconds',o.PairwiseMaxSeconds);
                    if ipw.timeLimited
                        f(a,m) = true; if isempty(msg{m}), msg{m} = 'time limit (PairwiseMaxSeconds)'; end
                        t(a,m) = toc(t0); continue
                    end
                case 10                                  % cost-aware selection with integer groups
                    pilot = max(fitCumulantNNLS(A,Y,min(size(A,1),4)),1e-3);
                    Csel = selectProjectionsByCost(A,Y,pilot,'CostExponent',o.CostExponent,'MaxGroup',o.CostMaxGroup);
                    est = estimateMuOptimizedProjections(A,Y,'Display','off','ProjectionMatrix',Csel,'Mu0',pilot);
            end
            if any(~isfinite(est)), f(a,m) = true; msg{m} = 'non-finite estimate'; else, e(a,m) = norm(est(:)-mu)/norm(mu); end
        catch ME
            f(a,m) = true; if isempty(msg{m}), msg{m} = sprintf('%s (%s)',ME.message,ME.identifier); end
        end
        t(a,m) = toc(t0);
    end
end
end

function u = maxUniquePairs(Y)
u = 0; m = size(Y,2);
for i = 1:m-1, for k = i+1:m, u = max(u,size(unique(Y(:,[i k]),'rows'),1)); end, end
end

function x = poissonDraw(rs,lam,n)
% Poisson draws by inversion from the replicate's own stream (vectorized and
% independent of the global generator, so replicates are reproducible in parfor).
x = poissinv(rand(rs,n,1),lam);
end

function q = quantileLocal(v,pr)
v = sort(v(isfinite(v))); n = numel(v); q = nan(size(pr));
if n == 0, return; end
for k = 1:numel(pr)
    h = 1+(n-1)*pr(k); lo = floor(h); hi = ceil(h); q(k) = v(lo)+(h-lo)*(v(hi)-v(lo));
end
end

function r = profileRates(cs,profile)
% Relative path rates (mean one) of a rate profile.
r = cs.relRates(:);
if strcmpi(profile,'heavyPairs')
    two = sum(cs.A,1).' == 2;            % paths crossing exactly two sensors
    r(two) = 5*r(two); r = r/mean(r);
end
end
