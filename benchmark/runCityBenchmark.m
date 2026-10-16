function R = runCityBenchmark(varargin)
%RUNCITYBENCHMARK Methods compared on city-like networks, with a robustness sweep.
%
%   R = runCityBenchmark()
%   R = runCityBenchmark('Replicates',20,'UseParallel',true)
%
%   Instances (cityRoutingInstance): a GridSize x GridSize two-way street grid,
%   random origin-destination routes with stochastic route choice, sensors on
%   busy links; paths with equal sensor signatures are merged into classes and
%   unobservable paths dropped. The error ||mu_hat - mu||/||mu|| is on the class
%   rates. Path rates are RateScale times an irregular sequence in [0.5,1.5].
%
%   One factor varies at a time around the baseline (Sensors 8, Paths 40,
%   N 2000, no contamination):
%     sensors        SensorLevels        [2 4 8 16]
%     paths          PathLevels          [10 20 40 80]
%     observations   NLevels             [100 500 2000 10000]  (nested samples)
%     contamination  ContaminationLevels [0 0.005 0.02 0.05]   share of intervals with a
%                    spike of SpikeFactor (10) times the sensor mean at one sensor
%   Each replicate draws a new instance; all methods see the same data.
%
%   Methods: NNLS on the means; cumulant matching and bootstrap GLS through the
%   identifying order (at most 5; if more is needed they are rank deficient);
%   projected CL; projected GMM; pairwise CL (skipped above PairwiseMaxUnique
%   and PairwiseMaxSensors). IncludeRobust (false) adds the projected CL with the
%   contaminated likelihood (RobustEpsilon) and the robust projected GMM.
%   Recorded per fit: error, seconds, failure, rank deficiency, number of classes.
%
%   Output R.raw (table), R.summary (medians per sweep, level and method), saved
%   to OutputDir (MAT, CSV). Name/value: Replicates 100, RateScale 2, GridSize 10,
%   SpikeFactor 10, RobustEpsilon 0.01, BootstrapResamples 100,
%   PairwiseMaxUnique 1500, PairwiseMaxSensors 10 (the pairwise CL is skipped above
%   these), PairwiseMaxSeconds 30 (time cap per pairwise fit; a capped fit counts as
%   failed), Seed 1, UseParallel false, OutputDir 'city_benchmark_results'.
p = inputParser;
addParameter(p,'Replicates',100); addParameter(p,'RateScale',2); addParameter(p,'GridSize',10);
addParameter(p,'Baseline',struct('Sensors',8,'Paths',40,'N',2000,'Contamination',0));
addParameter(p,'SensorLevels',[2 4 8 16]); addParameter(p,'PathLevels',[10 20 40 80]);
addParameter(p,'NLevels',[100 500 2000 10000]); addParameter(p,'ContaminationLevels',[0 0.005 0.02 0.05]);
addParameter(p,'SpikeFactor',10); addParameter(p,'RobustEpsilon',0.01);
addParameter(p,'IncludeRobust',false);
addParameter(p,'BootstrapResamples',100); addParameter(p,'PairwiseMaxUnique',1500);
addParameter(p,'PairwiseMaxSensors',10); addParameter(p,'PairwiseMaxSeconds',30);
addParameter(p,'Seed',1); addParameter(p,'UseParallel',false);
addParameter(p,'OutputDir','city_benchmark_results');
parse(p,varargin{:}); o = p.Results;
% All output paths lie outside the library (see safeOutputPath).
o = guardOutputFields(o,{'OutputDir'});
methods = {'NNLS (means)','Cumulant matching','Bootstrap GLS','Projected CL','Projected GMM','Pairwise CL'};
if o.IncludeRobust          % robust variants, not used in the paper
    methods = [methods(1:4),{'Projected CL, robust'},methods(5),{'Projected GMM, robust'},methods(6)];
end
o.MethodNames = methods;
sweeps = {'sensors',o.SensorLevels; 'paths',o.PathLevels; 'observations',o.NLevels; ...
          'contamination',o.ContaminationLevels};
workers = 0; if o.UseParallel, workers = Inf; end
nM = numel(methods); rows = {}; firstErr = containers.Map();
for s = 1:size(sweeps,1)
    levels = sweeps{s,2};
    for l = 1:numel(levels)
        fprintf('City benchmark: %s = %g\n',sweeps{s,1},levels(l));
        cfgL = o.Baseline;
        switch sweeps{s,1}
            case 'sensors', cfgL.Sensors = levels(l);
            case 'paths', cfgL.Paths = levels(l);
            case 'observations', cfgL.N = levels(l);
            case 'contamination', cfgL.Contamination = levels(l);
        end
        nR = o.Replicates;
        E = nan(nR,nM); T = E; F = false(nR,nM); D = false(nR,nM); NC = nan(nR,1); IO = nan(nR,1); MSG = cell(nR,nM);
        parfor (r = 1:nR, workers)
            [E(r,:),T(r,:),F(r,:),D(r,:),NC(r),IO(r),MSG(r,:)] = cityReplicate(cfgL,r,o,s,l); %#ok<PFBNS>
        end
        for k = 1:nM                                     % first error message of each method
            mk = MSG(:,k); mk = mk(~cellfun(@isempty,mk));
            if ~isempty(mk) && ~isKey(firstErr,methods{k}), firstErr(methods{k}) = mk{1}; end
        end
        for r = 1:nR
            for k = 1:nM
                rows(end+1,:) = {sweeps{s,1},levels(l),r,methods{k},E(r,k),T(r,k),F(r,k),D(r,k),NC(r),IO(r)}; %#ok<AGROW>
            end
        end
    end
end
raw = cell2table(rows,'VariableNames',{'Sweep','Level','Replicate','Method','RelError', ...
    'Seconds','Failed','RankDeficient','Classes','IdentifyingOrder'});
raw.Sweep = string(raw.Sweep); raw.Method = string(raw.Method);
G = findgroups(raw.Sweep,raw.Level,raw.Method); nG = max(G); S = cell(nG,10);
for g = 1:nG
    sel = G==g; ok = sel & ~raw.Failed & isfinite(raw.RelError); f = find(sel,1);
    q = quantileLocal(raw.RelError(ok),[0.25 0.5 0.75]);
    S(g,:) = {raw.Sweep(f),raw.Level(f),raw.Method(f),q(2),q(1),q(3),median(raw.Seconds(ok)), ...
        mean(raw.Failed(sel)),mean(raw.RankDeficient(sel)),median(raw.Classes(sel))};
end
summary = cell2table(S,'VariableNames',{'Sweep','Level','Method','MedianError','Q1Error','Q3Error', ...
    'MedianSeconds','FailureRate','RankDeficientShare','MedianClasses'});
R = struct('methods',{methods},'sweeps',{sweeps},'raw',raw,'summary',summary,'settings',o);
R.errorMessages = firstErr;
reportErrors(firstErr);
folder = char(o.OutputDir); if ~isfolder(folder), mkdir(folder); end
save(fullfile(folder,'city_benchmark_results.mat'),'R','-v7.3');
writetable(raw,fullfile(folder,'city_benchmark_raw.csv'));
writetable(summary,fullfile(folder,'city_benchmark_summary.csv'));
fprintf('City benchmark saved to %s\n',folder);
end

% -------------------------------------------------------------------------
function [e,t,f,d,nc,io,msg] = cityReplicate(c,r,o,s,l)
names = o.MethodNames; nM = numel(names);
e = nan(1,nM); t = nan(1,nM); f = false(1,nM); d = false(1,nM); io = NaN; msg = cell(1,nM);
% Data seed: in the observation sweep the samples are nested (one draw of the
% largest N), and in the contamination sweep the same clean data are
% contaminated at each level, so these two sweeps do not depend on the level.
if s >= 3, seed = o.Seed*100003 + 1009*s + r; else, seed = o.Seed*100003 + 1009*s + 101*l + r; end
inst = cityRoutingInstance(c.Paths,c.Sensors,seedForInstance(o,c,r),'GridSize',o.GridSize);
A = inst.A; n = size(A,2); nc = n;
if n == 0, f(:) = true; msg(:) = {'no observable class'}; return; end
nP = inst.nPaths; relPath = 0.5 + mod((1:nP)'*0.6180339887 + 0.37*r,1);
muPath = o.RateScale*relPath/mean(relPath);
obs = inst.pathToClass > 0;
mu = accumarray(inst.pathToClass(obs),muPath(obs),[n 1]);
rs = RandStream('mt19937ar','Seed',seed);
Nmax = c.N; if isfield(o,'NLevels') && s == 3, Nmax = max(o.NLevels); end
X = zeros(Nmax,n); for j = 1:n, X(:,j) = poissinv(rand(rs,Nmax,1),mu(j)); end
Y = X*A.'; Y = Y(1:c.N,:);
if c.Contamination > 0
    Y = contaminateCounts(Y,c.Contamination,A*mu,o.SpikeFactor,rs);
end
ord = []; for k = 1:min(5,size(A,1)), if rank(cumulantSystem(A,k)) == n, ord = k; break; end, end
if isempty(ord), ordUse = min(5,size(A,1)); io = Inf; else, ordUse = ord; io = ord; end
rsB = RandStream('mt19937ar','Seed',seed+7);
try, fitCumulantNNLS(A,Y,1); catch, end      % untimed warm-up: first-call costs are not charged to method 1
rankDef = rank(A) < n;                          % rank check outside the timed block
for k = 1:nM
    t0 = tic;
    try
        switch names{k}
            case 'NNLS (means)', est = fitCumulantNNLS(A,Y,1); d(k) = rankDef;
            case 'Cumulant matching', est = fitCumulantNNLS(A,Y,ordUse); d(k) = isempty(ord);
            case 'Bootstrap GLS', est = fitBootstrapGLS(A,Y,ordUse,o.BootstrapResamples,rsB); d(k) = isempty(ord);
            case 'Projected CL'
                [est,info] = estimateMuOptimizedProjections(A,Y,'Display','off');
                d(k) = info.projectionDesignRank < n;
            case 'Projected CL, robust'
                [est,info] = estimateMuOptimizedProjections(A,Y,'Display','off','Contamination',o.RobustEpsilon);
                d(k) = info.projectionDesignRank < n;
            case 'Projected GMM', est = estimateMuEfficientGMM(A,Y);
            case 'Projected GMM, robust', est = estimateMuEfficientGMM(A,Y,'Robust',true,'Epsilon',o.RobustEpsilon);
            case 'Pairwise CL'
                if size(A,1) < 2, error('city:pairwise','needs two sensors'); end
                if size(A,1) > o.PairwiseMaxSensors, error('city:pairwise','skipped: more than PairwiseMaxSensors sensors'); end
                if maxUnique(Y) > o.PairwiseMaxUnique, error('city:pairwise','skipped: too many distinct count pairs'); end
                [est,ipw] = fitPairwiseCompositeLikelihood(A,Y,fitCumulantNNLS(A,Y,min(2,size(A,1))),'MaxSeconds',o.PairwiseMaxSeconds);
                if ipw.timeLimited, error('city:pairwise','time limit (PairwiseMaxSeconds)'); end
                d(k) = rank(cumulantSystem(A,2)) < n;
        end
        if any(~isfinite(est)), f(k) = true; msg{k} = 'non-finite estimate'; else, e(k) = norm(est(:)-mu)/norm(mu); end
    catch ME
        f(k) = true; msg{k} = sprintf('%s (%s)',ME.message,ME.identifier);
    end
    t(k) = toc(t0);
end
end

function sd = seedForInstance(o,c,r)
% The instance depends on the replicate and on the network factors only, so the
% observation and contamination sweeps reuse the baseline instances.
sd = o.Seed*7919 + 97*c.Sensors + 13*c.Paths + r;
end

function u = maxUnique(Y)
u = 0; m = size(Y,2);
for i = 1:m-1, for k = i+1:m, u = max(u,size(unique(Y(:,[i k]),'rows'),1)); end, end
end

function q = quantileLocal(v,pr)
v = sort(v(isfinite(v))); n = numel(v); q = nan(size(pr));
if n == 0, return; end
for k = 1:numel(pr)
    h = 1+(n-1)*pr(k); lo = floor(h); hi = ceil(h); q(k) = v(lo)+(h-lo)*(v(hi)-v(lo));
end
end

function reportErrors(firstErr)
% Print the first error message of every method that failed at least once.
k = keys(firstErr);
if isempty(k), return; end
fprintf('First error of each method that failed:\n');
for i = 1:numel(k), fprintf('  %-24s %s\n',k{i},firstErr(k{i})); end
end
