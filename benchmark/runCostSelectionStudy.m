function R = runCostSelectionStudy(varargin)
%RUNCOSTSELECTIONSTUDY Cost-aware projection selection against the default design and the GMM.
%
%   R = runCostSelectionStudy('Design','corridor','Replicates',20)
%   R = runCostSelectionStudy('Design','synthetic','SyntheticResults',file,'Replicates',20)
%
%   Designs:
%     'corridor'  m = 8 sensors along a corridor, n = 16 distinct path classes,
%                 each path crossing 1-4 consecutive sensors; all rates positive,
%                 drawn from Uniform(2,6) per replicate (instance fixed by Seed).
%     'synthetic' the class matrix and the true class rates of the synthetic study
%                 (results.A and results.trueClassRates(:,r) of synthetic_results.mat;
%                 most rates are zero); fresh Poisson data per replicate.
%   Methods (same data, same pilot = cumulant NNLS floored at 1e-3):
%     'default'   estimateMuOptimizedProjections with its own binary design
%     'alpha=a'   selectProjectionsByCost(...,'CostExponent',a,'MaxGroup',MaxGroup),
%                 then the composite likelihood of the selected projections,
%                 for every a in CostExponents (1: gain per cost, 0: cost ignored)
%     'GMM'       estimateMuEfficientGMM (default binary design)
%   Recorded per fit: relative error, seconds (selection included), number of
%   projections, number of integer projections, total relative cost.
%
%   Name/value (defaults): Design 'corridor', Replicates 20, N 1000,
%   CostExponents [1 0.5 0], MaxGroup 2, CriterionSamples 500, Seed 1,
%   SyntheticResults '' (path to synthetic_results.mat), UseParallel true,
%   OutputDir 'cost_selection_results'.
p = inputParser;
addParameter(p,'Design','corridor',@(x)any(strcmpi(x,{'corridor','synthetic'})));
addParameter(p,'Replicates',20); addParameter(p,'N',1000);
addParameter(p,'CostExponents',[1 0.5 0]); addParameter(p,'MaxGroup',2);
addParameter(p,'CriterionSamples',500); addParameter(p,'Seed',1);
addParameter(p,'SyntheticResults',''); addParameter(p,'UseParallel',true);
addParameter(p,'OutputDir','cost_selection_results');
parse(p,varargin{:}); o = p.Results;
o = guardOutputFields(o,{'OutputDir'});
outDir = char(o.OutputDir); if ~isfolder(outDir), mkdir(outDir); end

% ---- design ----------------------------------------------------------------
switch lower(o.Design)
    case 'corridor'
        rs = RandStream('mt19937ar','Seed',o.Seed); m = 8; n = 16; A = zeros(m,n); j = 0;
        while j < n
            Lp = randi(rs,4); s = randi(rs,m-Lp+1); col = zeros(m,1); col(s:s+Lp-1) = 1;
            if ~any(all(A(:,1:j) == col,1)), j = j+1; A(:,j) = col; end
        end
        rateFun = @(r) 2 + 4*rand(RandStream('mt19937ar','Seed',o.Seed+10*r),n,1);
    case 'synthetic'
        assert(~isempty(o.SyntheticResults) && isfile(o.SyntheticResults), ...
            'runCostSelectionStudy:SyntheticResults','Give the path of synthetic_results.mat.');
        d = load(o.SyntheticResults,'results'); A = double(d.results.A~=0);
        TR = d.results.trueClassRates; n = size(A,2);
        rateFun = @(r) TR(:,1+mod(r-1,size(TR,2)));
end
methods = [{'default'}, arrayfun(@(a)sprintf('alpha=%g',a),o.CostExponents,'UniformOutput',false), {'GMM'}];
nM = numel(methods); R0 = o.Replicates;
err = nan(R0,nM); secs = nan(R0,nM); nProj = nan(R0,nM); nInt = nan(R0,nM); cost = nan(R0,nM);
fprintf('Cost-selection study, %s design: %d sensors, %d classes, %d replicates, N = %d\n', ...
    o.Design,size(A,1),n,R0,o.N);
% True rates of every replicate, computed before the loop: a function handle
% cannot be called as rateFun(r) inside parfor (it would be read as slicing).
RATES = zeros(n,R0);
for r = 1:R0, RATES(:,r) = rateFun(r); end
useParfor = o.UseParallel && ~isempty(ver('parallel'));
if useParfor
    parfor r = 1:R0
        [err(r,:),secs(r,:),nProj(r,:),nInt(r,:),cost(r,:)] = oneReplicate(A,RATES(:,r),o,r); %#ok<PFBNS>
    end
else
    for r = 1:R0
        [err(r,:),secs(r,:),nProj(r,:),nInt(r,:),cost(r,:)] = oneReplicate(A,RATES(:,r),o,r);
        fprintf('  replicate %d: errors %s\n',r,mat2str(round(err(r,:)*1000)/1000));
    end
end
med = @(x) median(x,1,'omitnan');
summary = table(methods(:),med(err).',med(secs).',med(nProj).',med(nInt).',med(cost).', ...
    'VariableNames',{'Method','MedianError','MedianSeconds','MedianProjections','MedianInteger','MedianRelativeCost'});
R = struct('design',o.Design,'A',A,'methods',{methods},'err',err,'secs',secs,'nProj',nProj, ...
    'nInt',nInt,'cost',cost,'summary',summary,'settings',o);
save(fullfile(outDir,sprintf('cost_selection_%s.mat',lower(o.Design))),'R');
writetable(summary,fullfile(outDir,sprintf('cost_selection_%s_summary.csv',lower(o.Design))));
disp(summary);
end

function [e,s,np,ni,cst] = oneReplicate(A,mu0,o,r)
nM = 2 + numel(o.CostExponents);
e = nan(1,nM); s = nan(1,nM); np = nan(1,nM); ni = nan(1,nM); cst = nan(1,nM);
rng(o.Seed + 1000*r,'twister');
Y = poissrnd(repmat(mu0(:).',o.N,1)) * A.';
err = @(mu) norm(double(mu(:))-mu0(:))/max(norm(mu0),eps);
pilot = max(fitCumulantNNLS(A,Y,min(size(A,1),4)),1e-3);
% default projected CL
try
    t = tic; [mu,info] = estimateMuOptimizedProjections(A,Y,'Display','off');
    s(1) = toc(t); e(1) = err(mu); np(1) = size(info.C,1); ni(1) = 0;
catch ME, warning('CostSelection:Default','%s',ME.message); end
% cost-aware selections
for q = 1:numel(o.CostExponents)
    try
        t = tic;
        [C,si] = selectProjectionsByCost(A,Y,pilot,'CostExponent',o.CostExponents(q), ...
            'MaxGroup',o.MaxGroup,'CriterionSamples',o.CriterionSamples);
        mu = estimateMuOptimizedProjections(A,Y,'Display','off','ProjectionMatrix',C,'Mu0',pilot);
        s(1+q) = toc(t); e(1+q) = err(mu); np(1+q) = size(C,1); ni(1+q) = sum(si.isInteger); cst(1+q) = si.totalCost;
    catch ME, warning('CostSelection:Select','alpha=%g: %s',o.CostExponents(q),ME.message); end
end
% GMM
try
    t = tic; mu = estimateMuEfficientGMM(A,Y,'Mu0',pilot);
    s(nM) = toc(t); e(nM) = err(mu);
catch ME, warning('CostSelection:GMM','%s',ME.message); end
end
