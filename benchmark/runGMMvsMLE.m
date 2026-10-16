function R = runGMMvsMLE(varargin)
%RUNGMMVSMLE How closely the estimators approximate the maximum-likelihood estimator.
%
%   R = runGMMvsMLE()
%   R = runGMMvsMLE('Replicates',100,'UseParallel',true)
%
%   On designs whose exact likelihood is computable (fitFullLikelihood): the easy
%   (2 sensors, 3 paths) and intermediate (3 sensors, 7 paths) benchmark matrices,
%   at rate scales RateScales, for sample sizes SampleSizes (nested within each
%   replicate). All methods see the same data:
%     MLE            exact maximum likelihood (started at the GMM estimate)
%     Projected GMM  estimateMuEfficientGMM
%     Projected CL   estimateMuOptimizedProjections
%     Bootstrap GLS  cumulants through the identifying order
%     Cumulants      NNLS on the same cumulants
%   Reported per case, scale, N and method:
%     RMSEratio      RMSE of the method / RMSE of the MLE (1 = as efficient)
%     DistanceRatio  median over replicates of ||mu_method - mu_MLE|| / ||mu_MLE - mu||
%                    (how far the method is from the MLE, relative to the MLE's own
%                    error; it tends to 0 for an estimator that is asymptotically
%                    equivalent to the MLE)
%   Saved to OutputDir (MAT, CSV). Name/value: Replicates 100, SampleSizes
%   [250 1000 4000 16000], RateScales [1 3], BootstrapResamples 200, Seed 1,
%   UseParallel false, OutputDir 'mle_comparison_results'.
p = inputParser;
addParameter(p,'Replicates',100); addParameter(p,'SampleSizes',[250 1000 4000 16000]);
addParameter(p,'RateScales',[1 3]); addParameter(p,'BootstrapResamples',200);
addParameter(p,'Seed',1); addParameter(p,'UseParallel',false);
addParameter(p,'OutputDir','mle_comparison_results');
parse(p,varargin{:}); o = p.Results;
% All output paths lie outside the library (see safeOutputPath).
o = guardOutputFields(o,{'OutputDir'});
cases = benchmarkRoutingCases(); cases = cases(1:2);
methods = {'MLE','Projected GMM','Projected CL','Bootstrap GLS','Cumulants'};
Ns = sort(o.SampleSizes(:).'); nM = numel(methods); nR = o.Replicates;
workers = 0; if o.UseParallel, workers = Inf; end
rows = {};
for c = 1:numel(cases)
    for s = 1:numel(o.RateScales)
        fprintf('MLE comparison: %s, rate scale %g\n',cases(c).name,o.RateScales(s));
        cs = cases(c); sc = o.RateScales(s);
        EST = cell(nR,1); MU = sc*cs.relRates(:);
        parfor (r = 1:nR, workers)
            EST{r} = mleReplicate(cs,MU,Ns,o,c,s,r); %#ok<PFBNS>
        end
        for a = 1:numel(Ns)
            E = nan(nR,nM); D = nan(nR,nM);
            for r = 1:nR
                est = EST{r}{a};
                if isempty(est) || any(~isfinite(est(:,1))), continue; end
                for k = 1:nM
                    if all(isfinite(est(:,k)))
                        E(r,k) = norm(est(:,k)-MU)/norm(MU);
                        D(r,k) = norm(est(:,k)-est(:,1))/max(norm(est(:,1)-MU),eps);
                    end
                end
            end
            rm = sqrt(mean(E.^2,1,'omitnan'));
            for k = 1:nM
                rows(end+1,:) = {cs.name,sc,Ns(a),methods{k},rm(k),rm(k)/rm(1), ...
                    median(D(:,k),'omitnan'),mean(isfinite(E(:,k)))}; %#ok<AGROW>
            end
        end
    end
end
summary = cell2table(rows,'VariableNames',{'Case','RateScale','N','Method','RMSE', ...
    'RMSEratio','DistanceRatio','SuccessRate'});
summary.Case = string(summary.Case); summary.Method = string(summary.Method);
R = struct('summary',summary,'methods',{methods},'settings',o);
folder = char(o.OutputDir); if ~isfolder(folder), mkdir(folder); end
save(fullfile(folder,'mle_comparison_results.mat'),'R','-v7.3');
writetable(summary,fullfile(folder,'mle_comparison_summary.csv'));
disp(summary(summary.Method~="MLE",:));
end

function out = mleReplicate(cs,mu,Ns,o,c,s,r)
% Estimates of all methods (columns: MLE, GMM, PCL, GLS, cumulants) for each N.
A = cs.A; n = numel(mu); rs = RandStream('mt19937ar','Seed',o.Seed + 1e6*c + 1e4*s + r);
X = poissinv(rand(rs,max(Ns),n),repmat(mu.',max(Ns),1)); Yall = X*A.';
out = cell(1,numel(Ns));
for a = 1:numel(Ns)
    Y = Yall(1:Ns(a),:); est = nan(n,5);
    try, est(:,2) = estimateMuEfficientGMM(A,Y); catch, end
    start = est(:,2); if any(~isfinite(start)), start = fitCumulantNNLS(A,Y,cs.identifyingOrder); end
    try, est(:,1) = fitFullLikelihood(A,Y,start); catch, end
    try, est(:,3) = estimateMuOptimizedProjections(A,Y,'Display','off'); catch, end
    try, est(:,4) = fitBootstrapGLS(A,Y,cs.identifyingOrder,o.BootstrapResamples,rs); catch, end
    try, est(:,5) = fitCumulantNNLS(A,Y,cs.identifyingOrder); catch, end
    out{a} = est;
end
end
