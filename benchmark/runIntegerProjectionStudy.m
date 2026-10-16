function R = runIntegerProjectionStudy(varargin)
%RUNINTEGERPROJECTIONSTUDY Integer projections between composite and full likelihood.
%
%   R = runIntegerProjectionStudy('Replicates',50)
%
%   For the easy and intermediate matrices of the identification benchmark, each
%   replicate draws N observations at rate scale RateScale and fits
%     - the exact maximum-likelihood estimate muML (fitFullLikelihood), and
%     - the composite likelihood of the single projection Z = c'Y with
%       c = (1, L, L^2, ..., L^(m-1)), for every L in Ls and for L* = max(Y)+1.
%       Every L >= 2 identifies mu (distinct base-L jump sizes). From L* on, no
%       two observed count vectors share a projected value, and the projected
%       probability of an observation exceeds the full one by at most
%       P(max_i Y_i >= L): the likelihood approaches the full likelihood.
%   Recorded per fit: relative error ||mu-mu0||/||mu0||, distance to muML in units
%   of the ML error ||mu-muML||/||muML-mu0||, loss in full log-likelihood
%   loglik(muML)-loglik(mu) (per observation), and seconds.
%
%   Name/value (defaults): Replicates 50, N 2000, RateScale 3, Ls [2 3 4 6 8 12 16],
%   Cases {'easy','intermediate'}, Seed 1, UseParallel true, OutputDir.
%   R.summary: one row per case and L (L = Inf denotes L*), median values.
p = inputParser;
addParameter(p,'Replicates',50); addParameter(p,'N',2000); addParameter(p,'RateScale',3);
addParameter(p,'Ls',[2 3 4 6 8 12 16]); addParameter(p,'Cases',{'easy','intermediate'});
addParameter(p,'Seed',1); addParameter(p,'UseParallel',true);
addParameter(p,'OutputDir','integer_projection_results');
parse(p,varargin{:}); o = p.Results;
o = guardOutputFields(o,{'OutputDir'});
outDir = char(o.OutputDir); if ~isfolder(outDir), mkdir(outDir); end
cases = benchmarkRoutingCases();
cases = cases(ismember(string({cases.name}),string(o.Cases)));
Ls = [sort(o.Ls(:)).', Inf];                 % Inf: the one-to-one L* of each replicate
nL = numel(Ls); R0 = o.Replicates; nC = numel(cases);
err = nan(nC,R0,nL); dist = nan(nC,R0,nL); loss = nan(nC,R0,nL); secs = nan(nC,R0,nL);
errML = nan(nC,R0); secML = nan(nC,R0); Lstar = nan(nC,R0);
for c = 1:nC
    A = double(cases(c).A); mu0 = o.RateScale*cases(c).relRates(:); m = size(A,1);
    e = nan(R0,nL); d = nan(R0,nL); l = nan(R0,nL); s = nan(R0,nL); eM = nan(R0,1); sM = nan(R0,1); Lst = nan(R0,1);
    useParfor = o.UseParallel && ~isempty(ver('parallel'));
    if useParfor
        parfor r = 1:R0
            [e(r,:),d(r,:),l(r,:),s(r,:),eM(r),sM(r),Lst(r)] = oneReplicate(A,mu0,m,Ls,o,c,r); %#ok<PFBNS>
        end
    else
        for r = 1:R0
            [e(r,:),d(r,:),l(r,:),s(r,:),eM(r),sM(r),Lst(r)] = oneReplicate(A,mu0,m,Ls,o,c,r);
        end
    end
    err(c,:,:) = e; dist(c,:,:) = d; loss(c,:,:) = l; secs(c,:,:) = s; errML(c,:) = eM; secML(c,:) = sM; Lstar(c,:) = Lst;
    fprintf('%s: ML error %.4f (%.2f s); one-to-one L* median %g\n',cases(c).name,median(eM,'omitnan'),median(sM,'omitnan'),median(Lst,'omitnan'));
end
% Summary: median over replicates
rows = {};
for c = 1:nC
    for k = 1:nL
        Lk = Ls(k); if isinf(Lk), Lk = median(Lstar(c,:),'omitnan'); end
        rows(end+1,:) = {cases(c).name, Ls(k), Lk, median(err(c,:,k),'omitnan'), median(dist(c,:,k),'omitnan'), ...
            median(loss(c,:,k),'omitnan'), median(secs(c,:,k),'omitnan'), median(errML(c,:),'omitnan'), median(secML(c,:),'omitnan')}; %#ok<AGROW>
    end
end
summary = cell2table(rows,'VariableNames',{'Case','L','LValue','MedianError','MedianDistanceToML', ...
    'MedianLossPerObs','MedianSeconds','MedianErrorML','MedianSecondsML'});
R = struct('cases',{{cases.name}},'Ls',Ls,'err',err,'dist',dist,'loss',loss,'secs',secs, ...
    'errML',errML,'secML',secML,'Lstar',Lstar,'summary',summary,'settings',o);
save(fullfile(outDir,'integer_projection_results.mat'),'R');
writetable(summary,fullfile(outDir,'integer_projection_summary.csv'));
disp(summary);
end

function [e,d,l,s,eM,sM,Lst] = oneReplicate(A,mu0,m,Ls,o,c,r)
nL = numel(Ls); e = nan(1,nL); d = nan(1,nL); l = nan(1,nL); s = nan(1,nL);
rng(o.Seed + 1000*c + r,'twister');
Y = poissrnd(repmat(mu0.',o.N,1)) * A.';
Lst = max(Y(:)) + 1;
start = max(fitCumulantNNLS(A,Y,min(m,4)),1e-3);       % common start for all fits
t = tic; [muML,iML] = fitFullLikelihood(A,Y,start); sM = toc(t);
eM = norm(muML-mu0)/norm(mu0); llML = iML.loglik;
for k = 1:nL
    Lk = Ls(k); if isinf(Lk), Lk = Lst; end
    cvec = Lk.^(0:m-1);
    try
        t = tic;
        mu = estimateMuOptimizedProjections(A,Y,'Display','off','ProjectionMatrix',cvec,'Mu0',start);
        s(k) = toc(t); mu = max(double(mu(:)),0);
        e(k) = norm(mu-mu0)/norm(mu0);
        d(k) = norm(mu-muML)/max(norm(muML-mu0),eps);
        [~,iE] = fitFullLikelihood(A,Y,mu,'EvaluateOnly',true);
        l(k) = (llML - iE.loglik)/o.N;
    catch ME
        warning('IntegerProjection:Failed','L=%g: %s',Lk,ME.message);
    end
end
end
