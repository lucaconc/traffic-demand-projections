function S = runMissingRoutesStudy(varargin)
%RUNMISSINGROUTESSTUDY  Estimation when the route dictionary misses active routes.
%
%   S = runMissingRoutesStudy('SyntheticResults',file, ...)
%
%   On the realistic network of the synthetic study, counts are generated from all
%   active classes (results.A, results.trueClassRates), but each estimator sees a
%   dictionary from which a fraction of the active classes has been removed.
%   Sensors left without any class are dropped. Reported, per method and
%   fraction: the median relative error on the classes that remain, and the median
%   share of the missing flow absorbed by them,
%       (sum(muHat) - sum(mu_kept)) / sum(mu_missing).
%
%   Options: 'SyntheticResults' (required), 'Fractions' ([0 .1 .2 .3]),
%   'Replicates' (20), 'N' (500), 'PriorCV' (0.5, for GLS-P), 'Seed',
%   'UseParallel' (true), 'OutputDir' ('missing_routes'), 'PaperResultsFile' (''),
%   'Methods' ({'PCL','GMM','GLS-P','NNLS'}).
ip = inputParser;
addParameter(ip,'SyntheticResults','');
addParameter(ip,'Fractions',[0 .1 .2 .3]);
addParameter(ip,'Replicates',20);
addParameter(ip,'N',500);
addParameter(ip,'PriorCV',0.5);
addParameter(ip,'Seed',20261016);
addParameter(ip,'UseParallel',true);
addParameter(ip,'OutputDir','missing_routes');
addParameter(ip,'PaperResultsFile','');
addParameter(ip,'Methods',{'PCL','GMM','GLS-P','NNLS'},@iscell);
parse(ip,varargin{:}); o = ip.Results;
assert(~isempty(o.SyntheticResults) && isfile(o.SyntheticResults), ...
    'runMissingRoutesStudy:SyntheticResults','Give the path of synthetic_results.mat (run the synthetic stage first).');
d = load(o.SyntheticResults,'results'); A = double(d.results.A ~= 0);
TR = d.results.trueClassRates; n = size(A,2);
R0 = o.Replicates; F = o.Fractions(:).'; nF = numel(F);
methods = o.Methods; nM = numel(methods);
% True rates of every replicate, as a matrix (sliced in parfor)
RATES = zeros(n,R0); for r = 1:R0, RATES(:,r) = TR(:,1+mod(r-1,size(TR,2))); end
err = nan(R0,nF,nM); absorbed = nan(R0,nF,nM);
fprintf('Missing-routes study: %d sensors, %d classes, %d replicates, N = %d, fractions %s\n', ...
    size(A,1),n,R0,o.N,mat2str(F));
useParfor = o.UseParallel && ~isempty(ver('parallel'));
if useParfor
    parfor r = 1:R0
        [err(r,:,:),absorbed(r,:,:)] = oneReplicate(A,RATES(:,r),F,methods,o,r); %#ok<PFBNS>
    end
else
    for r = 1:R0
        [err(r,:,:),absorbed(r,:,:)] = oneReplicate(A,RATES(:,r),F,methods,o,r);
    end
end
S = struct('methods',{methods},'fractions',F,'err',err,'absorbed',absorbed, ...
    'medianError',medianOverReplicates(err),'medianAbsorbed',medianOverReplicates(absorbed), ...
    'settings',o);
if ~isfolder(o.OutputDir), mkdir(o.OutputDir); end
save(fullfile(o.OutputDir,'missing_routes_results.mat'),'S');
fprintf('\nMedian relative error on the remaining classes\n%10s',''); fprintf('%10s',methods{:}); fprintf('\n');
for f = 1:nF, fprintf('%9.0f%%',100*F(f)); fprintf('%10.3f',S.medianError(f,:)); fprintf('\n'); end
if ~isempty(o.PaperResultsFile), writeMissingRoutesResults(S,o.PaperResultsFile); end   % paper macros (paperdata)
end

function [e,a] = oneReplicate(A,truth,F,methods,o,r)
% One replicate: the same counts for every fraction and method.
n = size(A,2); nF = numel(F); nM = numel(methods); e = nan(nF,nM); a = nan(nF,nM);
rs = RandStream('mt19937ar','Seed',o.Seed+97*r);
X = poissrnd(repmat(truth.',o.N,1)); Y = X*A.';
act = find(truth > 0); perm = act(randperm(rs,numel(act)));
for f = 1:nF
    nMiss = round(F(f)*numel(act)); miss = perm(1:nMiss);
    keep = true(n,1); keep(miss) = false;
    Ak = A(:,keep); rows = any(Ak,2); Ak = Ak(rows,:); Yk = Y(:,rows);
    tk = truth(keep); missFlow = sum(truth(miss));
    for k = 1:nM
        try
            switch methods{k}
                case 'PCL',   est = estimateMuOptimizedProjections(Ak,Yk,'Display','off');
                case 'GMM',   est = estimateMuEfficientGMM(Ak,Yk);
                case 'NNLS',  est = fitCumulantNNLS(Ak,Yk,1);
                case 'GLS-P'
                    rsP = RandStream('mt19937ar','Seed',o.Seed+97*r+5);
                    base = tk; pos = base > 0; if any(pos), base(~pos) = 0.1*mean(base(pos)); else, base(:) = 1; end
                    prior = base.*exp(o.PriorCV*randn(rsP,numel(tk),1) - o.PriorCV^2/2);
                    est = fitPriorGLS(Ak,Yk,prior,o.PriorCV);
            end
            est = double(est(:));
            e(f,k) = norm(est - tk)/max(norm(tk),eps);
            if missFlow > 0, a(f,k) = (sum(est) - sum(tk))/missFlow; end
        catch ME
            warning('runMissingRoutesStudy:Fit','%s, fraction %.2f: %s',methods{k},F(f),ME.message);
        end
    end
end
end

function M = medianOverReplicates(X)
% Median over the first dimension, ignoring NaN (failed fits); returns fractions x methods.
sz = size(X); M = nan(sz(2),sz(3));
for f = 1:sz(2)
    for k = 1:sz(3)
        v = X(:,f,k); v = v(isfinite(v));
        if ~isempty(v), M(f,k) = median(v); end
    end
end
end
