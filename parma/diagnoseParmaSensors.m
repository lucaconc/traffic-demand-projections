function D = diagnoseParmaSensors(G,measuredEdgeSets,sig,meanY,estMask,segmentKey,names,Ypool,varargin)
%DIAGNOSEPARMASENSORS Flag sensors with a likely wrong direction or faulty data.
%
%   D = diagnoseParmaSensors(G,measuredEdgeSets,sig,meanY,estMask,segmentKey,names,Ypool)
%
%   Uses the dictionary-wide NNLS fit of the pooled sensor means (all candidate
%   path classes, estimation sensors only) and tests the sensors with the
%   largest relative residuals (predicted - observed)/observed:
%
%   1. Leave-one-out: refit without the sensor and predict it. A sensor that
%      the other sensors cannot predict is inconsistent with the path set:
%      faulty detector, wrong road segment, or traffic outside the network.
%   2. Direction reversal: replace the sensor's row by the opposite directed
%      edges of the same road (fraction of each class's sampled member paths,
%      weighted by their hits, that use them) and refit. A clear drop of the
%      overall misfit suggests the direction is reversed.
%   3. Pair swap: for two sensors on the same road segment (same SegmentKey),
%      swap their observed means and refit. A clear drop suggests that the two
%      directions are exchanged.
%   4. Data quality, from the pooled complete intervals: share of zero counts
%      and longest run of identical nonzero counts (a stuck detector).
%
%   Name/value: TopK (20) sensors tested; RelImprovement (0.10) relative drop
%   of the overall misfit that counts as clear; LOOThreshold (0.5) relative
%   leave-one-out error; StuckRun (12) identical consecutive counts; ZeroShare
%   (0.3). Validation sensors are not fitted: for them the residual is already
%   a leave-one-out prediction, and the reversal test compares their prediction
%   errors.
%
%   D is a table, one row per measurement, with the test results, a Flag and
%   an Advice column.
%
%   Limits (checked on synthetic networks with planted faults): the tests rank
%   suspicious sensors well, but the suggested cause is only a hint. A faulty
%   detector can look like a reversed direction and vice versa, healthy
%   sensors coupled to a faulty one can be flagged as inconsistent, and a
%   sensor whose paths cross no other sensor cannot be checked at all
%   (Redundancy near 0). Inspect flagged sensors on the map before excluding
%   or correcting them.

p = inputParser;
addParameter(p,'TopK',20,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
addParameter(p,'RelImprovement',0.10,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'LOOThreshold',0.5,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'StuckRun',12,@(x)isnumeric(x)&&isscalar(x)&&x>=2);
addParameter(p,'ZeroShare',0.3,@(x)isnumeric(x)&&isscalar(x)&&x>0&&x<1);
parse(p,varargin{:}); o = p.Results;

m = numel(measuredEdgeSets);
Aall = double(sig.Aall);
meanY = double(meanY(:)); estMask = logical(estMask(:));
est = find(estMask);
opt = optimset('Display','off');
nnlsFit = @(Ar,yr) lsqnonneg(Ar,yr,opt);
mis = @(Ar,yr,mu) norm(Ar*mu-yr)/max(norm(yr),eps);

% Base fit on the estimation sensors
mu0 = nnlsFit(Aall(est,:),meanY(est));
f0 = mis(Aall(est,:),meanY(est),mu0);
% If the dictionary reproduces the means almost exactly, it can absorb any
% inconsistency, and relative misfit drops are meaningless: skip those tests.
canCompare = f0 >= 0.01;
pred0 = Aall*mu0;
res0 = (pred0-meanY)./max(meanY,1);   % denominator >= 1 vehicle: no blow-up for dead sensors

% Data quality on the pooled complete intervals
zeroShare = nan(m,1); longestRun = nan(m,1);
if ~isempty(Ypool) && size(Ypool,2) == m
    for i = 1:m
        y = double(Ypool(:,i)); y = y(isfinite(y));
        if isempty(y), continue; end
        zeroShare(i) = mean(y==0);
        longestRun(i) = longestConstantRun(y);
    end
end

% Candidates: largest |relative residual| among sensors with traffic
score = abs(res0); score(~(meanY>0) | ~isfinite(score)) = -Inf;
[~,ord] = sort(score,'descend');
cand = ord(1:min(o.TopK,nnz(isfinite(score) & score>-Inf)));

% Redundancy: share of the sensor's fitted flow on classes that also cross
% another estimation sensor (structural share if the fitted flow is zero). A
% sensor without shared classes is fitted exactly whatever it measures, so it
% cannot be checked against the others.
redundancy = nan(m,1);
for i = 1:m
    cross = Aall(i,:)>0;
    others = setdiff(est,i);
    shared = cross & any(Aall(others,:)>0,1);
    flow = Aall(i,:).*mu0(:).';
    if sum(flow) > 0
        redundancy(i) = sum(flow(shared))/sum(flow);
    else
        redundancy(i) = nnz(shared)/max(nnz(cross),1);
    end
end
looErr = nan(m,1); dropRev = nan(m,1); revAvail = false(m,1);
pairWith = nan(m,1); dropSwap = nan(m,1);
for i = cand(:).'
    % 1. Leave-one-out (estimation sensors); validation sensors are already out
    if estMask(i)
        rows = setdiff(est,i);
        mu = nnlsFit(Aall(rows,:),meanY(rows));
        looErr(i) = (Aall(i,:)*mu-meanY(i))/max(meanY(i),1);
    else
        looErr(i) = res0(i);
    end
    % 2. Direction reversal
    rowRev = reversedRow(G,measuredEdgeSets{i},sig);
    if ~isempty(rowRev)
        revAvail(i) = true;
        Ar = Aall; Ar(i,:) = rowRev;
        if estMask(i) && canCompare
            mu = nnlsFit(Ar(est,:),meanY(est));
            dropRev(i) = (f0-mis(Ar(est,:),meanY(est),mu))/f0;
        elseif ~estMask(i)
            e0 = abs(pred0(i)-meanY(i)); e1 = abs(Ar(i,:)*mu0-meanY(i));
            dropRev(i) = (e0-e1)/max(e0,eps);
        end
    end
    % 3. Pair swap with the other direction of the same road segment
    if ~isempty(segmentKey) && strlength(string(segmentKey(i)))>0
        j = find(string(segmentKey)==string(segmentKey(i)) & estMask);
        j = setdiff(j,i);
        if estMask(i) && ~isempty(j) && canCompare
            j = j(1); pairWith(i) = j;
            ys = meanY; ys([i j]) = meanY([j i]);
            mu = nnlsFit(Aall(est,:),ys(est));
            dropSwap(i) = (f0-mis(Aall(est,:),ys(est),mu))/f0;
        end
    end
end

% Flags
flag = strings(m,1); advice = strings(m,1);
for i = 1:m
    parts = strings(0,1); adv = strings(0,1);
    if dropSwap(i) >= o.RelImprovement
        parts(end+1) = "directions swapped?"; %#ok<AGROW>
        adv(end+1) = sprintf("swapping with measurement %d reduces the misfit by %.0f%%",pairWith(i),100*dropSwap(i)); %#ok<AGROW>
    end
    if dropRev(i) >= o.RelImprovement
        parts(end+1) = "direction reversed?"; %#ok<AGROW>
        adv(end+1) = sprintf("the opposite direction reduces the misfit by %.0f%%",100*dropRev(i)); %#ok<AGROW>
    end
    if isempty(parts) && ismember(i,cand) && redundancy(i) < 0.2
        parts(end+1) = "not checkable"; %#ok<AGROW>
        adv(end+1) = "its paths cross no other sensor, so the counts cannot be compared"; %#ok<AGROW>
    elseif isempty(parts) && abs(looErr(i)) >= o.LOOThreshold
        parts(end+1) = "inconsistent with other sensors"; %#ok<AGROW>
        adv(end+1) = sprintf("predicted from the others with relative error %+.2f: check detector, matching, direction, boundary",looErr(i)); %#ok<AGROW>
    end
    if zeroShare(i) >= o.ZeroShare
        parts(end+1) = "many zeros"; %#ok<AGROW>
        adv(end+1) = sprintf("%.0f%% of counts are zero",100*zeroShare(i)); %#ok<AGROW>
    end
    if longestRun(i) >= o.StuckRun
        parts(end+1) = "stuck values"; %#ok<AGROW>
        adv(end+1) = sprintf("%d identical consecutive nonzero counts",longestRun(i)); %#ok<AGROW>
    end
    flag(i) = strjoin(parts,"; "); advice(i) = strjoin(adv,"; ");
end
role = repmat("estimation",m,1); role(~estMask) = "validation";
if isempty(names), names = strings(m,1); end
D = table((1:m).',string(names(:)),role,meanY,pred0,res0,redundancy,looErr,revAvail,dropRev,pairWith,dropSwap, ...
    zeroShare,longestRun,flag,advice, ...
    'VariableNames',{'Measurement','SectionName','Role','MeanCount','PredictedNNLS','RelResidual', ...
    'Redundancy','LeaveOneOutError','ReverseAvailable','MisfitDropIfReversed','PairedWith','MisfitDropIfSwapped', ...
    'ZeroShare','LongestConstantRun','Flag','Advice'});
D.Properties.UserData = struct('baseMisfit',f0,'tested',cand(:));
end

% -------------------------------------------------------------------------
function row = reversedRow(G,edgeIds,sig)
% Class coefficients for the opposite directed edges of the sensor's road:
% hit-weighted share of each class's sampled member paths that use them.
row = [];
ids = double(edgeIds(:)); ids = ids(isfinite(ids) & ids>=1 & ids<=numedges(G));
if isempty(ids), return; end
en = G.Edges.EndNodes(ids,:);
if ~isnumeric(en), en = [findnode(G,en(:,1)) findnode(G,en(:,2))]; end
rev = findedge(G,en(:,2),en(:,1)); rev = rev(rev>0);
if isempty(rev), return; end
Q = size(sig.Aall,2); row = zeros(1,Q);
for q = 1:Q
    paths = sig.memberEdgePaths{q};
    if isempty(paths), continue; end
    if ~iscell(paths), paths = {paths}; end
    hits = ones(numel(paths),1);
    if isfield(sig,'memberPathHits') && numel(sig.memberPathHits)>=q && numel(sig.memberPathHits{q})==numel(paths)
        hits = double(sig.memberPathHits{q}(:));
    end
    use = cellfun(@(pe)any(ismember(pe,rev)),paths(:));
    row(q) = sum(hits(use))/max(sum(hits),eps);
end
end

function r = longestConstantRun(y)
% Longest run of identical consecutive nonzero values.
r = 0; cur = 1;
for k = 2:numel(y)
    if y(k)==y(k-1) && y(k)>0, cur = cur+1; else, cur = 1; end
    r = max(r,cur);
end
if numel(y)==1 && y(1)>0, r = 1; end
end
