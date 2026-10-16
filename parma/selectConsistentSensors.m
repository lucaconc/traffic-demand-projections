function S = selectConsistentSensors(A,meanY,names,varargin)
%SELECTCONSISTENTSENSORS Greedy removal of sensors inconsistent with the others.
%
%   S = selectConsistentSensors(A,meanY,names)
%   S = selectConsistentSensors(A,meanY,names,'MaxRemoveShare',0.2,'MinRelGain',0.02)
%
%   A       sensors x classes signature matrix of the candidate path dictionary
%           (estimation sensors only), meanY the pooled mean count of each sensor,
%           names a cell array (or string array) of sensor names.
%
%   The NRMSE of the NNLS fit of the sensor means, sqrt(mean(r.^2))/sqrt(mean(y.^2))
%   over the sensors kept, measures how consistent they are with the path
%   dictionary. At each step every candidate is deleted in turn and the sensor
%   whose deletion lowers the NRMSE of the OTHER sensors the most is removed: a
%   faulty or mismatched sensor distorts the fit of all the others, a healthy
%   one barely affects it. The search stops when no deletion lowers the NRMSE
%   by at least MinRelGain (5%, relative) and MinAbsGain (0.02), when MaxRemoveShare of the sensors has
%   been removed, or when fewer than MinSensors would remain.
%
%   Sensors whose classes cross no other sensor (Redundancy < MinRedundancy)
%   are never removed and do not enter the NRMSE: they cannot be predicted from
%   the others, so a large error there says nothing about their data.
%
%   S.steps        table: Step, Removed, Name, RelResidual (before removal), NRMSEAfter
%   S.baseNRMSE    NRMSE of the NNLS fit before any removal
%   S.recommended  indices (into the rows of A) of the sensors to exclude
%   S.recommendedNames, S.finalNRMSE

p = inputParser;
addParameter(p,'LOOThreshold',0.3,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'MaxRemoveShare',0.2,@(x)isnumeric(x)&&isscalar(x)&&x>0&&x<1);
addParameter(p,'MinSensors',10,@(x)isnumeric(x)&&isscalar(x)&&x>=2);
addParameter(p,'MinAbsGain',0.02,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
addParameter(p,'MinRelGain',0.05,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
addParameter(p,'MinRedundancy',0.2,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
parse(p,varargin{:}); o = p.Results;

A = double(A); y = double(meanY(:)); m = size(A,1);
if isstring(names), names = cellstr(names); end
if isempty(names), names = arrayfun(@(k)sprintf('sensor %d',k),(1:m).','UniformOutput',false); end
opt = optimset('Display','off');
active = y > 0;
% Structural redundancy (does any class through the sensor cross another sensor?)
checkable = false(m,1);
for i = 1:m
    through = A(i,:) > 0;
    others = active; others(i) = false;
    shared = through & any(A(others,:) > 0,1);
    checkable(i) = nnz(through) > 0 && nnz(shared)/nnz(through) >= o.MinRedundancy;
end

maxRemove = floor(o.MaxRemoveShare*nnz(active));
removed = false(m,1);
nrmse0 = fitMisfit(A,y,active&~removed,opt);
stepRemoved = zeros(0,1); stepErr = zeros(0,1); stepNRMSE = zeros(0,1);
cur = nrmse0;
while nnz(removed) < maxRemove && nnz(active&~removed) > o.MinSensors
    cand = find(active & ~removed & checkable);
    if isempty(cand), break; end
    % Misfit of the remaining sensors after deleting each candidate
    f = inf(numel(cand),1); relRes = nan(numel(cand),1);
    use = active & ~removed;
    mu = lsqnonneg(A(use,:),y(use),opt); res = (A*mu-y)./max(y,1);
    for k = 1:numel(cand)
        u = use; u(cand(k)) = false;
        f(k) = fitMisfit(A,y,u,opt);
        relRes(k) = res(cand(k));
    end
    [fbest,k] = min(f);
    % Stop unless the deletion clearly improves the others, in relative AND absolute terms
    if fbest > (1-o.MinRelGain)*cur || cur-fbest < o.MinAbsGain, break; end
    i = cand(k);
    removed(i) = true; cur = fbest;
    stepRemoved(end+1,1) = i; stepErr(end+1,1) = relRes(k); stepNRMSE(end+1,1) = fbest; %#ok<AGROW>
end
nRec = numel(stepRemoved);
S = struct();
S.steps = table((1:numel(stepRemoved)).',stepRemoved,names(stepRemoved),stepErr,stepNRMSE, ...
    'VariableNames',{'Step','Removed','Name','RelResidual','NRMSEAfter'});
S.baseNRMSE = nrmse0;
S.recommended = stepRemoved(1:nRec);
S.recommendedNames = names(S.recommended);
if nRec > 0, S.finalNRMSE = stepNRMSE(nRec); else, S.finalNRMSE = nrmse0; end
S.checkable = checkable;
end

function f = fitMisfit(A,y,use,opt)
% NRMSE of the NNLS fit of the means of the sensors in USE (in-sample).
mu = lsqnonneg(A(use,:),y(use),opt);
r = A(use,:)*mu-y(use);
f = sqrt(mean(r.^2))/max(sqrt(mean(y(use).^2)),eps);
end

function [nrmse,err] = looErrors(A,y,use,checkable,opt)
% Leave-one-out predictions of the sensors in USE, from the other sensors in USE.
m = size(A,1); err = nan(m,1); pred = nan(m,1);
idx = find(use);
for i = idx(:).'
    rows = idx(idx~=i);
    mu = lsqnonneg(A(rows,:),y(rows),opt);
    pred(i) = A(i,:)*mu;
    err(i) = (pred(i)-y(i))/max(y(i),1);
end
ev = use & checkable & isfinite(pred);
nrmse = sqrt(mean((pred(ev)-y(ev)).^2))/max(sqrt(mean(y(ev).^2)),eps);
end
