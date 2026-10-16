function [T,excluded] = optimizeParmaSensorSet(compactFile,varargin)
%OPTIMIZEPARMASENSORSET Exclude inconsistent sensors and check the gain.
%
%   [T,excluded] = optimizeParmaSensorSet(compactFile,<options of runParmaTrafficTomography>)
%   [T,excluded] = optimizeParmaSensorSet(compactFile,'BaselineFile','parma_tomography_results.mat',...)
%
%   1. Baseline: run the Parma study (or load BaselineFile, a saved run of the
%      modified pipeline). Its sensor selection (selectConsistentSensors) lists
%      the sensors whose deletion most improves the NNLS fit of the others.
%   2. Rerun with those sensors excluded (ExcludeSectionNames).
%   3. Compare the two runs on the temporal test (held-out days) and print T,
%      also written to parma_sensorset_comparison.csv.
%
%   Use the excluded names in the build with
%     cfg.parma.options = {'ExcludeSectionNames',excluded};
%   The comparison uses the same strata definition only approximately: removing
%   sensors changes the summed counts, and hence the time windows found.

ip = inputParser; ip.KeepUnmatched = true;
addParameter(ip,'BaselineFile','',@(x)ischar(x)||isstring(x));
parse(ip,varargin{:});
rest = unmatchedPairs(ip.Unmatched);
if strlength(string(ip.Results.BaselineFile)) > 0
    d = load(char(ip.Results.BaselineFile),'results'); R0 = d.results;
else
    R0 = runParmaTrafficTomography(compactFile,rest{:},'Plot',false,'SaveFile','parma_sensorset_baseline.mat');
end
assert(isfield(R0,'sensorSelection') && isfield(R0.sensorSelection,'recommendedNames'), ...
    'optimizeParmaSensorSet:NoSelection','The baseline has no sensor selection: rerun the modified pipeline.');
excluded = R0.sensorSelection.recommendedNames(:);
fprintf('\nExcluding %d sensors: %s\n',numel(excluded),strjoin(excluded.','; '));
if isempty(excluded)
    T = table(); fprintf('Nothing to exclude.\n'); return
end
R1 = runParmaTrafficTomography(compactFile,rest{:},'ExcludeSectionNames',excluded, ...
    'Plot',false,'SaveFile','parma_sensorset_excluded.mat');
T = [summaryRow(R0,"all sensors",0); summaryRow(R1,"recommended excluded",numel(excluded))];
disp(T);
writetable(T,safeOutputPath('parma_sensorset_comparison.csv'));
end

function r = summaryRow(R,label,nExcluded)
valid = R.validFitByStratum(:);
med = @(x) median(x(isfinite(x)));
nrmse = NaN;
if isfield(R,'strataTable') && ismember('TestNRMSE',R.strataTable.Properties.VariableNames)
    nrmse = med(R.strataTable.TestNRMSE(valid));
end
meas = R.testMeasuredMeanByStratum; pred = R.predictedMeanByStratum;
minObs = repmat(30*[R.strata.ObservationMinutes]/60,size(meas,1),1);
ok = logical(R.sensorUsedByStratum) & repmat(valid.',size(meas,1),1) & isfinite(meas) & meas>=minObs & meas>0;
loc = abs(100*(pred-meas)./meas);
r = table(string(label),nExcluded,numel(valid),nnz(valid),nrmse,med(loc(ok)), ...
    med(R.nnlsMisfitByStratum(valid)),med(R.trainMisfitByStratum(valid)), ...
    'VariableNames',{'Run','ExcludedSensors','Strata','ReliableStrata','MedianTestNRMSE', ...
    'MedianTestLocalErrorPct','MedianNNLSMisfit','MedianTrainMisfit'});
end

function c = unmatchedPairs(u)
f = fieldnames(u); c = cell(1,2*numel(f));
for k = 1:numel(f), c{2*k-1} = f{k}; c{2*k} = u.(f{k}); end
end
