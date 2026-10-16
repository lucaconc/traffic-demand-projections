function T = compareParmaAggregation(compactFile,minutesList,varargin)
%COMPAREPARMAAGGREGATION Run the Parma study at several time aggregations.
%
%   T = compareParmaAggregation(compactFile)                 % 10, 20 and 30 minutes
%   T = compareParmaAggregation(compactFile,[10 20 30 60],<options of runParmaTrafficTomography>)
%
%   Each run uses the same options except AggregationMinutes, saves its results
%   to parma_aggregation_<k>min.mat and writes no plots. The table T, also
%   written to parma_aggregation_comparison.csv and printed, compares:
%     Strata, ReliableStrata     strata built and strata with a reliable estimate
%     MedianClassesFitted        classes fitted per stratum
%     MedianNNLSMisfit           training-mean misfit of the NNLS fit (path-set adequacy)
%     MedianTrainMisfit          same for the composite estimate
%     MedianTestLocalError       median |local relative error| on the held-out days, in %
%                                (pairs with >= 30 vehicles/hour), composite estimate
%     MedianTestLocalErrorNNLS   same for the NNLS baseline
%     MedianTestNRMSE            per-stratum NRMSE of the held-out-day sensor means
%     MedianFano                 median sensor Fano factor in the strata
%     DownweightedShare          median share of projected counts with weight < 0.5
%     Seconds                    run time
%   Errors and misfits are dimensionless, so they are comparable across aggregations.
%   Use the same seed and options for all runs (the default), so that only the
%   aggregation changes.

if nargin<2 || isempty(minutesList), minutesList=[10 20 30]; end
rows=cell(numel(minutesList),1);
for k=1:numel(minutesList)
    mins=minutesList(k);
    fprintf('\n===== Aggregation %g minutes =====\n',mins);
    t0=tic;
    agg={};
    if mins>0, agg={'AggregationMinutes',mins}; end
    saveFile=sprintf('parma_aggregation_%gmin.mat',mins);
    try
        R=runParmaTrafficTomography(compactFile,varargin{:},agg{:},'Plot',false,'SaveFile',saveFile);
        rows{k}=summaryRow(R,mins,toc(t0),"");
    catch ME
        fprintf(2,'Aggregation %g min failed: %s\n',mins,ME.message);
        rows{k}=summaryRow([],mins,toc(t0),string(ME.message));
    end
end
T=vertcat(rows{:});
disp(T);
writetable(T,safeOutputPath('parma_aggregation_comparison.csv'));
fprintf('Comparison written to parma_aggregation_comparison.csv\n');
end

function r=summaryRow(R,mins,secs,failure)
nanRow=@() table(mins,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,secs,failure, ...
    'VariableNames',names());
if isempty(R), r=nanRow(); return; end
valid=R.validFitByStratum(:);
S=numel(valid);
med=@(x) median(x(isfinite(x)));
[e,eN]=testLocalErrors(R,valid);
nrmse=NaN;
if isfield(R,'strataTable') && ismember('TestNRMSE',R.strataTable.Properties.VariableNames)
    nrmse=med(R.strataTable.TestNRMSE(valid));
end
dw=NaN; if isfield(R,'downweightedFractionByStratum'), dw=med(R.downweightedFractionByStratum(valid)); end
r=table(mins,S,nnz(valid),med(R.reducedClassCountFitted(valid)), ...
    med(R.nnlsMisfitByStratum(valid)),med(R.trainMisfitByStratum(valid)),e,eN,nrmse, ...
    med(R.medianFanoByStratum(valid)),dw,secs,failure,'VariableNames',names());
end

function n=names()
n={'AggregationMinutes','Strata','ReliableStrata','MedianClassesFitted','MedianNNLSMisfit', ...
   'MedianTrainMisfit','MedianTestLocalError','MedianTestLocalErrorNNLS','MedianTestNRMSE', ...
   'MedianFano','DownweightedShare','Seconds','Failure'};
end

function [e,eN]=testLocalErrors(R,valid)
% Median |100*(predicted-observed)/observed| on the held-out days, over
% sensor-stratum pairs of reliable strata observing >= 30 vehicles/hour.
meas=R.testMeasuredMeanByStratum; pred=R.predictedMeanByStratum;
predN=nan(size(pred)); if isfield(R,'predictedMeanNNLSByStratum'), predN=R.predictedMeanNNLSByStratum; end
minutes=[R.strata.ObservationMinutes];
minObs=repmat(30*minutes/60,size(meas,1),1);
ok=logical(R.sensorUsedByStratum) & repmat(valid(:).',size(meas,1),1) & isfinite(meas) & meas>=minObs & meas>0;
loc=abs(100*(pred-meas)./meas); locN=abs(100*(predN-meas)./meas);
e=median(loc(ok & isfinite(loc))); eN=median(locN(ok & isfinite(locN)));
if isempty(e), e=NaN; end
if isempty(eN), eN=NaN; end
end
