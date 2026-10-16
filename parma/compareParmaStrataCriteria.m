function T = compareParmaStrataCriteria(compactFile,varargin)
%COMPAREPARMASTRATACRITERIA Compare the BIC and cross-validated stratifications.
%
%   T = compareParmaStrataCriteria(compactFile,<options of runParmaTrafficTomography>)
%
%   Runs the Parma study twice with identical options, StrataCriterion 'bic'
%   and 'cv', saving parma_strata_bic.mat and parma_strata_cv.mat. Both runs
%   use the same test days (fixed before stratification), and the cv run
%   selects its strata on separate validation days, so the test errors are
%   directly comparable. T, also written to parma_strata_comparison.csv, has:
%     Strata, ReliableStrata
%     PooledTestNRMSE        NRMSE of all sensor-stratum mean counts on the test days
%                            (reliable strata), the overall prediction error
%     MedianTestNRMSE        median over strata of the per-stratum test NRMSE
%     MedianTestLocalErrorPct, MedianTestLocalErrorNNLSPct
%                            median |local relative error| (pairs >= 30 veh/h),
%                            composite estimate and NNLS baseline
%     MedianNNLSMisfit, MedianTrainMisfit, Seconds

crit = ["bic","cv"]; rows = cell(2,1);
for k = 1:2
    fprintf('\n===== Strata criterion: %s =====\n',crit(k));
    t0 = tic;
    try
        R = runParmaTrafficTomography(compactFile,varargin{:},'StrataCriterion',char(crit(k)), ...
            'Plot',false,'SaveFile',sprintf('parma_strata_%s.mat',crit(k)));
        rows{k} = summaryRow(R,crit(k),toc(t0),"");
    catch ME
        fprintf(2,'Criterion %s failed: %s\n',crit(k),ME.message);
        rows{k} = summaryRow([],crit(k),toc(t0),string(ME.message));
    end
end
T = vertcat(rows{:});
disp(T);
writetable(T,safeOutputPath('parma_strata_comparison.csv'));
fprintf('Comparison written to parma_strata_comparison.csv\n');
end

function r = summaryRow(R,label,secs,failure)
names = {'Criterion','Strata','ReliableStrata','PooledTestNRMSE','MedianTestNRMSE', ...
    'MedianTestLocalErrorPct','MedianTestLocalErrorNNLSPct','MedianNNLSMisfit','MedianTrainMisfit', ...
    'Seconds','Failure'};
if isempty(R)
    r = table(label,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,secs,failure,'VariableNames',names); return
end
valid = R.validFitByStratum(:);
med = @(x) median(x(isfinite(x)));
meas = R.testMeasuredMeanByStratum; pred = R.predictedMeanByStratum;
predN = nan(size(pred)); if isfield(R,'predictedMeanNNLSByStratum'), predN = R.predictedMeanNNLSByStratum; end
use = logical(R.sensorUsedByStratum) & repmat(valid.',size(meas,1),1) & isfinite(meas) & isfinite(pred);
pooled = sqrt(mean((pred(use)-meas(use)).^2))/max(sqrt(mean(meas(use).^2)),eps);
minObs = repmat(30*[R.strata.ObservationMinutes]/60,size(meas,1),1);
ok = use & meas>=minObs & meas>0;
loc = abs(100*(pred-meas)./meas); locN = abs(100*(predN-meas)./meas);
nrmse = NaN;
if isfield(R,'strataTable') && ismember('TestNRMSE',R.strataTable.Properties.VariableNames)
    nrmse = med(R.strataTable.TestNRMSE(valid));
end
r = table(label,numel(valid),nnz(valid),pooled,nrmse,med(loc(ok)),med(locN(ok)), ...
    med(R.nnlsMisfitByStratum(valid)),med(R.trainMisfitByStratum(valid)),secs,failure, ...
    'VariableNames',names);
end
