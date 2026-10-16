function [T,P] = compareParmaEstimators(compactFile,varargin)
%COMPAREPARMAESTIMATORS Parma study with the PCL, the GMM and the robust GMM.
%
%   [T,P] = compareParmaEstimators(compactFile,<options of runParmaTrafficTomography>)
%   [T,P] = compareParmaEstimators(compactFile,...,'OutputDir',dir,'MacroFile',file)
%
%   OutputDir (current folder) receives the three result files and the CSVs;
%   MacroFile (OutputDir/paper_results_parma_estimators.tex) the paper macros
%   (writeParmaEstimatorResults). Both lie outside the library (safeOutputPath).
%
%   Runs the study three times with identical options and seeds, so that the
%   strata, test days, class sets and projection designs coincide, changing only
%   StratumEstimator ('pcl','gmm','robustGMM'). Results are saved as
%   parma_estimator_<name>.mat. Outputs (also written as CSV):
%
%   T  one row per estimator:
%      ReliableStrata             strata passing all reliability checks
%      PooledTestNRMSE            NRMSE of all sensor-stratum mean counts on the test
%                                 days (reliable strata)
%      MedianTestNRMSE            median of the per-stratum test NRMSE
%      MedianTestLocalErrorPct    median |local relative error| on the test days
%                                 (pairs with >= 30 veh/h, strata with >= 2 test days)
%      MedianTrainMisfit          median relative misfit of the training means
%      GMMFailures                strata in which the GMM failed (PCL estimate kept)
%      MedianSecondsPerStratum, TotalSeconds
%   P  paired comparison on the strata reliable under all three estimators: for
%      the GMM and the robust GMM, the share of strata whose test NRMSE is lower
%      than with the PCL, and the median ratio of test NRMSE to the PCL's.
ip = inputParser; ip.KeepUnmatched = true;
addParameter(ip,'OutputDir','.'); addParameter(ip,'PaperText',true,@(x)islogical(x)&&isscalar(x));   % false: no paper macros
addParameter(ip,'MacroFile',''); addParameter(ip,'IncludeRobust',false);
addParameter(ip,'PCLResultFile','');   % reuse the parma stage's PCL results when their settings match
parse(ip,varargin{:});
outDir = char(safeOutputPath(ip.Results.OutputDir)); if ~isfolder(outDir), mkdir(outDir); end
macroFile = ip.Results.MacroFile;
if isempty(macroFile), macroFile = fullfile(outDir,'paper_results_parma_estimators.tex'); end
u = ip.Unmatched; fn = fieldnames(u); varargin = cell(1,2*numel(fn));
for k = 1:numel(fn), varargin{2*k-1} = fn{k}; varargin{2*k} = u.(fn{k}); end
est = {'pcl','gmm'}; if ip.Results.IncludeRobust, est{end+1} = 'robustGMM'; end
nE = numel(est); R = cell(1,nE); secs = zeros(1,nE);
for k = 1:nE
    if strcmp(est{k},'pcl')
        [Rk,why] = reusablePCL(ip.Results.PCLResultFile,varargin);
        if ~isempty(Rk)
            fprintf('\n===== Parma with StratumEstimator = pcl: reusing %s =====\n',ip.Results.PCLResultFile);
            R{k} = Rk; secs(k) = NaN; continue
        elseif ~isempty(ip.Results.PCLResultFile)
            fprintf('\nPCL results not reused (%s); refitting.\n',why);
        end
    end
    fprintf('\n===== Parma with StratumEstimator = %s =====\n',est{k});
    t0 = tic;
    R{k} = runParmaTrafficTomography(compactFile,varargin{:},'StratumEstimator',est{k}, ...
        'Plot',false,'SaveFile',fullfile(outDir,sprintf('parma_estimator_%s.mat',est{k})));
    secs(k) = toc(t0);
end
rows = cell(nE,10);
for k = 1:nE, rows(k,:) = summaryRow(R{k},est{k},secs(k)); end
T = cell2table(rows,'VariableNames',{'Estimator','Strata','ReliableStrata','PooledTestNRMSE', ...
    'MedianTestNRMSE','MedianTestLocalErrorPct','MedianTrainMisfit','GMMFailures', ...
    'MedianSecondsPerStratum','TotalSeconds'});
% Paired comparison on the strata reliable under all three
common = true(size(R{1}.validFitByStratum(:)));
for k = 1:nE, common = common & R{k}.validFitByStratum(:); end
nr = cellfun(@(r) stratumNRMSE(r),R,'UniformOutput',false);
prow = cell(nE-1,5);
for k = 2:nE
    a = nr{k}(common); b = nr{1}(common); ok = isfinite(a) & isfinite(b) & b > 0;
    prow(k-1,:) = {est{k},nnz(ok),mean(a(ok) < b(ok)),median(a(ok)./b(ok)),median(a(ok)-b(ok))};
end
P = cell2table(prow,'VariableNames',{'Estimator','CommonStrata','ShareBetterThanPCL', ...
    'MedianNRMSERatioToPCL','MedianNRMSEDifference'});
disp(T); disp(P);
writetable(T,fullfile(outDir,'parma_estimator_comparison.csv'));
writetable(P,fullfile(outDir,'parma_estimator_paired.csv'));
if nE >= 3, Rrob = R{3}; else, Rrob = []; end
if ip.Results.PaperText, writeParmaEstimatorResults(R{1},R{2},Rrob,macroFile); end
end

function row = summaryRow(R,name,total)
valid = R.validFitByStratum(:);
med = @(x) median(x(isfinite(x)));
meas = R.testMeasuredMeanByStratum; pred = R.predictedMeanByStratum;
use = logical(R.sensorUsedByStratum) & repmat(valid.',size(meas,1),1) & isfinite(meas) & isfinite(pred);
if isfield(R,'testDaysByStratum'), use = use & repmat(R.testDaysByStratum(:).' >= 2,size(meas,1),1); end
pooled = sqrt(mean((pred(use)-meas(use)).^2))/max(sqrt(mean(meas(use).^2)),eps);
minObs = repmat(30*[R.strata.ObservationMinutes]/60,size(meas,1),1);
ok = use & meas >= minObs & meas > 0;
loc = abs(100*(pred-meas)./meas);
fails = 0;
for q = 1:numel(R.infoByStratum)
    iq = R.infoByStratum{q};
    if isstruct(iq) && isfield(iq,'gmmError'), fails = fails + 1; end
end
row = {string(name),numel(valid),nnz(R.reliableByStratum),pooled,med(stratumNRMSE(R)), ...
    med(loc(ok)),med(R.trainMisfitByStratum(valid)),fails,med(R.estimationTimeByStratum(valid)),total};
end

function v = stratumNRMSE(R)
% Per-stratum test NRMSE of the sensor mean counts.
S = numel(R.validFitByStratum); v = nan(S,1);
for q = 1:S
    if ~R.validFitByStratum(q), continue; end
    m = R.testMeasuredMeanByStratum(:,q); p = R.predictedMeanByStratum(:,q);
    u = logical(R.sensorUsedByStratum(:,q)) & isfinite(m) & isfinite(p);
    if any(u), v(q) = sqrt(mean((p(u)-m(u)).^2))/max(sqrt(mean(m(u).^2)),eps); end
end
end

function [R,why] = reusablePCL(file,opts)
% Load the PCL results of the parma stage if they were produced with the same
% result-relevant options as this run (output paths, caches, plotting and
% parallelism are ignored); otherwise return [] and the reason.
R = []; why = '';
if isempty(file), why = 'no file given'; return; end
if ~isfile(file), why = 'file not found'; return; end
d = load(file); fn = fieldnames(d);
if isfield(d,'results'), S = d.results; elseif isfield(d,'R'), S = d.R; else, S = d.(fn{1}); end
if ~isfield(S,'settings'), why = 'no saved settings'; return; end
st = S.settings;
if isfield(st,'StratumEstimator') && ~strcmpi(string(st.StratumEstimator),'pcl')
    why = 'saved results are not from the PCL'; return
end
ignore = {'SaveFile','CacheDir','OfficialDataCacheDir','Plot','UseParallel','AnimationVideoFile'};
for k = 1:2:numel(opts)
    name = char(opts{k});
    isOutputCSV = endsWith(name,'CSV') && ~strcmpi(name,'EdgeOverrideCSV');   % EdgeOverrideCSV is an input
    if any(strcmpi(name,ignore)) || isOutputCSV, continue; end
    if ~isfield(st,name), why = sprintf('option %s not in the saved settings',name); return; end
    a = st.(name); b = opts{k+1};
    if ischar(a) || isstring(a), same = strcmpi(string(a),string(b)); else, same = isequaln(a,b); end
    if ~same, why = sprintf('option %s differs',name); return; end
end
R = S;
end
