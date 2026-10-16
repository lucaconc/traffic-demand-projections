function report = buildPaperArtifacts(cfg,stages)
%BUILDPAPERARTIFACTS Generate every figure and data file read by the paper.
%
%   setupPaperLib
%   cfg = paperConfig('/path/to/paper');          % folder of paper_scl.tex
%   cfg.parma.compactFile = '/path/to/parma_archive';
%   report = buildPaperArtifacts(cfg);             % all stages
%   report = buildPaperArtifacts(cfg,{'figures','captions','check'});
%
%   Stages (run in the order given):
%     'synthetic'  Monte Carlo study on the OSM network (Section 6). Saves the
%                  trials and writes paper_results_synthetic.tex.
%     'parma'      Parma case study (Section 7): main run and, if
%                  cfg.parma.screeningCap is set (default []: not), a second run with a larger
%                  class cap. Writes paper_results_parma.tex.
%     'macros'     Rewrite both result-macro files from the saved MAT files,
%                  without re-running the studies.
%     'figures'    Synthetic figures (Fig. 2), Parma figures (Figs. 3-5) and the
%                  computation-time figure (Fig. 6) from the saved MAT files.
%     'missing'    Estimation with missing routes (after 'synthetic').
%     'captions'   Caption data and caption commands from the saved MAT files.
%     'check'      List the files read by the paper and write a manifest.
%
%   All paths are relative to cfg.outputRoot; the current folder is
%   restored on exit.  A log of the command window is written to
%   cfg.paths.logs.  report.stages lists status, time and messages.

if nargin < 2 || isempty(stages), stages = cfg.stages; end
if ischar(stages) || isstring(stages), stages = cellstr(stages); end
known = {'synthetic','missing','benchmark','city','mle','integer','costsel','parma','estimators','macros','figures','captions','check'};
bad = setdiff(stages,known);
assert(isempty(bad),'buildPaperArtifacts:Stage','Unknown stage(s): %s.',strjoin(bad,', '));
libRoot = fileparts(mfilename('fullpath'));
addpath(fullfile(libRoot,'util'));                    % safeOutputPath is needed before path setup
cfg.outputRoot = char(safeOutputPath(cfg.outputRoot));   % never inside the library
if ~isfolder(cfg.outputRoot), mkdir(cfg.outputRoot); end
oldDir = pwd; oldPath = path;
restoreState = onCleanup(@()restoreSession(oldDir,oldPath));
% Use only this library: MathWorks default path + library folders, and a
% current folder without code (MATLAB searches the current folder first).
if ~isfield(cfg,'isolatePath') || cfg.isolatePath
    restoredefaultpath;
end
addpath(fullfile(libRoot,'util'));      % safeOutputPath (used below) must survive the path reset
% Settings missing from cfg (e.g. a cfg created by an older paperConfig)
% take their default values, so old configurations keep working.
cfg = withDefaults(cfg,paperConfigDefaults(cfg,libRoot));
P = absolutePaths(cfg);
ensureFolder(P.work);
cd(P.work);                 % first, so that nothing is resolved from the old folder
addpath(libRoot);
rehash;
setupPaperLib();
code = [dir(fullfile(P.work,'*.m')); dir(fullfile(P.work,'*.mlx')); dir(fullfile(P.work,'*.p'))];
assert(isempty(code),'buildPaperArtifacts:CodeInWork', ...
    'The work folder %s must not contain code files.',P.work);
conflicts = checkLibraryShadowing();
if ~isempty(conflicts)
    error('buildPaperArtifacts:Shadowed', ['These files shadow the library and would be ' ...
        'used instead of it; remove or rename them:\n  %s'],strjoin(conflicts,'\n  '));
end
fprintf('Library: %s\n',libRoot);
ensureFolder(P.captionData); ensureFolder(P.logs);
logFile = fullfile(P.logs,['build_' datestr(now,'yyyymmdd_HHMMSS') '.log']); %#ok<TNOW1,DATST>
cfg.paths = P;            % stages below use absolute paths
diary(logFile); stopDiary = onCleanup(@()diary('off'));
fprintf('Building paper artifacts in %s (work folder %s)\n',cfg.outputRoot,pwd);
fprintf('Stages: %s\n',strjoin(stages,', '));

report = struct('stage',{},'ok',{},'seconds',{},'message',{},'skipped',{});
for k = 1:numel(stages)
    name = stages{k};
    fprintf('\n===== Stage %d/%d: %s =====\n',k,numel(stages),name);
    t0 = tic; ok = true; msg = '';
    why = skipReason(name,cfg);
    if ~isempty(why)
        fprintf('Skipped: %s\n',why);
        report(end+1) = struct('stage',name,'ok',true,'seconds',0,'message',why,'skipped',true); %#ok<AGROW>
        continue
    end
    try
        switch name
            case 'synthetic', stageSynthetic(cfg);
            case 'benchmark', stageBenchmark(cfg);
            case 'city', stageCity(cfg);
            case 'mle', stageMLE(cfg);
            case 'integer', stageInteger(cfg);
            case 'costsel', stageCostSelection(cfg);
            case 'missing', stageMissingRoutes(cfg);
            case 'estimators', stageEstimators(cfg);
            case 'parma',     stageParma(cfg);
            case 'macros',    stageMacros(cfg);
            case 'figures',   stageFigures(cfg);
            case 'captions',  stageCaptions(cfg);
            case 'check',     stageCheck(cfg,report);
        end
    catch ME
        ok = false; msg = ME.message;
        try, msg = getReport(ME,'extended','hyperlinks','off'); catch, end
        fprintf(2,'Stage %s failed:\n%s\n',name,msg);
    end
    report(end+1) = struct('stage',name,'ok',ok,'seconds',toc(t0),'message',msg,'skipped',false); %#ok<AGROW>
    close all force
    if ~ok && ~cfg.continueOnError
        fprintf(2,'Stopping (set cfg.continueOnError=true to continue).\n');
        break
    end
end
fprintf('\n===== Summary =====\n');
for k = 1:numel(report)
    if report(k).skipped, tag = 'skipped'; elseif report(k).ok, tag = 'ok     '; else, tag = 'FAILED '; end
    fprintf('  %-10s %s  %8.1f s\n',report(k).stage,tag,report(k).seconds);
end
fprintf('Log: %s\n',logFile);
end

% =========================================================================
function stageSynthetic(cfg)
P = cfg.paths; S = cfg.synthetic;
opts = {'OutputDir',P.syntheticRun, ...
    'PaperResults',paperText(cfg),'PaperResultsFile',fullfile(P.captionData,'paper_results_synthetic.tex'), ...
    'FullLikelihoodReplicates',S.fullLikelihoodReplicates, ...
    'FullLikelihoodN',S.fullLikelihoodN,'Visible',cfg.visible, ...
    'UseParallel',logical(cfg.useParallel)};
if isfield(S,'numTruePaths') && ~isempty(S.numTruePaths), opts = [opts,{'NumTruePaths',S.numTruePaths}]; end
if S.useParmaGraph
    if isfile(P.parmaMat)
        R = loadResults(P.parmaMat);
        opts = [opts,{'Graph',R.G}];
        fprintf('Using the road graph saved by the Parma run.\n');
    else
        warning('buildPaperArtifacts:NoParmaGraph', ...
            'useParmaGraph is set but %s does not exist; downloading OSM roads instead.',P.parmaMat);
    end
end
if ~any(strcmp(opts(1:2:end),'Graph')) && isfield(S,'graphFile') && ~isempty(S.graphFile) && isfile(S.graphFile)
    g = load(S.graphFile,'G'); opts = [opts,{'Graph',g.G}];
    fprintf('Using the road graph bundled with the library: %s\n',S.graphFile);
end
runSyntheticOSMTomographyTests(opts{:},S.options{:});
end

function stageParma(cfg)
P = cfg.paths; C = cfg.parma;
assert(~isempty(C.compactFile),'buildPaperArtifacts:NoData', ...
    'Set cfg.parma.compactFile to the Parma traffic archive.');
ensureFolder(P.diagnostics);
caches = {'CacheDir',P.parmaOSMCache,'OfficialDataCacheDir',P.parmaOfficialData, ...
    'UseParallel',logical(cfg.useParallel), ...
    'SectionAuditCSV',fullfile(P.diagnostics,'parma_section_audit.csv'), ...
    'SensorResidualCSV',fullfile(P.diagnostics,'parma_sensor_residuals.csv'), ...
    'SensorDiagnosticsCSV',fullfile(P.diagnostics,'parma_sensor_diagnostics.csv'), ...
    'SensorSelectionCSV',fullfile(P.diagnostics,'parma_sensor_selection.csv'), ...
    'ValidationOutlierCSV',fullfile(P.diagnostics,'parma_validation_outliers.csv'), ...
    'IndistinguishableCSV',fullfile(P.diagnostics,'parma_indistinguishable_sensors.csv')};
res = runParmaTrafficTomography(C.compactFile,caches{:},C.options{:}, ...
    'MaxPathClasses',C.classCap,'SaveFile',P.parmaMat,'Plot',false);
res2 = [];
if ~isempty(C.screeningCap)
    fprintf('\nSecond Parma run with MaxPathClasses=%d (screening sensitivity).\n',C.screeningCap);
    res2 = runParmaTrafficTomography(C.compactFile,caches{:},C.options{:}, ...
        'MaxPathClasses',C.screeningCap,'SaveFile',P.parmaMatDoubled,'Plot',false);
end
if paperText(cfg), writeParmaPaperResults(res,fullfile(P.captionData,'paper_results_parma.tex'),res2); end
end

function stageBenchmark(cfg)
% Identification benchmark: run, then figure and macros (also redone by the
% figures and macros stages from the saved results).
P = cfg.paths; B = cfg.benchmark;
runIdentificationBenchmark('Replicates',B.replicates,'SampleSizes',B.sampleSizes, ...
    'RateScales',B.rateScales,'UseParallel',cfg.useParallel,'OutputDir',P.benchmarkRun,B.options{:});
benchmarkOutputs(cfg);
end

function stageCity(cfg)
% City-like benchmark: no longer in the paper; run only on request.
P = cfg.paths; C = cfg.city; cfg.city.enabled = true;
runCityBenchmark('Replicates',C.replicates,'UseParallel',cfg.useParallel,'OutputDir',P.cityBenchmarkRun,C.options{:});
benchmarkOutputs(cfg);
end

function stageEstimators(cfg)
% Parma with the PCL, the GMM and the robust GMM on identical strata.
P = cfg.paths;
C = cfg.parma;
% Same caches and class cap as the parma stage; its PCL results are reused when
% their settings match, so only the GMM is refitted.
compareParmaEstimators(C.compactFile,'PaperText',paperText(cfg),'CacheDir',P.parmaOSMCache,'OfficialDataCacheDir',P.parmaOfficialData, ...
    C.options{:},'MaxPathClasses',C.classCap,'UseParallel',cfg.useParallel, ...
    'PCLResultFile',P.parmaMat,'OutputDir',P.parmaEstimatorsRun, ...
    'MacroFile',fullfile(P.captionData,'paper_results_parma_estimators.tex'));
end

function stageMLE(cfg)
P = cfg.paths; M = cfg.mle;
runGMMvsMLE('Replicates',M.replicates,'UseParallel',cfg.useParallel,'OutputDir',P.mleComparisonRun,M.options{:});
benchmarkOutputs(cfg);
end

function stageMissingRoutes(cfg)
% Estimation with a route dictionary that misses active routes (Section 5.3).
P = cfg.paths; K = cfg.missing;
assert(isfile(P.syntheticMat),'buildPaperArtifacts:NoSynthetic','%s not found: run the synthetic stage first.',P.syntheticMat);
ensureFolder(P.missingRun);
runMissingRoutesStudy('SyntheticResults',P.syntheticMat,'Replicates',K.replicates,'Fractions',K.fractions, ...
    'UseParallel',logical(cfg.useParallel),'OutputDir',P.missingRun, ...
    'PaperResultsFile',textFile(cfg,fullfile(P.captionData,'paper_results_missing.tex')),K.options{:});
end

function stageCostSelection(cfg)
% Cost-aware projection selection on the corridor and the synthetic network (Section 5.4).
P = cfg.paths; K = cfg.costsel;
runCostSelectionStudy('Design','corridor','Replicates',K.replicates,'UseParallel',cfg.useParallel, ...
    'OutputDir',P.costSelectionRun,K.options{:});
if isfile(P.syntheticMat)
    runCostSelectionStudy('Design','synthetic','SyntheticResults',P.syntheticMat,'Replicates',K.replicates, ...
        'UseParallel',cfg.useParallel,'OutputDir',P.costSelectionRun,K.options{:});
else
    warning('buildPaperArtifacts:NoSynthetic','%s not found: run the synthetic stage first; cost selection on the synthetic network skipped.',P.syntheticMat);
end
benchmarkOutputs(cfg);
end

function stageInteger(cfg)
% Integer projections between composite and full likelihood (paper Section 5.2).
P = cfg.paths; I = cfg.integer;
runIntegerProjectionStudy('Replicates',I.replicates,'UseParallel',cfg.useParallel,'OutputDir',P.integerProjectionRun,I.options{:});
benchmarkOutputs(cfg);
end

function benchmarkOutputs(cfg)
% Figures and macros of both benchmarks from their saved results; each output is
% produced independently, so that one failure does not prevent the others.
P = cfg.paths;
mat = fullfile(P.benchmarkRun,'benchmark_results.mat');
if isfile(mat)
    d = load(mat,'R'); close all force
    try
        generateBenchmarkFigure(d.R,'RateScale',cfg.benchmark.figureScale,'OutputDir',P.benchmarkFigures,'Visible',cfg.visible);
    catch ME, warning('buildPaperArtifacts:BenchmarkFigure','Benchmark figure skipped: %s',ME.message); end
    if paperText(cfg)
        try
            writeBenchmarkResults(d.R,fullfile(P.captionData,'paper_results_benchmark.tex'),'RateScale',cfg.benchmark.figureScale);
        catch ME, warning('buildPaperArtifacts:BenchmarkMacros','Benchmark macros skipped: %s',ME.message); end
    end
end
ef = fullfile(P.parmaEstimatorsRun,{'parma_estimator_pcl.mat','parma_estimator_gmm.mat','parma_estimator_robustGMM.mat'});
if paperText(cfg) && all(cellfun(@isfile,ef(1:2)))
    efRob = []; if isfile(ef{3}), efRob = ef{3}; end      % robust run optional
    try
        writeParmaEstimatorResults(ef{1},ef{2},efRob,fullfile(P.captionData,'paper_results_parma_estimators.tex'));
    catch ME, warning('buildPaperArtifacts:EstimatorMacros','Parma estimator macros skipped: %s',ME.message); end
end
cf = {fullfile(P.costSelectionRun,'cost_selection_corridor.mat'), fullfile(P.costSelectionRun,'cost_selection_synthetic.mat')};
if any(cellfun(@isfile,cf))
    try
        generateCostSelectionFigure(cf,P.benchmarkFigures); close all force
    catch ME, warning('buildPaperArtifacts:CostSelFigure','Cost selection figure skipped: %s',ME.message); end
    if paperText(cfg)
        try
            writeCostSelectionResults(cf,fullfile(P.captionData,'paper_results_costsel.tex'));
        catch ME, warning('buildPaperArtifacts:CostSelMacros','Cost selection macros skipped: %s',ME.message); end
    end
end
imat = fullfile(P.integerProjectionRun,'integer_projection_results.mat');
if isfile(imat)
    di = load(imat,'R'); close all force
    try
        generateIntegerProjectionFigure(di.R,P.benchmarkFigures);
    catch ME, warning('buildPaperArtifacts:IntegerFigure','Integer projection figure skipped: %s',ME.message); end
    if paperText(cfg)
        try
            writeIntegerProjectionResults(di.R,fullfile(P.captionData,'paper_results_integer.tex'));
        catch ME, warning('buildPaperArtifacts:IntegerMacros','Integer projection macros skipped: %s',ME.message); end
    end
end
mmat = fullfile(P.mleComparisonRun,'mle_comparison_results.mat');
if paperText(cfg) && isfile(mmat)
    dm = load(mmat,'R');
    try
        writeMLEComparisonResults(dm.R,fullfile(P.captionData,'paper_results_mle.tex'));
    catch ME, warning('buildPaperArtifacts:MLEMacros','MLE comparison macros skipped: %s',ME.message); end
end
cmat = fullfile(P.cityBenchmarkRun,'city_benchmark_results.mat');
if isfile(cmat) && isfield(cfg,'city') && isfield(cfg.city,'enabled') && cfg.city.enabled   % not in the paper
    dc = load(cmat,'R'); close all force
    try
        generateCityBenchmarkFigure(dc.R,'OutputDir',P.benchmarkFigures,'Visible',cfg.visible);
    catch ME, warning('buildPaperArtifacts:CityFigure','City benchmark figure skipped: %s',ME.message); end
    if paperText(cfg)
        try
            writeCityBenchmarkResults(dc.R,fullfile(P.captionData,'paper_results_city.tex'));
        catch ME, warning('buildPaperArtifacts:CityMacros','City benchmark macros skipped: %s',ME.message); end
    end
end
end

function stageMacros(cfg)
P = cfg.paths;
benchmarkOutputs(cfg);
if isfile(P.syntheticMat)
    R = loadResults(P.syntheticMat);
    if isfield(R,'paperResults')
        S = R.paperResults;
    else
        S = collectSyntheticPaperResults(R.raw,R.settings,R.A,R.sensorSweepModels);
        if isfield(R,'sensorSpacingKm'), S.sensorSpacingKm = R.sensorSpacingKm; end
    end
    if (~isfield(S,'fullLik') || isempty(S.fullLik)) && cfg.synthetic.fullLikelihoodReplicates > 0
        S.fullLik = runningExampleFullLikelihood('N',cfg.synthetic.fullLikelihoodN, ...
            'R',cfg.synthetic.fullLikelihoodReplicates,'Seed',R.settings.Seed+999);
    end
    writeSyntheticPaperResults(S,fullfile(P.captionData,'paper_results_synthetic.tex'));
else
    warning('buildPaperArtifacts:NoSynthetic','%s not found; synthetic macros not written.',P.syntheticMat);
end
if isfile(P.parmaMat)
    R = loadResults(P.parmaMat);
    R2 = [];
    if isfile(P.parmaMatDoubled), R2 = loadResults(P.parmaMatDoubled); end
    writeParmaPaperResults(R,fullfile(P.captionData,'paper_results_parma.tex'),R2);
else
    warning('buildPaperArtifacts:NoParma','%s not found; Parma macros not written.',P.parmaMat);
end
end

function stageFigures(cfg)
P = cfg.paths;
benchmarkOutputs(cfg);
force = isfield(cfg.figures,'force') && cfg.figures.force;
synOut = fullfile(P.syntheticFigures,{'01_error_observations.pdf','02_error_outliers.pdf','03_error_impurity.pdf', ...
    '04_error_sensors.pdf','05_error_traveltime.pdf'});
if force || needsRedraw(synOut,{P.syntheticMat, which('generateSyntheticTomographyFigures')})
    generateSyntheticTomographyFigures(loadResults(P.syntheticMat), ...
        'OutputDir',P.syntheticFigures,'Visible',cfg.visible, ...
        'ShowNNLS',logical(cfg.figures.syntheticNNLS));
    close all force
else
    fprintf('Synthetic figures up to date: not redrawn.\n');
end
parOut = fullfile(P.parmaFigures,{'01_sensors_directions.pdf','02_peak_flows_by_day_type.pdf', ...
    '03_all_strata.pdf','04_validation_violins.pdf'});
if ~isfile(P.parmaMat)
    fprintf('Parma results not available (%s): Parma figures skipped.\n',P.parmaMat);
elseif force || needsRedraw(parOut,{P.parmaMat, which('generateParmaRequestedFigures')})
    generateParmaRequestedFigures(loadResults(P.parmaMat), ...
        'OutputDir',P.parmaFigures,'Visible',cfg.visible, ...
        'ShowNNLS',logical(cfg.figures.parmaNNLS),cfg.figures.parmaOptions{:});
    close all force
else
    fprintf('Parma figures up to date: not redrawn.\n');
end
if ~(isfield(cfg.figures,'diagnostics') && cfg.figures.diagnostics) || ~isfile(P.parmaMat)
    fprintf('Diagnostic figures skipped (cfg.figures.diagnostics = false).\n');
    return
end
generateComputationTimeFigure(P.syntheticMat,P.parmaMat, ...
    'OutputDir',P.computeFigures,'Visible',cfg.visible);
close all force
try
    generateSensorResidualFigure(P.parmaMat,'OutputDir',P.diagnostics,'Visible',cfg.visible);
catch ME
    warning('buildPaperArtifacts:ResidualMap','Sensor-residual map skipped: %s',ME.message);
end
close all force
try
    generateNRMSEvsSizeFigure(P.parmaMat,'OutputDir',P.diagnostics,'Visible',cfg.visible);
catch ME
    warning('buildPaperArtifacts:NRMSEvsSize','NRMSE-vs-size figure skipped: %s',ME.message);
end
close all force
try
    generateSensorEdgeMaps(P.parmaMat,'OutputDir',fullfile(P.diagnostics,'edge_maps'),'Visible',cfg.visible);
catch ME
    warning('buildPaperArtifacts:EdgeMaps','Sensor edge maps skipped: %s',ME.message);
end
end

function stageCaptions(cfg)
P = cfg.paths;
generatePaperCaptionData(P.parmaMat,P.syntheticMat,'OutputDir',P.captionData, ...
    'ParmaNNLSMarkers',logical(cfg.figures.parmaNNLS), ...
    'SyntheticNNLSMarkers',logical(cfg.figures.syntheticNNLS));
end

function stageCheck(cfg,report)
status = checkPaperArtifacts(cfg);   % accepts absolute or relative cfg.paths
writeBuildManifest(cfg,status,report);
missing = ~status.Present & ~status.Skipped;
if any(missing)
    error('buildPaperArtifacts:Missing','%d file(s) read by the paper are missing.',nnz(missing));
end
end

% =========================================================================
function R = loadResults(file)
assert(isfile(file),'buildPaperArtifacts:NoFile','Results file not found: %s',file);
data = load(file,'results');
assert(isfield(data,'results'),'buildPaperArtifacts:NoResults','%s has no variable results.',file);
R = data.results;
end

function ensureFolder(folder)
if ~isfolder(folder), mkdir(folder); end
end

function P = absolutePaths(cfg)
% All configured paths, made absolute with respect to cfg.outputRoot.
P = cfg.paths; f = fieldnames(P);
for k = 1:numel(f)
    v = char(P.(f{k}));
    if ~isAbsolute(v), v = fullfile(cfg.outputRoot,v); end
    P.(f{k}) = char(safeOutputPath(v));                 % error if inside the library
end
end

function tf = isAbsolute(p)
tf = startsWith(p,'/') || startsWith(p,'\\') || ~isempty(regexp(p,'^[A-Za-z]:[\\/]','once'));
end

function d = paperConfigDefaults(cfg,libRoot)
% Defaults from the paperConfig.m of THIS library, for the same paper folder.
here = pwd; back = onCleanup(@()cd(here));
cd(libRoot);                       % resolve paperConfig from the library
if isfield(cfg,'paperRoot') && ~isempty(cfg.paperRoot)
    d = paperConfig(cfg.paperRoot);
else
    d = paperConfig();
end
end

function s = withDefaults(s,d)
% Recursively add the fields of d that are missing in s.
f = fieldnames(d);
for k = 1:numel(f)
    if ~isfield(s,f{k})
        s.(f{k}) = d.(f{k});
    elseif isstruct(d.(f{k})) && isscalar(d.(f{k})) && isstruct(s.(f{k})) && isscalar(s.(f{k}))
        s.(f{k}) = withDefaults(s.(f{k}),d.(f{k}));
    end
end
end

function restoreSession(oldDir,oldPath)
diary('off');
path(oldPath);
cd(oldDir);
end

function tf = needsRedraw(outFiles,srcFiles)
% True if any output is missing or older than any source (results or figure code).
tf = false;
src = srcFiles(cellfun(@(f) ~isempty(f) && isfile(f), srcFiles));
if isempty(src), tf = true; return; end
tSrc = max(cellfun(@(f) dir(f).datenum, src));
for k = 1:numel(outFiles)
    if ~isfile(outFiles{k}) || dir(outFiles{k}).datenum < tSrc, tf = true; return; end
end
end

function why = skipReason(name,cfg)
% Reason to skip a stage whose input data are not available ('' to run it).
% The Parma counts are not public: without them, the stages that fit them are
% skipped, and the other stages run as usual.
why = '';
if any(strcmp(name,{'captions','macros'})) && ~paperText(cfg)
    why = sprintf('stage ''%s'' writes the paper''s text, not part of this version (cfg.paperText = false).',name); return
end
if any(strcmp(name,{'parma','estimators'}))
    f = '';
    if isfield(cfg,'parma') && isfield(cfg.parma,'compactFile'), f = char(cfg.parma.compactFile); end
    if isempty(f)
        why = 'Parma counts not configured (cfg.parma.compactFile is empty); they are not public.';
    elseif ~isfile(f)
        why = sprintf('Parma counts not found: %s.',f);
    end
end
end

function tf = paperText(cfg)
% True if the paper's result macros and captions are written (cfg.paperText).
tf = ~isfield(cfg,'paperText') || logical(cfg.paperText);
end

function f = textFile(cfg,f)
% The macro file to write, or '' when the paper's text is not written.
if ~paperText(cfg), f = ''; end
end
