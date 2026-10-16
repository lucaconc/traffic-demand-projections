function cfg = paperConfig(paperRoot,mode)
%PAPERCONFIG Default configuration for buildPaperArtifacts.
%
%   cfg = paperConfig()            library inside <paper>/matlab/: paper found automatically
%   cfg = paperConfig(paperRoot)
%   cfg = paperConfig([],'quick')   quick preview: few replicates, separate output folder
%
%   paperRoot is the folder that contains paper_scl.tex. When omitted, and
%   the library sits in <paper>/matlab/tomography_paper_lib, <paper> is used.
%   Every output is
%   written below cfg.outputRoot = fullfile(paperRoot,'matlab'), in the
%   layout read by the paper (see cfg.paths).  Edit the returned struct to
%   change any setting, then call buildPaperArtifacts(cfg).
%
%   Settings that MUST be provided:
%     cfg.parma.compactFile   the Parma traffic archive (data owners: the
%                             Municipality of Parma and itcity)
%
%   All other values reproduce the paper.  Seeds are fixed, so a rebuild
%   with the same data, code and MATLAB release gives the same numbers.

libRoot = fileparts(mfilename('fullpath'));
outName = 'matlab';                                     % outputs in <paper>/matlab
if nargin < 1 || isempty(paperRoot)
    [matlabDir,~] = fileparts(libRoot);
    [candidate,~] = fileparts(matlabDir);
    if ~isempty(dir(fullfile(candidate,'*.tex')))
        paperRoot = candidate;                          % library in <paper>/matlab/
    else
        % No paper next to the library (e.g. a clone of the public repository):
        % results go to a fixed folder beside the library, never inside it.
        paperRoot = fileparts(libRoot);
        outName = 'tomography_results';
        fprintf('No paper found next to the library: results go to %s\n',fullfile(paperRoot,outName));
    end
end
paperRoot = char(paperRoot);

cfg = struct();
cfg.paperRoot  = paperRoot;
cfg.outputRoot = fullfile(paperRoot,outName);
cfg.stages     = {'synthetic','missing','parma','figures','captions','check'};
cfg.continueOnError = false;     % true: log a failed stage and go on
cfg.visible    = 'off';          % figure visibility while building
% Use ONLY the files of this library: during the build the MATLAB path is
% reset to the MathWorks default plus the library folders, and the current
% folder is cfg.paths.work, which contains no code. Both are restored after.
cfg.isolatePath = true;
% Fit the Parma strata in parallel, and parallelize the synthetic estimator
% loops, with parfor (Parallel Computing Toolbox). Results are identical.
cfg.useParallel = false;

% ---- File layout read by paper_scl.tex (relative to cfg.outputRoot) ----
cfg.paths = struct( ...
    'syntheticRun',     'synthetic_osm_tests', ...
    'syntheticMat',     fullfile('synthetic_osm_tests','synthetic_results.mat'), ...
    'syntheticFigures', 'synthetic_paper_figures', ...
    'parmaMat',         'parma_tomography_results.mat', ...
    'parmaMatDoubled',  'parma_tomography_results_cap200.mat', ...
    'parmaFigures',     fullfile('parma_realdata_tomography','parma_requested_figures'), ...
    'computeFigures',   'computation_figures', ...
    'benchmarkRun',     'benchmark_results', ...         % identification benchmark (MAT, CSV)
    'benchmarkFigures', 'benchmark_figures', ...
    'cityBenchmarkRun', 'city_benchmark_results', ...
    'mleComparisonRun', 'mle_comparison_results', ...
    'integerProjectionRun','integer_projection_results', ...
    'missingRun',       'missing_routes', ...   % estimation with missing routes (stage 'missing')
    'costSelectionRun', 'cost_selection_results', ...   % cost-aware projection selection (stage 'costsel')   % integer projections vs full likelihood (stage 'integer')
    'parmaEstimatorsRun','parma_estimator_comparison', ...  % PCL vs GMM vs robust GMM (stage 'estimators')
    'captionData',      'paper_caption_data', ...
    'logs',             'logs', ...
    'diagnostics',      'diagnostics', ...               % audits, residual tables and maps
    'work',             'work', ...                  % current folder during the build (CSV tables)
    'parmaOSMCache',    'parma_osm_cache', ...       % OpenStreetMap cache of the Parma run
    'parmaOfficialData','parma_official_data');      % official street layer cache

% ---- Identification benchmark (Section 6): run only when the stage is requested ----
cfg.benchmark = struct();
cfg.benchmark.replicates = 100;               % paired Monte Carlo replicates
cfg.benchmark.sampleSizes = [100 500 2000 10000];
cfg.benchmark.rateScales = [2 10];            % overall rate scales (relative rates fixed)
cfg.benchmark.figureScale = 10;               % rate scale shown in the figure and table
cfg.benchmark.options = {};                   % extra options of runIdentificationBenchmark
cfg.city = struct('replicates',100,'options',{{}},'enabled',false);   % city-like benchmark: not in the paper (stage 'city' on request)
cfg.integer = struct('replicates',50,'options',{{}});
cfg.missing = struct('replicates',20,'fractions',[0 .1 .2 .3],'options',{{}});   % missing-routes study (stage 'missing')
cfg.costsel = struct('replicates',20,'options',{{}});  % cost-aware selection study (stage 'costsel')  % integer projections (stage 'integer')
cfg.mle = struct('replicates',100,'options',{{}});    % comparison with the exact MLE (stage 'mle')

% ---- Synthetic study (Section 6) ----
cfg.synthetic = struct();
cfg.synthetic.useParmaGraph = false;   % true: reuse the road graph saved by the Parma run
cfg.synthetic.graphFile = fullfile(fileparts(mfilename('fullpath')),'data','central_parma_osm_graph.mat');   % bundled graph (exportSyntheticGraph); used if present
cfg.synthetic.fullLikelihoodReplicates = 0;     % exact-likelihood comparison on the running example: not in the paper (200 to run it)
cfg.synthetic.fullLikelihoodN = 500;
cfg.synthetic.numTruePaths = 36;     % active paths among the 60 candidates (12: sparse demand, which favours PAIR)
cfg.synthetic.options = {'OutlierCounts',0};   % extra options of runSyntheticOSMTomographyTests; the spike
                                       % experiment is not in the paper (K = 0 only: the clean baseline, computed anyway)

% ---- Parma case study (Section 7) ----
cfg.parma = struct();
cfg.parma.compactFile = '';
cfg.parma.classCap = Inf;              % MaxPathClasses of the main run (Inf: all classes with positive NNLS pilot)
cfg.parma.screeningCap = [];           % class cap of an optional second run (screening sensitivity); [] = no second run
cfg.parma.options = {};                % extra name/value options of runParmaTrafficTomography

% ---- Figures and captions ----
cfg.paperText = false;   % write the paper's result macros and captions; false: results and figures only (public repository)

cfg.figures = struct();
cfg.figures.syntheticNNLS = false;     % NNLS markers in the synthetic figure
cfg.figures.parmaNNLS = true;          % NNLS markers in the validation figure
cfg.figures.parmaOptions = {};         % extra options of generateParmaRequestedFigures
cfg.figures.diagnostics = false;       % diagnostic figures (computation times, residual and edge maps): not in the paper
cfg.figures.force = false;             % true: redraw paper figures even if they are newer than their results and code
% ---- Quick preview ----
% Checks that every stage runs and the paper builds, in a fraction of the time.
% Outputs go to <paper>/matlab_quick, so they never replace the paper's results;
% the numbers are noisy and must not be quoted.
if nargin >= 2 && strcmpi(mode,'quick')
    cfg.outputRoot = [cfg.outputRoot '_quick'];
    cfg.benchmark.replicates = 10;  cfg.benchmark.sampleSizes = [100 1000];
    cfg.mle.replicates = 10;
    cfg.integer.replicates = 5;
    cfg.costsel.replicates = 3;
    cfg.missing.replicates = 3;
    cfg.synthetic.options = [cfg.synthetic.options, {'Replicates',5}];
    cfg.continueOnError = true;
end

end
