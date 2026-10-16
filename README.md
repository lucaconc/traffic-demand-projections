# Tomography paper library

MATLAB code that produces every figure and every number of the paper
*Traffic demand estimation from flow measurements: a projection-based approach*
(Ardizzoni, Consolini, Laurini, Locatelli, Naz).

A single call runs the synthetic study, the Parma case study, the figures and
the LaTeX data files, and writes them in the folder layout that
`paper_scl.tex` reads.

## Requirements

- MATLAB R2020b or newer
- Optimization Toolbox (`lsqnonneg`, `fmincon`)
- Statistics and Machine Learning Toolbox (`poissrnd`, `gamrnd`, `prctile`)
- Mapping Toolbox (`shaperead`, `projcrs`, `projinv`), optional: needed to match sensors with the official street layer of the Municipality, as described in the paper; without it the matching uses OpenStreetMap only
- Internet access the first time only, to download OpenStreetMap roads (cached afterwards)
- The Parma traffic archive (provided by the Municipality of Parma and itcity; not distributed)

`setupPaperLib` checks these requirements and reports what is missing.

## Quick start

With the library in `<paper>/matlab/tomography_paper_lib`, the paper folder is found automatically:

```matlab
addpath('/path/to/paper/matlab/tomography_paper_lib')
cfg = paperConfig();                                % or paperConfig('/path/to/paper')
cfg.parma.compactFile = '/path/to/parma_archive';
report = buildPaperArtifacts(cfg);                  % all stages, several hours
```

`buildPaperArtifacts` uses **only the files of this library**: during the build
the MATLAB path is reset to the MathWorks default plus the library folders,
and the current folder is `<paper>/matlab/work`, which contains no code
(MATLAB searches the current folder first). The previous path and current
folder are restored afterwards. Set `cfg.isolatePath = false` to keep your path.
The existing Parma caches (`parma_osm_cache`, `parma_official_data`) in
`<paper>/matlab` are reused.

Then compile `paper_scl.tex`. Every output is written below
`/path/to/paper/matlab/`, which is where the paper looks for it.

## Stages

`buildPaperArtifacts(cfg,stages)` runs the listed stages in order.

| Stage | What it does | Needs |
|---|---|---|
| `synthetic` | Monte Carlo study on the OSM network (Section 6); saves the trials and writes `paper_results_synthetic.tex` | internet or cached OSM roads |
| `parma` | Parma case study (Section 7); main run plus a second run with a larger class cap; writes `paper_results_parma.tex` | `cfg.parma.compactFile` |
| `macros` | rewrites both result-macro files from the saved MAT files | saved results |
| `figures` | Fig. 2 (synthetic) and Figs. 3–5 (Parma) from the saved MAT files | saved results |
| `captions` | caption data and caption commands from the saved MAT files | saved results |
| `check` | lists the files read by the paper; writes `build_manifest.json` | — |

The default is `{'synthetic','parma','figures','captions','check'}`.
Typical partial rebuilds:

```matlab
buildPaperArtifacts(cfg,{'figures','captions','check'});   % after changing a plot
buildPaperArtifacts(cfg,{'macros','check'});               % after changing a summary
buildPaperArtifacts(cfg,{'parma','figures','captions','check'});
```

## Configuration

`paperConfig` returns a struct; edit it before calling `buildPaperArtifacts`.

| Field | Default | Meaning |
|---|---|---|
| `outputRoot` | `<paper>/matlab` | root of all outputs |
| `stages` | all but `macros` | default stages |
| `continueOnError` | `false` | continue after a failed stage |
| `useParallel` | `false` | fit the Parma strata in parallel, and parallelize the estimator's inner loops in the synthetic study, with `parfor` (Parallel Computing Toolbox); identical results |
| `synthetic.useParmaGraph` | `false` | reuse the road graph of the Parma run instead of downloading |
| `synthetic.fullLikelihoodReplicates` | `200` | replicates of the composite vs exact likelihood comparison (`0` skips it) |
| `synthetic.options` | `{}` | extra options of `runSyntheticOSMTomographyTests` |
| `parma.compactFile` | `''` | Parma traffic archive (required for the `parma` stage) |
| `parma.classCap`, `parma.screeningCap` | `50`, `[]` | class cap of the Parma run; class cap of an optional second run for the screening sensitivity (`[]`: no second run; e.g. `100` enables it) |
| `parma.options` | `{}` | extra options of `runParmaTrafficTomography` |
| `figures.syntheticNNLS`, `figures.parmaNNLS` | `false`, `true` | NNLS baseline markers (captions follow automatically) |
| `figures.parmaOptions` | `{}` | extra options of `generateParmaRequestedFigures`, e.g. `{'FlowBands',8,'FlowBandScale','quantile'}` for the peak-flow maps |

## Outputs and their place in the paper

Paths are relative to `outputRoot`.

| File | Used for |
|---|---|
| `synthetic_paper_figures/01..05_*.pdf` | Fig. 2(a)–(e) |
| `synthetic_paper_figures/simulation_settings.tex` | settings of Section 6 |
| `parma_realdata_tomography/parma_requested_figures/01_sensors_directions.pdf` | Fig. 3 |
| `parma_realdata_tomography/parma_requested_figures/03_all_strata.pdf` | Fig. 4 |
| `parma_realdata_tomography/parma_requested_figures/02_peak_flows_by_day_type.pdf` | Fig. 5(a) |
| `parma_realdata_tomography/parma_requested_figures/04_validation_violins.pdf` | Fig. 5(b) |
| `computation_figures/computation_times.pdf` | Fig. 6 |
| `paper_caption_data/paper_caption_data.tex`, `paper_figure_captions.tex` | caption text and numbers |
| `paper_caption_data/paper_results_synthetic.tex` | `\ResSyn...` macros of Section 6 |
| `paper_caption_data/paper_results_parma.tex` | `\ResParma...` macros of Section 7 |
| `paper_caption_data/parma_multistart_report.csv` | per-stratum sensitivity to the initial point (not read by the paper) |

Diagnostics not read by the paper go to `diagnostics/`: the section audit
(`parma_section_audit.csv`), the per-sensor NNLS residuals (`parma_sensor_residuals.csv`),
the sensor checks for faulty detectors and wrong directions (`parma_sensor_diagnostics.csv`,
from `diagnoseParmaSensors`), and the residual map with flagged sensors ringed (`sensor_residuals.pdf`).

Intermediate results (`synthetic_osm_tests/synthetic_results.mat`,
`parma_tomography_results*.mat`), CSV tables behind every figure, OSM caches
and the build logs (`logs/`) are kept in the same folder.

## Library layout

```
setupPaperLib.m          path and requirement checks
paperConfig.m            default configuration
buildPaperArtifacts.m    single entry point (stages, logging, error handling)
runPaperLibTests.m       unit tests
estimator/   estimateMuOptimizedProjections     composite-likelihood estimator (Section 4)
synthetic/   runSyntheticOSMTomographyTests     Monte Carlo study (Section 6)
             syntheticPaperMetrics, collectSyntheticPaperResults, runningExampleFullLikelihood
parma/       runParmaTrafficTomography          case study (Section 7)
             runParmaScreeningSensitivity       standalone two-run helper
             findStableParmaInterval            automatic calibration period (>= 90 days)
             diagnoseParmaSensors               flags faulty sensors and wrong directions (diagnostic)
             compareParmaAggregation            reruns the study at 10/20/30-minute aggregation and compares
             selectConsistentSensors            sensors whose deletion most improves the NNLS fit of the others
             optimizeParmaSensorSet             excludes them and compares the results with all sensors
             compareParmaStrataCriteria         BIC vs cross-validated stratification on the same test days
             diagnoseParmaFanoDecomposition     optional Fano diagnostic (not used by the paper)
             animateParmaTrafficFlows           optional animation (not used by the paper)
figures/     generateSyntheticTomographyFigures, generateParmaRequestedFigures,
             generateComputationTimeFigure, generateSensorResidualFigure (diagnostic),
             generateNRMSEvsSizeFigure (diagnostic), generateSensorEdgeMaps (diagnostic)
paperdata/   generatePaperCaptionData, writeSyntheticPaperResults,
             writeParmaPaperResults, writePaperMacros,
             multistartSummaryText, writeParmaMultistartReport
util/        checkPaperArtifacts, writeBuildManifest, machineSummary, checkLibraryShadowing,
             benchmarkEstimator (timing of the subproblem solvers on your machine)
tests/       testEstimator, testPaperData
```

Each function can also be called on its own; see its help text.

## Reproducibility

- All random number generators are seeded; a rebuild with the same data, code
  and MATLAB release reproduces the numbers.
- `build_manifest.json` records the date, the MATLAB release and toolboxes,
  the configuration, the stages run, and the size and SHA-256 hash of every
  file read by the paper.
- Each build writes a log of the command window to `logs/`.

## Tests

```matlab
results = runPaperLibTests();
```

The tests check the estimator on the running example of the paper:
recovery of the intensities, the rank condition, the signature fallback,
the sandwich standard errors against their theoretical value, and the
rejection of duplicate columns. They also check the macro writers. They take
about a minute and need no data or internet.

## Notes

- Old copies of library functions elsewhere (e.g. in `<paper>/matlab`) are
  ignored by `buildPaperArtifacts`; `setupPaperLib` warns about them because
  they would be used when calling functions interactively from that folder.
- CSV tables written by the pipelines with relative names go to `<paper>/matlab/work`.

- The weekday-profile figure (`05_weekday_sensor_profile`) is still produced
  by `generateParmaRequestedFigures` but is not used by the paper.
- The numbers in `paper_caption_data.tex` and in the result macros come only
  from saved results; nothing is typed by hand.


## Identification benchmark (Section 6)

Run explicitly (long: 4 matrices x 3 rate scales x 4 sample sizes x 200 replicates x 8 methods).
Matrices: all nonzero binary columns of {0,1}^k for k = 3, 4, 5, plus the running example.

    cfg = paperConfig();
    cfg.useParallel = true;                          % strongly recommended
    buildPaperArtifacts(cfg,{'benchmark','captions','check'});

Quick trial: `cfg.benchmark.replicates = 20;`. Results: `<paper>/matlab/benchmark_results/`
(MAT, raw and summary CSV); figure `benchmark_figures/benchmark_error.pdf`; table and settings
`paper_caption_data/paper_results_benchmark.tex`.

## City-like benchmark (Section 6)

Random street-grid instances; one factor at a time (sensors, paths, observations,
contamination); errors on class rates and solution times of all methods.

    cfg = paperConfig(); cfg.useParallel = true;
    cfg.city.replicates = 20;                       % trial; 100 for the paper
    buildPaperArtifacts(cfg,{'city','captions','check'});

Results: `<paper>/matlab/city_benchmark_results/`; figure `benchmark_figures/city_benchmark.pdf`;
macros `paper_caption_data/paper_results_city.tex`.

## Comparison with the exact maximum-likelihood estimator (Section 6)

On the easy and intermediate benchmark matrices the exact likelihood is computed by
successive convolution (fitFullLikelihood); runGMMvsMLE compares every method with it.

    cfg = paperConfig(); cfg.useParallel = true;
    cfg.mle.replicates = 20;                        % trial; 100 for the paper
    buildPaperArtifacts(cfg,{'mle','captions','check'});

Results: `<paper>/matlab/mle_comparison_results/`; table macros `paper_caption_data/paper_results_mle.tex`.

## Parma with the weighted projected scores

`runParmaTrafficTomography(...,'StratumEstimator','gmm')` (or `'robustGMM'`) fits every
stratum with the GMM (or its robust version) after the projected-CL fit, which still
supplies the rank and conditioning checks. To compare the three estimators on identical
strata, test days and class sets:

    [T,P] = compareParmaEstimators('/path/to/parma_archive');

T summarizes reliable strata, test NRMSE, local errors, misfit, failures and time;
P compares the per-stratum test NRMSE of each GMM version with the PCL's.
To use one of them in the paper build: `cfg.parma.options = {'StratumEstimator','gmm'};`

## Output location: never inside the library

No function writes inside the library folder tree. Every output path (options such as
`OutputDir`, `SaveFile`, `*CSV`, `CacheDir`, the writers' file arguments, the drivers'
result files and the build's configured paths) passes through `util/safeOutputPath`:

* relative paths are resolved in the current folder, or, if the current folder is inside
  the library, in the sibling folder `<parent of the library>/tomography_paper_output`
  (with a warning);
* absolute paths inside the library raise the error `tomography:OutputInsideLibrary`.

`paperConfig` recognizes the paper folder by any `.tex` file in `<paper>/` when the library
sits in `<paper>/matlab/`, and never falls back to a folder inside the library.

## Parma with corrected sensor matching

Some detectors are matched to the wrong road edge (opposite directions on one edge,
intersection approaches on one junction segment, turn lanes on the through road).
They are corrected with a table (`EdgeOverrideCSV`: columns SectionName, Action =
reverse | exclude | edges, EdgeIDs). The library holds a reference copy in
`data/parma_edge_overrides.csv`; your own copy in `<paper>/matlab/` takes precedence.

    out = runParmaWithOverrides();          % or scripts/run_parma_overrides.m

checks the setup, runs Parma with the table (results in a new file), verifies the
corrections, and writes the reliability tally, the error decomposition, the edge
candidates (`listEdgeCandidates`) and the edge maps to `<paper>/matlab/parma_override_run/`.
To find edge IDs for new corrections: `listEdgeCandidates(resultsFile, names)` lists the
road edges near each sensor with street name and direction; the edge maps print the IDs.

## Optional MEX acceleration

About half of the GMM's run time is the compound Poisson recursion. A C version
(`estimator/mex/compoundLogPmfMex.c`) computes the same terms in the same order:

    buildPaperMex();        % compiles into <parent of the library>/tomography_paper_mex
    clear functions         % so that loaded functions see the MEX file

The binary is written outside the library; `setupPaperLib` adds that folder to the path
when it exists. Without it, the MATLAB code is used and results are the same. The build
checks the MEX against the MATLAB recursion and fails if they differ. A C compiler must be
configured (`mex -setup C`; on macOS, Xcode's command line tools).


## Running without the Parma data

The Parma traffic counts were provided by the Municipality of Parma and itcity
and are not public; they may be obtained from the data owners subject to their
conditions. Without them, everything else runs:

```matlab
addpath('/path/to/tomography_paper_lib'); setupPaperLib();
cfg = paperConfig();          % leave cfg.parma.compactFile empty
report = buildPaperArtifacts(cfg,{'synthetic','missing','benchmark','mle','integer','costsel', ...
    'figures','captions','check'});
```

- The stages `parma` and `estimators` are reported as `skipped`, not failed.
- `figures` and `captions` produce everything except the Parma figures; the
  Parma captions are placeholders.
- `check` marks the Parma files as `SKIP` and passes if everything else is present.
- The synthetic study uses the road graph bundled in
  `data/central_parma_osm_graph.mat` (OpenStreetMap data, ODbL) when present,
  so no internet access is needed; otherwise it downloads OpenStreetMap roads.

## Sharing the library

```matlab
exportSyntheticGraph(cfg);                         % bundle the road graph (once)
makeShareableLibrary('tomography_paper_lib.zip');  % zip without restricted data
```

`makeShareableLibrary` leaves out `data/parma_edge_overrides.csv` (corrections to
the matching of the Municipality's counting sections; add
`'IncludeEdgeOverrides',true` only if the data agreement allows), compiled MEX
binaries, and every `.mat` file except the bundled graph, and warns about
absolute personal paths in the code.

## Licence

The code is released under the MIT License (see `LICENSE`). The bundled road
graph is OpenStreetMap data under the Open Database License (ODbL): road data
(c) OpenStreetMap contributors. The Parma traffic counts are not part of the
library and are not covered by this licence.

## Use of AI assistance

Parts of this code were developed with the assistance of AI tools, Claude
(Anthropic) and ChatGPT (OpenAI), which were used to draft, refactor and debug
functions and documentation. All code was reviewed, tested and validated by the
authors, who take full responsibility for it.

## About this version

This version reproduces the experiments of the paper and their figures. The code that wrote the paper's numerical sentences and captions is not included: cfg.paperText is false, and the stages 'captions' and 'macros' are reported as skipped. Numerical results are in the results files (.mat) of each study.
