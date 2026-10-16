function results = runParmaTrafficTomography(compactFile,varargin)
%RUNPARMATRAFFICTOMOGRAPHY Apply optimized-projection Poisson tomography
%to the real Parma traffic-section archive.
%
%   results = runParmaTrafficTomography(compactFile)
%   results = runParmaTrafficTomography(compactFile,'Name',value,...)
%
% INPUT
%   compactFile  Compact archive produced by downloadTrafficParmaCompact.m
%                or packLegacyTrafficParma.m.
%
% PIPELINE
%   1. Read section names and GPS coordinates from the compact archive.
%   2. Download/cache a drivable OpenStreetMap road graph around the sensors.
%   3. Match each traffic section to a physical road segment using both GPS
%      distance and fuzzy road-name agreement. Direction words in the section
%      name (ingresso/uscita/centro/nord/sud/est/ovest) are used when present.
%   4. Use the native traffic time series to identify Poisson strata as
%      (time span) x (Workday/Saturday/Sunday-Holiday) x (calendar date span),
%      with configurable year/weekday pooling. Robust outlier removal and Fano-aware refinement
%      are applied before fitting the final stationary strata.
%   5. Generate many plausible shortest paths and collapse paths with the same
%      measured-section signature. Each retained latent variable therefore
%      represents an identifiable path-signature class, not an arbitrary OD pair.
%   6. Use NNLS only to screen the candidate signature classes, then estimate
%      their Poisson intensities with estimateMuOptimizedProjections.m.
%
% The paper model is
%       X_j ~ Poisson(mu_j), independent,
%       Y = A X.
% Here the columns of A are distinct nonzero measured-section signatures.
%
% IMPORTANT MODELING NOTE
%   The Parma archive does not contain OD labels or the exact paths followed by
%   individual vehicles. Consequently, without an external OD-zone/path model,
%   individual OD flows are not identifiable. This implementation estimates
%   aggregate intensities of distinguishable path-signature classes. That is the
%   strongest interpretation directly supported by the observed sections.
%
% SLICE-SPECIFIC CLASS REDUCTION
% SliceClassHandling='merge' (default): use the reduction described below.
% SliceClassHandling='strict': original behavior; skip strata containing any
% unobserved or duplicate class signature after sensor exclusion.
% Example: runParmaTrafficTomography(compactFile,'SliceClassHandling','strict')
% After sensor exclusion, identical NONZERO signatures are merged for fitting.
% Each reduced-class intensity is divided equally among its original classes.
% Zero-signature classes are omitted from fitting and assigned zero by convention,
% NOT estimated to have zero traffic. Reduction metadata flags them explicitly.
% Predictions at excluded/validation sensors therefore depend on this allocation.
% Original class indexing and member paths remain unchanged in saved outputs.
%
% MAIN OPTIONS
%   TrafficType            'FT' (total traffic count by default)
%   DateRange              [] -> automatically select a stable interval (below), or [start end]
%   AutoStableDateRange    true -> when DateRange=[], choose the best stable interval automatically
%   StableIntervalMinDays  90 minimum contiguous calibration interval
%   StableIntervalMinCoverage 0.90 overall sensor-day coverage target
%   StableIntervalMinSensorCoverage 0.80 lower-quartile per-sensor coverage target
%   StableIntervalMaxRelativeDeviation 0.30 max robust weekday-adjusted relative spread
%   StableIntervalMaxTrendFraction 0.20 max robust end-to-end weekday-adjusted drift
%   StableIntervalMinPositiveFraction 0.70 minimum fraction of finite daily counts > 0
%   StableIntervalRequireThresholds true -> stop if no >=90-day candidate passes quality thresholds
%   StableIntervalCSV      'parma_stable_interval_candidates.csv'
%   TimeWindowHours        [0 24] -> daily span used for Poisson segmentation
%   WeekdaysOnly           false (legacy option; all day types are used by default)
%   MergeYears             false -> fit each year separately; true -> pool all years by month/day
%   MergeWeekdays          true -> {Workday, Saturday, Sunday/Holiday}; false -> individual weekdays
%   MergeCorsie            true -> merge colocated sections differing only by a 'corsia' lane qualifier
%   LaneMergeMaxDistanceM  100 maximum GPS separation within one lane group
%   LaneMergeRequireAllLanes true -> merged count is missing unless every lane is measured at that interval
%   LaneMergeCSV           'parma_lane_groups.csv' audit table for lane aggregation
%   MinWindowCoverage      0.75
%   SignalReadChunkMB      128 -> batch signalData HDF5 reads in bounded chunks
%   MinPoissonSegmentMinutes 30 (minimum constant-intensity regime length)
%   MaxPoissonSegmentsPerDayType [44 11 11] time windows for workdays, Saturdays and
%                          Sundays/holidays (a scalar applies to every day type),
%                          roughly proportional to the number of days of each type
%   MaxPoissonSegmentsPerWeekday 11 legacy scalar name, used when the above is empty
%   MinSeasonSegmentDays   14 minimum calendar span length (two of each weekday)
%   MaxSeasonSegmentsPerTimeInterval 6 BIC date spans per time/day-type interval
%                          (each at least MinSeasonSegmentDays=14 long); a vector gives one
%                          value per day type, e.g. [6 3 3] for workdays/Saturdays/Sundays
%   (with these defaults there are at most (44+11+11) windows x 6 periods
%   = 396 strata;
%   with summed counts the BIC almost always selects the maximum number of
%   windows, so MaxPoissonSegmentsPerDayType effectively sets the count)
%   SamplesPerFittedClass  Inf -> every stratum fits all dictionary classes crossing its
%                          sensors (no per-stratum selection). On Parma, with a cap of one
%                          class per 2 intervals the test error was 0.085, without it 0.079,
%                          at equal reliability, multistart agreement and run time; a cap of
%                          one per 5 gave 0.140. With a finite value k, at most floor(nTrain/k)
%                          classes are fitted:
%                          the classes with the largest stratum NNLS pilot intensities,
%                          plus, for coverage, the crossing class with the largest pilot of
%                          every sensor left uncovered; the other classes are set to zero.
%                          Removing classes keeps the projections identifying.
%   MinSamplesPerClass     1 training interval per fitted class, else not fitted
%   NumPathSamples         300000 accepted paths (SignaturePatience, PathMaxTries scaled too)
%   PathCostNoise          0.5 -> candidate paths are shortest paths under edge costs
%                          multiplied by exp(0.5*N(0,1)), redrawn every PathNoiseBatch=100
%                          paths: a simple stochastic route-choice model that adds the
%                          alternative routes drivers actually use (0: exact shortest paths)
%   PathsPerSource         10 -> one shortest-path tree per sampled origin serves this many
%                          random destinations (one search instead of ten); 1 reproduces
%                          the previous sampling, with one search per path
%   SensorAnchoredShare    0.3 -> share of paths that start or end at a sensor: the path runs
%                          through one of the sensor's edges and continues to (or comes from)
%                          a random node. These represent trips that enter or leave the
%                          network right before or after a sensor (side streets, parking).
%   PathTerminals          'all' -> origins and destinations drawn among all nodes;
%                          'boundary+hubs' -> only boundary nodes and nodes of degree >= 4
%   ScreeningMethod        'nnls' -> the dictionary is chosen from the data: the classes with
%                          positive NNLS rate on the pooled sensor means (see MaxPathClasses);
%                          'frequency' -> the same number of classes, chosen without the data:
%                          those the route-choice model samples most often (for comparison)
%   MaxPathClasses         Inf -> the NNLS screening keeps every class with positive pilot
%                          intensity; NNLS rarely has more positive entries than estimation
%                          sensors, so the dictionary stays small (56 classes on Parma). A
%                          finite value keeps only the largest pilots: with 50 on Parma the
%                          test error was 0.079, without a cap 0.075, at equal reliability,
%                          multistart agreement and run time
%   EnsureSensorCoverage   true -> in addition, for every estimation sensor with positive
%                          mean that no kept class crosses, keep the crossing class with the
%                          largest pilot intensity (otherwise that sensor's predicted mean is
%                          zero, whatever the fit)
%   UseParallel            false -> true fits the strata in parallel with parfor (Parallel
%                          Computing Toolbox), and also parallelizes the cross-validated
%                          stratification. Every random draw uses a fixed seed per stratum,
%                          so results are identical to a serial run.
%   IncludeSingleSensors   true -> every stratum design includes the single-sensor projection
%                          of each sensor, with the same weight as the other projections,
%                          so that each sensor's own law enters the fit (without them the
%                          fitted means drift and many strata fail the misfit check)
%   EdgeOverrideCSV        '' -> optional table correcting the sensor-to-edge matching, one
%                          row per sensor, columns SectionName, Action and EdgeIDs:
%                            reverse  use the opposite direction: the reverse edge, or else
%                                     the nearest edge within 60 m with opposite bearing
%                                     (other carriageway of a dual road)
%                            exclude  do not use the sensor
%                            edges    use the edge IDs in EdgeIDs (space separated; IDs
%                                     are printed on the edge maps of the diagnostics).
%                                     Edge IDs can change when the network is rebuilt:
%                                     prefer street, which does not depend on them
%                            street   use the nearest edge (within 200 m of the sensor)
%                                     at least 10 m long (junction stubs are skipped)
%                                     whose street name contains Street and whose
%                                     direction of travel is within 50 deg of Direction
%                                     (N, NE, E, SE, S, SW, W, NW); columns Street and
%                                     Direction
%                          SectionName is matched to the cleaned section name (case and
%                          surrounding spaces ignored). Applied after the automatic
%                          matching and the name-based exclusions; each change is logged.
%                          Rows with an empty Action are skipped (other columns, such as
%                          notes, are ignored).
%   LocalSensorClasses     false -> if true, the class dictionary also contains, for every
%                          sensor, a class crossing only that sensor: traffic entering and
%                          leaving next to it on roads outside the network. Sensors that
%                          every sampled path crosses together (indistinguishable groups)
%                          can then have different means. Local classes go through the
%                          usual NNLS screening and class caps.
%   MeanAnchorWeight       0 -> no mean anchoring (default); w > 0 adds the single-sensor
%                          projections with total weight w (see the estimator help)
%   NumStarts              5 initial points per stratum fit; the per-stratum multi-start
%                          report (results.multistart*) assesses nonconvexity
%   Contamination          1e-10 weight of the background law in each scalar likelihood:
%                          counts impossible under Y=AX get weight zero instead of making
%                          the likelihood infinite. Keep it tiny: a larger weight (e.g. 0.01)
%                          also downweights the long tails of overdispersed counts and biases
%                          the intensities, and their fitted sensor means, downward.
%   RequireRankCondition   true -> fits failing the jump-class rank condition are unreliable
%   MaxSensitivityCondition 1e8 -> fits whose criterion Hessian is worse conditioned are unreliable
%   MaxExcessMisfit        0.10 -> fits whose training-mean misfit exceeds that of the NNLS
%                          fit of the same classes by more than 0.10 are unreliable. The
%                          misfit is ||A*mu - mean(Ytrain)||/||mean(Ytrain)||; NNLS gives
%                          the smallest misfit the retained classes allow, so a large NNLS
%                          misfit means the candidate paths cannot explain the sensor means
%                          (a modelling limit), while a large excess means a poor fit.
%   MaxTrainMeanRelError   Inf -> optional absolute bound on the misfit
%   TreatHolidaysAsSunday  true
%   IncludeParmaPatronHoliday true (13 January)
%   AdditionalHolidayDates [] exact extra holiday dates
%   MinSamplesPerStratum   30 minimum needed for train/test estimation (strata are still retained)
%   MinStratumIntervals    100 -> strata are built so that each holds at least this many
%                          intervals (days x time slots): a time window must reach it over all
%                          days of its type, and a calendar period over its days. This avoids
%                          strata too small to estimate; missing data can still reduce them.
%   ExcludeMissingSensorsPerSlice true -> never pass missing values to the estimator; retain sensors meeting SliceSensorMinCoverage, then drop incomplete rows jointly
%   SliceSensorMinCoverage 0.95 -> per-slice completeness threshold; set 1 for old strict behavior
%   AbnormalDayFactor      3 -> a sensor-day whose total count is below 1/3 or above 3 times
%                          the sensor's median daily total for the same day type (workday,
%                          Saturday, Sunday) is set missing: partial failures and gross
%                          miscounts that are not zero (Inf disables the rule)
%   ValidationOutlierCSV   'parma_validation_outliers.csv' -> the largest test errors, each
%                          classified as a data change on the test days or a fit error
%   ExcludeSectionNames    {} -> sections with exactly these names (as printed in the log
%                          and in parma_sensor_selection.csv) are not used
%   SensorSelection        true -> report the sensors whose deletion most improves the NNLS
%                          fit of the others (selectConsistentSensors); they are only listed,
%                          use ExcludeSectionNames or optimizeParmaSensorSet to exclude them
%   StratumEstimator       'pcl' -> each stratum is fitted with the projected composite
%                          likelihood; 'gmm' -> the efficiently weighted projected scores
%                          (estimateMuEfficientGMM), started at the PCL estimate; 'robustGMM'
%                          -> their robust version with contamination GMMEpsilon (0.01). The
%                          PCL fit always runs first: it provides the rank check and the
%                          conditioning diagnostic; all predictions, the misfit check and the
%                          validation use the final estimate. A failed GMM keeps the PCL
%                          estimate and is recorded in info.gmmError.
%   GMMProjections         'pclDesign' -> the GMM uses the projections chosen for the PCL
%                          fit of the stratum (D-optimal, full rank) plus every single
%                          sensor, with mean and variance conditions and no second start:
%                          about n+m projections instead of up to 300, much faster;
%                          'covered' -> all covered sensor subsets up to the identifying
%                          order (the estimator's default; slow for Parma-size strata)
%   GMMWeightIterations    2 -> weight matrix estimated twice; 1 halves the GMM time
%   MaxFlowPerHour         5000 -> a raw detector whose median nonzero count exceeds this many
%                          vehicles per hour is physically implausible (e.g. a cumulative
%                          counter) and is treated as missing
%   MaxSectionFlowPerHour  8000 -> same check after lanes are summed into one measurement
%                          (catches lane groups that sum too many detectors)
%   ExcludeSectionSuffixes ["(M)","(CL)"] and ExcludeSectionContains ["CICLI MILAZZO",
%                          "VEICOLO LUNGO","VEICOLI LUNGHI"] -> sections whose cleaned name
%                          ends with / contains these texts (case-insensitive) are not used.
%                          Plain text matching, applied together with ExcludeSectionNamePattern.
%   MaxSensorZeroShare     0.9 -> a raw detector (lane or section) whose counts are zero
%                          at least 90% of the time is treated as dead (all missing) before
%                          lanes are summed; a multi-lane measurement with a dead lane is
%                          then missing too, because its total cannot be known
%   ExcludeSectionNamePattern -> sections whose name matches are not used. Default:
%                          long-vehicle counts ('veicolo lungo'), sections whose name ends in
%                          '(M)' or '(CL)' (they duplicate other detectors on the same roads),
%                          and the bicycle counters 'CICLI MILAZZO'. Set '' to keep all.
%   AggregationMinutes     [] -> native interval (10 min). A multiple of it (e.g. 20 or 30)
%                          sums consecutive native counts after cleaning; an aggregated
%                          interval is missing if any of its native intervals is. Coarser
%                          intervals reduce inconsistencies from vehicles crossing successive
%                          sensors in different intervals, but leave fewer observations.
%   CleanFaultySensorData  true -> before stratification, mark as missing the data of a
%                          sensor that is off or stuck, so that each stratum keeps only
%                          sensors with reliable measurements (coverage rule above):
%                          * off day: daily total <= OffDayFraction (0.05) of the sensor's
%                            median daily total (if that median is >= OffMinDailyCount=20);
%                          * stuck day: >= StuckRunMinutes (90) of identical consecutive
%                            counts >= StuckMinValue (5);
%                          * zero run: >= ZeroRunMinutes (120) of zeros in time slots where
%                            the sensor's median count is >= ZeroRunMinTypical (5);
%                          and, within each stratum, a sensor with mean >= 5 whose Fano factor
%                          is below MinSensorFanoPerStratum (0.1) is dropped (nearly constant).
%   OutlierRemoval         true -> robustly remove anomalous native sensor counts
%   OutlierMADThreshold    3.5 robust MAD standard deviations on log1p(count)
%                          (the usual robust cut-off; 6 keeps only gross outliers)
%   OutlierMinSamples      20 minimum values per sensor/preliminary stratum
%   OutlierMaxFractionPerSensor 0.3 safety cap on removed values per sensor/stratum
%   FanoRefinement         false -> (true) recursively split overdispersed strata;
%                          the Fano factors are computed and reported in any case
%   MaxMedianFanoPerStratum 1.5 target median sensor Fano factor
%   FanoRefinementMaxPasses 1
%   FanoRefinementMinSamples 100 complete native intervals per child (parent >= 200)
%   FanoRefinementMinImprovement 0.05 fractional dispersion improvement required
%   FanoDecomposition      true -> year/weekday variance-decomposition diagnostic
%   FanoDecompositionCSV   'parma_fano_decomposition.csv'
%   FanoDecompositionSensorCSV 'parma_fano_decomposition_sensor.csv'
%   FanoDecompositionMinGroupSamples 5
%   UseAllSensors          true -> retain every accepted sensor with adequate coverage
%   MaxObservedSections    12 (used only when UseAllSensors=false)
%   MinObservedSections    4
%   MinCompleteDays        80
%   MinSectionCoverage     0.30 broad-window preselection only; stricter slice checks happen later
%   MaxSensorDistanceM     150
%   MinSensorMatchScore    0.45
%   OfficialDataAssistedMatching true -> use Comune official street geometry/names and traffic direction evidence
%   OfficialRoadAutoDownload true -> cache/download Comune Parma 2025 street-centerline SHP when possible
%   OfficialRoadShapefile  '' -> auto-discover/download; may point to an already extracted .shp
%   OfficialRoadSearchRadiusM 120
%   RejectAmbiguousMatches false; ambiguous assignments are retained but flagged for audit
%   MinSensorMatchMargin 0.02
%   MapMarginKm            2
%   RoadClasses            major/collector roads only by default; local
%                          residential, living_street, and unclassified roads
%                          are intentionally excluded from the OSM download
%   NumPathSamples         300000 accepted nonzero-signature shortest paths
%   PathMaxTries           3000000
%   SignaturePatience      300000 attempts without a new signature
%   MaxPathClasses         Inf (all classes with positive NNLS pilot)
%   TrainFraction          0.9 share of the days used for fitting (the others are test days)
%   TestDaySampling        'systematic' -> within each day type (working day, Saturday,
%                          Sunday/holiday), every k-th day in date order from a random start,
%                          k = round(1/(1-TrainFraction)) for working days and
%                          k = round(1/WeekendTestFraction) for Saturdays and Sundays/holidays:
%                          test days are spread evenly over the year, so short strata are
%                          rarely left without test days;
%   WeekendTestFraction    0.2 -> test share of weekend day types under systematic sampling
%                          (weekend strata have few days: with 10%, a third of them had none);
%                          'random' -> a simple random sample of all days. Either way the
%                          test days are drawn before any learning step.
%   SectionHoldoutValidation false -> (true) also reserve physical sections for validation
%   SectionHoldoutValidation false -> all sections are used for estimation; validation is
%                          temporal: in every stratum a share 1-TrainFraction (10%) of the
%                          data is held out at every sensor (true: also hold out sensors)
%   StrataCriterion        'bic' -> time windows and calendar periods maximize a piecewise-
%                          constant Poisson likelihood of the summed counts, with the number
%                          of segments chosen by BIC. 'cv' -> they minimize the squared error,
%                          on validation days, of the sensor means predicted by an NNLS fit on
%                          the training days (the number of segments is chosen by that error).
%                          Validation days for 'cv' are a share CVValidationShare (0.2) of the
%                          days not held out for the final test, so the test stays unbiased.
%                          CVTolerance (0.01): the fewest segments within 1% of the best error.
%                          For speed, boundaries lie on a grid of CVTimeBlockMinutes (30) in the
%                          day and CVCalendarBlockDays (7) in the calendar.
%   TemporalTestBlocks     'days' -> the held-out data are randomly chosen whole days of the
%                          stratum (seeded), so test intervals span the whole calendar period
%                          and are not adjacent to training intervals of the same day;
%                          'last' -> the last 10% of the stratum's intervals
%   SectionValidationFraction 0.125 target fraction of selected measurements held out
%   MinValidationCoupling  0.8 -> a sensor can be held out only if at least 80% of the
%                          sampled paths through it also cross an estimation sensor, so
%                          that its traffic is informed by the estimation data. A sensor
%                          whose traffic bypasses all estimation sensors cannot be
%                          predicted and would only test the path model, not the estimator.
%   SectionValidationMinCount 1 minimum held-out measurement when feasible
%   SectionValidationRows  [] optional explicit measured-row indices to hold out
%   SectionSplitSeed       [] -> Seed+1701
%   SectionSplitCSV        'parma_section_holdout_split.csv'
%   SectionValidationCSV   'parma_section_holdout_validation.csv'
%   Seed                   1234
%   Plot                   true
%   AnimateTraffic         false -> estimate/save only; animate later from saved MAT
%   AnimationDayType      'Workday' (or 'all')
%   AnimationSeasonDate    'auto' or e.g. '02-14'
%   AnimationYear          'auto' (latest separated year / pooled AllYears) or numeric year
%   AnimationPoissonSlice  'auto' or numeric Stratum ID; GUI chooses one slice at a time
%   AnimationWeekday       legacy alias; leave empty by default
%   AnimationDurationSeconds 30
%   AnimationFrameRate     30
%   AnimationMaxCars       250
%   AnimationVehicleSpeedKmH 120 (visual speed; faster for readable motion)
%   AnimationMinEdgeTraversalsPerSecond 0.35 minimum visible edge-fraction speed
%   AnimationVideoFile     '' (set an .mp4 filename to save video)
%   CacheDir               'parma_osm_cache'
%   UseCache               true
%   SaveFile               'parma_tomography_results.mat'
%   MappingCSV             'parma_section_edge_matches.csv'
%   ClassCSV               'parma_path_signature_classes.csv' (long form)
%   StrataCSV              'parma_poisson_strata.csv'
%   FlowErrorCSV           'parma_measured_predicted_flow_errors.csv'
%   SensorErrorCSV         'parma_sensor_error_summary.csv'
%   EstimatorFastMode      true -> common projection design + faster exact derivatives
%   ReuseProjectionDesign  true -> design C once from pooled cleaned data and reuse it
%   EstimatorFastMaxCandidates 400 candidate budget in fast mode
%   EstimatorOptions       {} overrides estimator defaults
%
% Requires estimateMuOptimizedProjections.m on the MATLAB path. Optimization
% Toolbox is required by that estimator. Mapping Toolbox is NOT required for the
% OSM download/matching; geographic plotting is optional.
%
% Map data: OpenStreetMap contributors, ODbL.

if nargin < 1 || isempty(compactFile)
    compactFile = fullfile('traffic_data_Parma','traffic_Parma_compact.mat');
end
compactFile = char(string(compactFile));
codeVersion = '2026-09-30-selectable-slice-class-handling-v41';
fprintf('runParmaTrafficTomography version: %s\n',codeVersion);
if exist(compactFile,'file') ~= 2
    error('runParmaTrafficTomography:MissingData', ...
        'Compact traffic archive not found: %s',compactFile);
end
if exist('estimateMuOptimizedProjections','file') ~= 2
    error('runParmaTrafficTomography:MissingEstimator', ...
        'estimateMuOptimizedProjections.m must be on the MATLAB path.');
end

%% Options
ip = inputParser;
ip.FunctionName = mfilename;
addParameter(ip,'TrafficType','FT',@(x)ischar(x)||isstring(x));
addParameter(ip,'DateRange',[],@(x)isempty(x)||(isdatetime(x)&&numel(x)==2));
addParameter(ip,'AutoStableDateRange',true,@isScalarLogicalLike);
addParameter(ip,'StableIntervalMinDays',90,@isPositiveInteger);
addParameter(ip,'StableIntervalCandidateLengthStepDays',7,@isPositiveInteger);
addParameter(ip,'StableIntervalMinCoverage',0.90,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=0&&x<=1);
addParameter(ip,'StableIntervalMinSensorCoverage',0.80,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=0&&x<=1);
addParameter(ip,'StableIntervalMaxRelativeDeviation',0.30,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=0);
addParameter(ip,'StableIntervalMaxTrendFraction',0.20,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=0);
addParameter(ip,'StableIntervalMinPositiveFraction',0.70,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=0&&x<=1);
addParameter(ip,'StableIntervalRequireThresholds',false,@isScalarLogicalLike);
addParameter(ip,'StableIntervalCSV','parma_stable_interval_candidates.csv',@(x)ischar(x)||isstring(x));
addParameter(ip,'TimeWindowHours',[0 24],@(x)isnumeric(x)&&numel(x)==2&&all(isfinite(x))&&x(1)>=0&&x(2)<=24&&x(2)>x(1));
addParameter(ip,'WeekdaysOnly',false,@isScalarLogicalLike);
addParameter(ip,'MergeYears',false,@isScalarLogicalLike);
addParameter(ip,'MergeWeekdays',true,@isScalarLogicalLike);
addParameter(ip,'MergeCorsie',true,@isScalarLogicalLike);
addParameter(ip,'LaneMergeMaxDistanceM',100,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0);
addParameter(ip,'LaneMergeRequireAllLanes',true,@isScalarLogicalLike);
addParameter(ip,'LaneMergeCSV','parma_lane_groups.csv',@(x)ischar(x)||isstring(x));
addParameter(ip,'MinWindowCoverage',0.75,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0&&x<=1);
addParameter(ip,'SignalReadChunkMB',128,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=1);
addParameter(ip,'MinPoissonSegmentMinutes',30,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0);
addParameter(ip,'MaxPoissonSegmentsPerWeekday',11,@isPositiveInteger);
addParameter(ip,'MaxPoissonSegmentsPerDayType',[44 11 11],@(x)isempty(x)||(isnumeric(x)&&isvector(x)&&all(isfinite(x))&&all(x>=1)&&all(x==floor(x))));
addParameter(ip,'MinSeasonSegmentDays',14,@isPositiveInteger);
addParameter(ip,'MaxSeasonSegmentsPerTimeInterval',6,@(x)isnumeric(x)&&isvector(x)&&all(isfinite(x))&&all(x>=1)&&all(x==floor(x)));
addParameter(ip,'TreatHolidaysAsSunday',true,@isScalarLogicalLike);
addParameter(ip,'IncludeParmaPatronHoliday',true,@isScalarLogicalLike);
addParameter(ip,'AdditionalHolidayDates',NaT(0,1),@(x)isempty(x)||isdatetime(x));
addParameter(ip,'MinSamplesPerStratum',30,@isPositiveInteger);
addParameter(ip,'MinStratumIntervals',100,@isPositiveInteger);
% Reliability of each stratum fit (strata failing a check are reported as
% not estimated, and their estimates are kept only for auditing).
addParameter(ip,'SamplesPerFittedClass',Inf,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(ip,'MinSamplesPerClass',1,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=0);
addParameter(ip,'RequireRankCondition',true,@isScalarLogicalLike);
addParameter(ip,'MaxSensitivityCondition',1e8,@(x)isnumeric(x)&&isscalar(x)&&x>1);
addParameter(ip,'MaxTrainMeanRelError',Inf,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(ip,'MaxExcessMisfit',0.10,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
addParameter(ip,'ExcludeMissingSensorsPerSlice',true,@isScalarLogicalLike);
addParameter(ip,'SliceClassHandling','merge', ...
    @(x)(ischar(x)||(isstring(x)&&isscalar(x)))&&any(strcmpi(x,{'merge','strict'})));
addParameter(ip,'SliceSensorMinCoverage',0.95,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0&&x<=1);
addParameter(ip,'OutlierRemoval',true,@isScalarLogicalLike);
addParameter(ip,'OutlierMADThreshold',3.5,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0);
addParameter(ip,'OutlierMinSamples',20,@isPositiveInteger);
addParameter(ip,'OutlierMaxFractionPerSensor',0.3,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=0&&x<0.5);
addParameter(ip,'FanoRefinement',false,@isScalarLogicalLike);
addParameter(ip,'MaxMedianFanoPerStratum',1.5,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=1);
addParameter(ip,'FanoRefinementMaxPasses',1,@isPositiveInteger);
addParameter(ip,'FanoRefinementMinSamples',100,@isPositiveInteger);
addParameter(ip,'FanoRefinementMinImprovement',0.05,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=0&&x<1);
addParameter(ip,'FanoMinMeanCount',1,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=0);
addParameter(ip,'FanoDecomposition',true,@isScalarLogicalLike);
addParameter(ip,'FanoDecompositionCSV','parma_fano_decomposition.csv',@(x)ischar(x)||isstring(x));
addParameter(ip,'FanoDecompositionSensorCSV','parma_fano_decomposition_sensor.csv',@(x)ischar(x)||isstring(x));
addParameter(ip,'FanoDecompositionMinGroupSamples',5,@isPositiveInteger);
addParameter(ip,'UseAllSensors',true,@isScalarLogicalLike);
addParameter(ip,'MaxObservedSections',12,@isPositiveInteger);
addParameter(ip,'MinObservedSections',4,@isPositiveInteger);
addParameter(ip,'MinCompleteDays',80,@isPositiveInteger);
addParameter(ip,'MinSectionCoverage',0.30,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=0&&x<=1);
addParameter(ip,'MaxSensorDistanceM',150,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0);
addParameter(ip,'DistanceScaleM',50,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0);
addParameter(ip,'MinSensorMatchScore',0.45,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=0&&x<=1);
addParameter(ip,'OfficialDataAssistedMatching',true,@isScalarLogicalLike);
addParameter(ip,'OfficialRoadAutoDownload',true,@isScalarLogicalLike);
addParameter(ip,'OfficialRoadShapefile','',@(x)ischar(x)||isstring(x));
addParameter(ip,'OfficialRoadNameField','auto',@(x)ischar(x)||isstring(x));
addParameter(ip,'OfficialDataCacheDir','parma_official_data',@(x)ischar(x)||isstring(x));
addParameter(ip,'OfficialRoadResourcePage','https://opendata.comune.parma.it/dataset/livello-cartografico-del-reticolo-stradale-di-parma-2025/resource/5da8981b-cb8d-4e28-960f',@(x)ischar(x)||isstring(x));
addParameter(ip,'OfficialRoadSearchRadiusM',120,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0);
addParameter(ip,'OfficialRoadDistanceScaleM',45,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0);
addParameter(ip,'RejectAmbiguousMatches',false,@isScalarLogicalLike);
addParameter(ip,'MinSensorMatchMargin',0.02,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=0&&x<=1);
addParameter(ip,'MapMarginKm',2,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0);
addParameter(ip,'RoadClasses', ...
    {'motorway','trunk','primary','secondary','tertiary', ...
     'motorway_link','trunk_link','primary_link','secondary_link','tertiary_link'}, ...
    @(x)iscellstr(x)||isstring(x)); %#ok<ISCLSTR>
addParameter(ip,'NumPathSamples',300000,@isPositiveInteger);
addParameter(ip,'PathCostNoise',0.5,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=0);
addParameter(ip,'PathsPerSource',10,@isPositiveInteger);
addParameter(ip,'SensorAnchoredShare',0.3,@(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=1);
addParameter(ip,'PathNoiseBatch',100,@isPositiveInteger);
addParameter(ip,'PathTerminals','all',@(x)any(strcmpi(x,{'all','boundary+hubs'})));
addParameter(ip,'PathMaxTries',3000000,@isPositiveInteger);
addParameter(ip,'SignaturePatience',300000,@isPositiveInteger);
addParameter(ip,'MaxPathClasses',Inf,@(x)isnumeric(x)&&isscalar(x)&&x>=1&&(isinf(x)||x==floor(x)));
addParameter(ip,'EnsureSensorCoverage',true,@isScalarLogicalLike);
addParameter(ip,'NumStarts',5,@isPositiveInteger);
addParameter(ip,'MeanAnchorWeight',0,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=0);
addParameter(ip,'IncludeSingleSensors',true,@(x)islogical(x)||isnumeric(x));
addParameter(ip,'LocalSensorClasses',false,@(x)islogical(x)||isnumeric(x));
addParameter(ip,'EdgeOverrideCSV','',@(x)ischar(x)||isstring(x));
addParameter(ip,'UseParallel',false,@isScalarLogicalLike);
addParameter(ip,'Contamination',1e-10,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=0&&x<1);
addParameter(ip,'MinPathClassPilot',1e-6,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=0);
addParameter(ip,'ScreeningMethod','nnls',@(x)any(strcmpi(x,{'nnls','frequency'})));
addParameter(ip,'TrainFraction',0.9,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0.5&&x<1);
addParameter(ip,'TestDaySampling','systematic',@(x)any(strcmpi(x,{'systematic','random'})));
addParameter(ip,'WeekendTestFraction',0.2,@(x)isnumeric(x)&&isscalar(x)&&x>0&&x<0.5);
addParameter(ip,'SectionHoldoutValidation',false,@isScalarLogicalLike);
addParameter(ip,'StrataCriterion','bic',@(x)any(strcmpi(x,{'bic','cv'})));
addParameter(ip,'CVTimeBlockMinutes',30,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(ip,'CVCalendarBlockDays',7,@isPositiveInteger);
addParameter(ip,'CVTolerance',0.01,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
addParameter(ip,'CVValidationShare',0.2,@(x)isnumeric(x)&&isscalar(x)&&x>0&&x<0.5);
addParameter(ip,'TemporalTestBlocks','days',@(x)any(strcmpi(x,{'days','last'})));
addParameter(ip,'SectionValidationFraction',0.125,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=0&&x<1);
addParameter(ip,'SectionValidationMinCount',1,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=0&&x==floor(x));
addParameter(ip,'SectionValidationRows',[],@(x)isempty(x)||(isnumeric(x)&&isvector(x)&&all(isfinite(x))&&all(x>=1)&&all(x==floor(x))));
addParameter(ip,'SectionSplitSeed',[],@(x)isempty(x)||(isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=0&&x==floor(x)));
addParameter(ip,'SectionSplitCSV','parma_section_holdout_split.csv',@(x)ischar(x)||isstring(x));
addParameter(ip,'SectionValidationCSV','parma_section_holdout_validation.csv',@(x)ischar(x)||isstring(x));
addParameter(ip,'Seed',1234,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=0&&x==floor(x));
addParameter(ip,'Plot',true,@isScalarLogicalLike);
addParameter(ip,'AnimateTraffic',false,@isScalarLogicalLike);
addParameter(ip,'AnimationDayType','Workday',@(x)ischar(x)||isstring(x));
addParameter(ip,'AnimationSeasonDate','auto',@(x)ischar(x)||isstring(x)||isdatetime(x));
addParameter(ip,'AnimationYear','auto',@(x)(ischar(x)||isstring(x))||(isnumeric(x)&&isscalar(x)&&isfinite(x)));
addParameter(ip,'AnimationPoissonSlice','auto',@(x)(ischar(x)||isstring(x))||(isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=1));
addParameter(ip,'AnimationWeekday','',@(x)ischar(x)||isstring(x));
addParameter(ip,'AnimationDurationSeconds',30,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0);
addParameter(ip,'AnimationFrameRate',30,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0);
addParameter(ip,'AnimationMaxCars',250,@isPositiveInteger);
addParameter(ip,'AnimationVehicleSpeedKmH',120,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0);
addParameter(ip,'AnimationMinEdgeTraversalsPerSecond',0.35,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=0);
addParameter(ip,'AnimationVideoFile','',@(x)ischar(x)||isstring(x));
addParameter(ip,'CacheDir','parma_osm_cache',@(x)ischar(x)||isstring(x));
addParameter(ip,'UseCache',true,@isScalarLogicalLike);
addParameter(ip,'OverpassURL','https://overpass-api.de/api/interpreter',@(x)ischar(x)||isstring(x));
addParameter(ip,'OverpassMirrors', ...
    {'https://overpass.private.coffee/api/interpreter', ...
     'https://maps.mail.ru/osm/tools/overpass/api/interpreter', ...
     'https://lz4.overpass-api.de/api/interpreter'}, ...
    @(x)iscellstr(x)||isstring(x)); %#ok<ISCLSTR>
addParameter(ip,'OverpassTimeout',180,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=30);
addParameter(ip,'OverpassTileKm',6,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0);
addParameter(ip,'OverpassRetries',3,@isPositiveInteger);
addParameter(ip,'OverpassRetryPause',4,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=0);
addParameter(ip,'UserAgent','MATLAB-Parma-NetworkTomography/1.0',@(x)ischar(x)||isstring(x));
addParameter(ip,'SaveFile','parma_tomography_results.mat',@(x)ischar(x)||isstring(x));
addParameter(ip,'MappingCSV','parma_section_edge_matches.csv',@(x)ischar(x)||isstring(x));
addParameter(ip,'MinValidationCoupling',0.8,@(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=1);
addParameter(ip,'SectionAuditCSV','parma_section_audit.csv',@(x)ischar(x)||isstring(x));
addParameter(ip,'SensorResidualCSV','parma_sensor_residuals.csv',@(x)ischar(x)||isstring(x));
addParameter(ip,'DiagnoseSensors',true,@isScalarLogicalLike);
addParameter(ip,'CleanFaultySensorData',true,@isScalarLogicalLike);
addParameter(ip,'AggregationMinutes',[],@(x)isempty(x)||(isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0));
addParameter(ip,'AbnormalDayFactor',3,@(x)isnumeric(x)&&isscalar(x)&&(x>1||isinf(x)));
addParameter(ip,'ValidationOutlierCSV','parma_validation_outliers.csv',@(x)ischar(x)||isstring(x));
addParameter(ip,'ExcludeSectionNames',{},@(x)iscell(x)||isstring(x));
addParameter(ip,'SensorSelection',true,@isScalarLogicalLike);
addParameter(ip,'IndistinguishableCSV','parma_indistinguishable_sensors.csv',@(x)ischar(x)||isstring(x));
addParameter(ip,'SensorSelectionCSV','parma_sensor_selection.csv',@(x)ischar(x)||isstring(x));
addParameter(ip,'StratumEstimator','pcl',@(x)any(strcmpi(x,{'pcl','gmm','robustGMM'})));
addParameter(ip,'GMMEpsilon',0.01,@(x)isnumeric(x)&&isscalar(x)&&x>0&&x<1);
addParameter(ip,'GMMWeightIterations',2,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
addParameter(ip,'GMMProjections','pclDesign',@(x)any(strcmpi(x,{'pclDesign','covered'})));
addParameter(ip,'MaxFlowPerHour',5000,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(ip,'MaxSectionFlowPerHour',8000,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(ip,'ExcludeSectionSuffixes',["(M)","(CL)"],@(x)isstring(x)||iscellstr(x)||ischar(x));
addParameter(ip,'ExcludeSectionContains',["CICLI MILAZZO","VEICOLO LUNGO","VEICOLI LUNGHI"],@(x)isstring(x)||iscellstr(x)||ischar(x));
addParameter(ip,'MaxSensorZeroShare',0.9,@(x)isnumeric(x)&&isscalar(x)&&x>0&&x<=1);
addParameter(ip,'ExcludeSectionNamePattern', ...
    '(?i)veicol\w*\s+lung|\(\s*m\s*\)\s*$|\(\s*cl\s*\)\s*$|(?<![A-Za-z0-9])cicli\s+milazzo(?![A-Za-z0-9])',@(x)ischar(x)||isstring(x));
addParameter(ip,'OffDayFraction',0.05,@(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<1);
addParameter(ip,'OffMinDailyCount',20,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
addParameter(ip,'StuckRunMinutes',90,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(ip,'StuckMinValue',5,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
addParameter(ip,'ZeroRunMinutes',120,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(ip,'ZeroRunMinTypical',5,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
addParameter(ip,'MinSensorFanoPerStratum',0.1,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
addParameter(ip,'DiagnoseSensorsTopK',20,@isPositiveInteger);
addParameter(ip,'SensorDiagnosticsCSV','parma_sensor_diagnostics.csv',@(x)ischar(x)||isstring(x));
addParameter(ip,'ClassCSV','parma_path_signature_classes.csv',@(x)ischar(x)||isstring(x));
addParameter(ip,'StrataCSV','parma_poisson_strata.csv',@(x)ischar(x)||isstring(x));
addParameter(ip,'FlowErrorCSV','parma_measured_predicted_flow_errors.csv',@(x)ischar(x)||isstring(x));
addParameter(ip,'SensorErrorCSV','parma_sensor_error_summary.csv',@(x)ischar(x)||isstring(x));
addParameter(ip,'EstimatorFastMode',true,@isScalarLogicalLike);
addParameter(ip,'ReuseProjectionDesign',true,@isScalarLogicalLike);
addParameter(ip,'EstimatorFastMaxCandidates',400,@isPositiveInteger);
addParameter(ip,'EstimatorOptions',{},@iscell);
parse(ip,varargin{:});
o = ip.Results;
% All output paths lie outside the library (see safeOutputPath).
o = guardOutputFields(o,{'StableIntervalCSV','LaneMergeCSV','FanoDecompositionCSV','FanoDecompositionSensorCSV','OfficialDataCacheDir','SectionSplitCSV','SectionValidationCSV','AnimationVideoFile','CacheDir','SaveFile','MappingCSV','SectionAuditCSV','SensorResidualCSV','ValidationOutlierCSV','IndistinguishableCSV','SensorSelectionCSV','SensorDiagnosticsCSV','ClassCSV','StrataCSV','FlowErrorCSV','SensorErrorCSV'});
o.TrafficType = string(o.TrafficType);
o.AutoStableDateRange = logical(o.AutoStableDateRange);
o.StableIntervalRequireThresholds = logical(o.StableIntervalRequireThresholds);
o.StableIntervalCSV = char(string(o.StableIntervalCSV));
o.WeekdaysOnly = logical(o.WeekdaysOnly);
o.MergeYears = logical(o.MergeYears);
o.MergeWeekdays = logical(o.MergeWeekdays);
o.MergeCorsie = logical(o.MergeCorsie);
o.LaneMergeRequireAllLanes = logical(o.LaneMergeRequireAllLanes);
o.LaneMergeCSV = char(string(o.LaneMergeCSV));
o.OfficialDataAssistedMatching = logical(o.OfficialDataAssistedMatching);
o.OfficialRoadAutoDownload = logical(o.OfficialRoadAutoDownload);
o.OfficialRoadShapefile = char(string(o.OfficialRoadShapefile));
o.OfficialRoadNameField = char(string(o.OfficialRoadNameField));
o.OfficialDataCacheDir = char(string(o.OfficialDataCacheDir));
o.OfficialRoadResourcePage = char(string(o.OfficialRoadResourcePage));
o.RejectAmbiguousMatches = logical(o.RejectAmbiguousMatches);
o.TreatHolidaysAsSunday = logical(o.TreatHolidaysAsSunday);
o.OutlierRemoval = logical(o.OutlierRemoval);
o.ExcludeMissingSensorsPerSlice = logical(o.ExcludeMissingSensorsPerSlice);
o.SliceClassHandling = lower(char(o.SliceClassHandling));
o.FanoRefinement = logical(o.FanoRefinement);
o.FanoDecomposition = logical(o.FanoDecomposition);
o.FanoDecompositionCSV = char(string(o.FanoDecompositionCSV));
o.FanoDecompositionSensorCSV = char(string(o.FanoDecompositionSensorCSV));
o.IncludeParmaPatronHoliday = logical(o.IncludeParmaPatronHoliday);
o.AdditionalHolidayDates = o.AdditionalHolidayDates(:);
if ~isempty(o.MaxPoissonSegmentsPerDayType)
    o.MaxPoissonSegmentsPerWeekday=o.MaxPoissonSegmentsPerDayType;
end
o.UseAllSensors = logical(o.UseAllSensors);
o.SectionHoldoutValidation = logical(o.SectionHoldoutValidation);
o.SectionValidationRows = unique(double(o.SectionValidationRows(:)),'stable');
if isempty(o.SectionSplitSeed), o.SectionSplitSeed=o.Seed+1701; end
o.SectionSplitCSV = char(string(o.SectionSplitCSV));
o.SectionValidationCSV = char(string(o.SectionValidationCSV));
o.EstimatorFastMode = logical(o.EstimatorFastMode);
o.ReuseProjectionDesign = logical(o.ReuseProjectionDesign);
o.Plot = logical(o.Plot);
o.AnimateTraffic = logical(o.AnimateTraffic);
o.AnimationDayType = char(string(o.AnimationDayType));
o.AnimationSeasonDate = o.AnimationSeasonDate;
o.AnimationYear = o.AnimationYear;
o.AnimationPoissonSlice = o.AnimationPoissonSlice;
o.AnimationWeekday = char(string(o.AnimationWeekday));
if ~isempty(o.AnimationWeekday), o.AnimationDayType=o.AnimationWeekday; end
if ~o.MergeWeekdays && strcmpi(o.AnimationDayType,'Workday'), o.AnimationDayType='Monday'; end
o.AnimationVideoFile = char(string(o.AnimationVideoFile));
o.UseCache = logical(o.UseCache);
o.CacheDir = char(string(o.CacheDir));
o.OverpassURL = char(string(o.OverpassURL));
o.OverpassMirrors = cellstr(string(o.OverpassMirrors));
o.UserAgent = char(string(o.UserAgent));
o.SaveFile = char(string(o.SaveFile));
o.MappingCSV = char(string(o.MappingCSV));
o.ClassCSV = char(string(o.ClassCSV));
o.StrataCSV = char(string(o.StrataCSV));
o.FlowErrorCSV = char(string(o.FlowErrorCSV));
o.SensorErrorCSV = char(string(o.SensorErrorCSV));
o.RoadClasses = cellstr(string(o.RoadClasses));
requestedDateRange=o.DateRange;
if o.AutoStableDateRange && isempty(o.DateRange) && exist('findStableParmaInterval','file')~=2
    error('runParmaTrafficTomography:MissingStableIntervalDiagnostic', ...
        ['AutoStableDateRange=true requires findStableParmaInterval.m on the MATLAB path. ', ...
         'Place it beside runParmaTrafficTomography.m or set AutoStableDateRange=false.']);
end

if ~o.UseAllSensors && o.MinObservedSections > o.MaxObservedSections
    error('runParmaTrafficTomography:Options', ...
        'MinObservedSections cannot exceed MaxObservedSections when UseAllSensors=false.');
end

oldRng = rng;
cleanupRng = onCleanup(@() rng(oldRng)); %#ok<NASGU>
rng(o.Seed,'twister');

fprintf('\n============================================================\n');
fprintf('Parma real-data optimized-projection network tomography\n');
fprintf('Traffic archive: %s\n',compactFile);
fprintf('Map data: OpenStreetMap contributors (ODbL)\n');
fprintf('OSM road classes: %s\n',strjoin(string(o.RoadClasses),', '));
fprintf('Use all usable sensors: %d\n',o.UseAllSensors);
fprintf('Official-data-assisted sensor matching: %d; reject ambiguous=%d; margin >= %.3f\n', ...
    o.OfficialDataAssistedMatching,o.RejectAmbiguousMatches,o.MinSensorMatchMargin);
fprintf('Path search: target=%d, max tries=%d, patience=%d, max classes=%d\n', ...
    o.NumPathSamples,o.PathMaxTries,o.SignaturePatience,o.MaxPathClasses);
fprintf('Estimator speed mode: %d; common projection design: %d; fast candidate budget=%d\n', ...
    o.EstimatorFastMode,o.ReuseProjectionDesign,o.EstimatorFastMaxCandidates);
fprintf('Section hold-out validation=%d; target fraction=%.1f%%; split seed=%d\n', ...
    o.SectionHoldoutValidation,100*o.SectionValidationFraction,o.SectionSplitSeed);
fprintf('Automatic stable date range=%d; minimum=%d days; coverage target=%.1f%%\n', ...
    o.AutoStableDateRange,o.StableIntervalMinDays,100*o.StableIntervalMinCoverage);
fprintf('signalData reader: chunked matfile I/O, target %.0f MB/read\n',o.SignalReadChunkMB);
fprintf('Poisson strata: time span x weekday grouping x calendar date span x year grouping\n');
fprintf('  Merge weekdays=%d; merge years=%d\n',o.MergeWeekdays,o.MergeYears);
fprintf('  Merge corsie=%d; lane radius=%.1f m; require all lanes=%d\n', ...
    o.MergeCorsie,o.LaneMergeMaxDistanceM,o.LaneMergeRequireAllLanes);
fprintf('  Slice class handling: %s\n',o.SliceClassHandling);
fprintf('  Missing-data protection=%d; per-slice sensor coverage >= %.1f%%, then complete-row filtering\n', ...
    o.ExcludeMissingSensorsPerSlice,100*o.SliceSensorMinCoverage);
fprintf('  Outlier removal=%d; robust log-MAD threshold=%.2f; min samples=%d; cap=%.1f%%/sensor-stratum\n', ...
    o.OutlierRemoval,o.OutlierMADThreshold,o.OutlierMinSamples,100*o.OutlierMaxFractionPerSensor);
fprintf('  Native-interval observations; Fano refinement=%d, target median Fano <= %.2f\n', ...
    o.FanoRefinement,o.MaxMedianFanoPerStratum);
fprintf('  Fano decomposition diagnostic=%d (year, weekday, year x weekday)\n',o.FanoDecomposition);
fprintf('  BIC: up to %s time intervals per day type, up to %s coarse date spans/time interval\n', ...
    mat2str(o.MaxPoissonSegmentsPerWeekday),mat2str(o.MaxSeasonSegmentsPerTimeInterval));
fprintf('  Calendar spans: minimum %d days (~%.1f months); Fano refinement obeys same minimum\n', ...
    o.MinSeasonSegmentDays,o.MinSeasonSegmentDays/30.44);
fprintf('Animate immediately after estimation: %d; day type=%s; max visual cars=%d\n', ...
    o.AnimateTraffic,o.AnimationDayType,o.AnimationMaxCars);
fprintf('============================================================\n');

%% Load compact metadata
meta = load(compactFile, ...
    'sectionId','sectionName','sectionLat','sectionLon', ...
    'typeNames','dateStart','dateEnd','nRecords');
sectionId = string(meta.sectionId(:));
sectionName = string(meta.sectionName(:));
sectionLat = double(meta.sectionLat(:));
sectionLon = double(meta.sectionLon(:));
nSections = numel(sectionId);

validGPS = isfinite(sectionLat) & isfinite(sectionLon) & ...
    sectionLat > 40 & sectionLat < 50 & sectionLon > 5 & sectionLon < 20;
% Reject isolated bad coordinates before constructing the download box.
if any(validGPS)
    latMed=median(sectionLat(validGPS));
    lonMed=median(sectionLon(validGPS));
    dMed=haversineKm(sectionLat,sectionLon,latMed,lonMed);
    validGPS = validGPS & dMed <= 50;
end
if nnz(validGPS) < o.MinObservedSections
    error('runParmaTrafficTomography:GPS', ...
        'Only %d sections have usable GPS coordinates.',nnz(validGPS));
end
fprintf('Traffic sections: %d total, %d with valid GPS.\n',nSections,nnz(validGPS));

%% Download/cache OSM graph using the actual sensor extent
bbox = sensorBoundingBox(sectionLat(validGPS),sectionLon(validGPS),o.MapMarginKm);
[osm,osmCacheFile] = getOSMRoadData(bbox,o);
G = buildDirectedRoadGraph(osm);
fprintf('OSM graph: %d nodes, %d directed edges.\n',numnodes(G),numedges(G));

%% Match sensor sections to physical OSM road segments
cityCenter = [median(sectionLat(validGPS)),median(sectionLon(validGPS))];
officialRoads=loadParmaOfficialRoadLayerLocal(o);
match = matchSectionsToRoads(G,sectionId,sectionName,sectionLat,sectionLon,validGPS,cityCenter,o,officialRoads);

% Merge multiple physical detector lanes ("corsia") into one measurement row.
% The representative keeps one directed OSM edge/orientation, while the native
% traffic data are summed across every member lane before tomography.
[match,laneMergeTable] = mergeCorsiaMatchesLocal(match,sectionId,sectionName,o);
% Names can carry invisible padding (null, control or non-breaking-space
% characters) after the text, which would defeat end-of-name patterns.
cleanNames=cleanSectionNamesLocal(match.SectionName);
if ~isempty(o.ExcludeSectionNames)
    byList=ismember(cleanNames,cleanSectionNamesLocal(o.ExcludeSectionNames(:)));
    if any(byList)
        match.Accepted(byList)=false;
        fprintf('Sections excluded by name list: %d\n',nnz(byList));
    end
    missingNames=setdiff(cleanSectionNamesLocal(o.ExcludeSectionNames(:)),cleanNames);
    if ~isempty(missingNames)
        warning('runParmaTrafficTomography:ExcludeNames','%d names in ExcludeSectionNames were not found.',numel(missingNames));
    end
end
% Plain-text exclusions (independent of regular-expression details)
upNames=upper(regexprep(cleanNames,'\s+',''));
byText=false(size(cleanNames));
for t=string(o.ExcludeSectionSuffixes(:)).'
    byText=byText | endsWith(upNames,upper(regexprep(t,'\s+','')));
end
for t=string(o.ExcludeSectionContains(:)).'
    byText=byText | contains(upNames,upper(regexprep(t,'\s+','')));
end
if any(byText)
    match.Accepted(byText)=false;
    fprintf('Sections excluded by name text: %d (%s)\n',nnz(byText),strjoin(cellstr(cleanNames(byText)),'; '));
end
if strlength(string(o.ExcludeSectionNamePattern))>0
    byName=~cellfun(@isempty,regexpi(cellstr(cleanNames),char(o.ExcludeSectionNamePattern),'once'));
    if any(byName)
        match.Accepted(byName)=false;
        fprintf('Sections excluded by name pattern: %d (%s)\n',nnz(byName), ...
            strjoin(cellstr(string(match.SectionName(byName))),'; '));
    end
end
if strlength(string(o.EdgeOverrideCSV))>0
    match=applyEdgeOverridesLocal(match,cleanNames,char(o.EdgeOverrideCSV),G);
end
if o.MergeCorsie && ismember('LaneCount',laneMergeTable.Properties.VariableNames)
    bigGroups=find(laneMergeTable.LaneCount>=4);
    if ~isempty(bigGroups)
        fprintf('Lane groups with 4 or more detectors (check that they are one carriageway): %d\n',numel(bigGroups));
        for gI=bigGroups(:).'
            nameCol=intersect({'SectionName','Label','GroupName'},laneMergeTable.Properties.VariableNames);
            if ~isempty(nameCol), nm=string(laneMergeTable.(nameCol{1})(gI)); else, nm=string(gI); end
            fprintf('  %-45s %d detectors\n',char(strtrim(nm)),laneMergeTable.LaneCount(gI));
        end
    end
end
if o.MergeCorsie
    nMergedGroups=nnz(laneMergeTable.LaneCount>1);
    nMergedRaw=sum(max(laneMergeTable.LaneCount-1,0));
    fprintf('Lane merging: %d multi-lane groups collapsed (%d redundant lane rows removed).\n', ...
        nMergedGroups,nMergedRaw);
else
    fprintf('Lane merging disabled.\n');
end
if ~isempty(o.LaneMergeCSV)
    writetable(laneMergeTable,o.LaneMergeCSV);
    fprintf('Lane grouping table written to: %s\n',o.LaneMergeCSV);
end

nGoodMatch = nnz(match.Accepted);
fprintf('Accepted sensor-road measurements after lane merging: %d/%d raw sections.\n',nGoodMatch,nSections);
if ismember('OfficialRoadUsed',match.Properties.VariableNames)
    fprintf('  Official centerline evidence used: %d; HIGH confidence: %d; ambiguous raw matches: %d.\n', ...
        nnz(match.OfficialRoadUsed),nnz(match.MatchConfidence=="HIGH"),nnz(match.Ambiguous));
end
if nGoodMatch < o.MinObservedSections
    error('runParmaTrafficTomography:Matching', ...
        ['Only %d section-to-road measurements passed the thresholds after lane merging. ', ...
         'Increase MaxSensorDistanceM or lower MinSensorMatchScore after inspecting the mapping table.'], ...
        nGoodMatch);
end

%% Build broad-window daily counts only for sensor quality/coverage selection
candidateSections = find(match.Accepted);
candidateGroups = match.LaneMemberSections(candidateSections);
[Ycand,obsDates,coverage,intervalMinutes,dataInfo] = aggregateDailyTraffic( ...
    compactFile,candidateGroups,o);

% Attach group-level coverage to every raw lane member for auditability.
match.DataCoverage(:) = NaN;
for ig=1:numel(candidateSections)
    members=match.LaneMemberSections{candidateSections(ig)};
    match.DataCoverage(members)=coverage(ig);
end

% Export a human-auditable mapping table with both geometric/name scores and
% actual data coverage plus lane-group membership.
matchExport = match;
matchExport.EdgeIDs = strings(height(match),1);
matchExport.LaneMemberSectionIndices = strings(height(match),1);
for i=1:height(match)
    if ~isempty(match.DirectedEdgeIDs{i})
        matchExport.EdgeIDs(i) = strjoin(string(match.DirectedEdgeIDs{i}),',');
    end
    if ~isempty(match.LaneMemberSections{i})
        matchExport.LaneMemberSectionIndices(i)=strjoin(string(match.LaneMemberSections{i}),',');
    end
end
% Explain why every raw section is or is not available for broad-window selection.
matchExport.BroadCoverageEligible=false(height(matchExport),1);
matchExport.PreselectionStatus=strings(height(matchExport),1);
for i=1:height(matchExport)
    if ismember('IsLaneRepresentative',matchExport.Properties.VariableNames) && ~matchExport.IsLaneRepresentative(i)
        matchExport.PreselectionStatus(i)="merged_into_lane_group";
    elseif ~isfinite(matchExport.Lat(i)) || ~isfinite(matchExport.Lon(i))
        matchExport.PreselectionStatus(i)="missing_gps";
    elseif isempty(match.DirectedEdgeIDs{i})
        matchExport.PreselectionStatus(i)="no_directed_osm_edge";
    elseif isfinite(matchExport.DistanceM(i)) && matchExport.DistanceM(i)>o.MaxSensorDistanceM
        matchExport.PreselectionStatus(i)="too_far_from_osm";
    elseif isfinite(matchExport.MatchScore(i)) && matchExport.MatchScore(i)<o.MinSensorMatchScore
        matchExport.PreselectionStatus(i)="low_match_score";
    elseif o.RejectAmbiguousMatches && matchExport.Ambiguous(i)
        matchExport.PreselectionStatus(i)="ambiguous_match";
    elseif ~matchExport.Accepted(i)
        matchExport.PreselectionStatus(i)="not_accepted";
    elseif ~isfinite(matchExport.DataCoverage(i))
        matchExport.PreselectionStatus(i)="no_traffic_data";
    elseif matchExport.DataCoverage(i)<o.MinSectionCoverage
        matchExport.PreselectionStatus(i)="low_broad_coverage";
    else
        matchExport.PreselectionStatus(i)="eligible_after_broad_coverage";
        matchExport.BroadCoverageEligible(i)=true;
    end
end
matchExport.DirectedEdgeIDs = [];
matchExport.LaneMemberSections = [];
matchExport.SegmentKey = string(matchExport.SegmentKey);
writetable(matchExport,o.MappingCSV);
fprintf('Section/edge mapping written to: %s\n',o.MappingCSV);
statusCats=unique(matchExport.PreselectionStatus);
fprintf('Sensor preselection audit:\n');
for is=1:numel(statusCats)
    fprintf('  %-30s %d\n',char(statusCats(is)),nnz(matchExport.PreselectionStatus==statusCats(is)));
end

fprintf('Native traffic aggregation interval: %.6g minutes.\n',intervalMinutes);
fprintf('Candidate days before sensor selection: %d.\n',numel(obsDates));

%% Select observed measurements
% Every candidate column is now one physical measurement. A multi-lane group
% already contains the SUM across its member corsie.
[selectedLocal,completeRows] = selectObservedSections( ...
    match(candidateSections,:),Ycand,o);
selectedSections = candidateSections(selectedLocal);   % representative raw row for metadata/map
selectedGroups = candidateGroups(selectedLocal);       % raw archive rows contributing to each measurement
Yselection = round(Ycand(:,selectedLocal));
datesSelection = obsDates;

if size(Yselection,2) < o.MinObservedSections
    error('runParmaTrafficTomography:TooFewObservedSections', ...
        'Only %d usable observed measurements remained after coverage/lane selection.',size(Yselection,2));
end

m = size(Yselection,2);
nBroadComplete = nnz(completeRows);
fprintf('\nSelected %d road measurements from %d accepted post-merge matches.\n', ...
    m,numel(candidateSections));
fprintf('Broad-window complete days across all selected measurements: %d/%d (diagnostic only).\n', ...
    nBroadComplete,numel(obsDates));
LaneCountSel=double(match.LaneCount(selectedSections));
disp(table((1:m).',match.SectionID(selectedSections),match.SectionName(selectedSections),LaneCountSel, ...
    match.RoadName(selectedSections),match.OfficialRoadName(selectedSections),match.DistanceM(selectedSections), ...
    match.MatchScore(selectedSections),match.MatchConfidence(selectedSections),match.DataCoverage(selectedSections), ...
    match.DirectionRule(selectedSections),match.AssignmentSource(selectedSections), ...
    'VariableNames',{'MeasuredRow','SectionID','SectionName','LaneCount','OSMRoad','OfficialRoad','DistanceM','MatchScore','Confidence','Coverage','Direction','AssignmentSource'}));

%% Generate topology signatures and choose the section hold-out BEFORE temporal stratification
% The split uses only OSM/path topology, never traffic values. This allows the
% subsequent BIC/Fano temporal partition itself to be learned from estimation
% sections only, avoiding indirect leakage from held-out section measurements.
measuredEdgeSets = match.DirectedEdgeIDs(selectedSections);
sig = generatePathSignatureClasses(G,measuredEdgeSets,o);
if isempty(sig.Aall)
    error('runParmaTrafficTomography:NoSignatures', ...
        'No shortest path crossed the selected observed road sections.');
end
fprintf('Generated %d distinct nonzero measured-edge signatures.\n',size(sig.Aall,2));
if logical(o.LocalSensorClasses)
    [sig,nLocal] = addLocalSensorClasses(sig,measuredEdgeSets,G);
    fprintf('Added %d local single-sensor classes (traffic entering and leaving next to a sensor).\n',nLocal);
end
[sectionEstimationMask,sectionValidationMask,sectionSplitInfo] = ...
    chooseSectionHoldoutLocal(sig.Aall,m,o,sig.hits);
estimationRows=find(sectionEstimationMask);
validationRows=find(sectionValidationMask);
mEst=numel(estimationRows); mVal=numel(validationRows);
fprintf('Section split: %d estimation / %d validation measurements (target validation=%d).\n', ...
    mEst,mVal,sectionSplitInfo.TargetValidationCount);
if mVal>0
    fprintf('  validation measured rows: %s\n',strjoin(string(validationRows.'),', '));
    if isfield(sectionSplitInfo,'Coupling')
    fprintf('  share of sampled paths through each validation sensor that also cross an estimation sensor:\n');
    for r=validationRows(:).'
        fprintf('    %4d %-40s %.2f\n',r,char(string(match.SectionName(selectedSections(r)))), ...
            sectionSplitInfo.Coupling(r));
    end
    end
end

%% Choose a stable >=90-day estimation interval from ESTIMATION sensors only
% This step occurs after the topology-only section split. Held-out validation
% flows therefore cannot influence the selected dates. The broad daily counts
% are used only for interval quality scoring; final tomography still uses native
% interval observations from the compact archive.
stableInterval=struct();
if o.AutoStableDateRange && isempty(requestedDateRange)
    stableInput=struct('dates',datesSelection,'Y',Yselection(:,sectionEstimationMask));
    stableInterval=findStableParmaInterval(stableInput, ...
        'MinDays',o.StableIntervalMinDays, ...
        'CandidateLengthStepDays',o.StableIntervalCandidateLengthStepDays, ...
        'MinCoverage',o.StableIntervalMinCoverage, ...
        'MinSensorCoverage',o.StableIntervalMinSensorCoverage, ...
        'MaxRelativeDeviation',o.StableIntervalMaxRelativeDeviation, ...
        'MaxTrendFraction',o.StableIntervalMaxTrendFraction, ...
        'MinPositiveFraction',o.StableIntervalMinPositiveFraction, ...
        'UseEstimationSensorsOnly',false, ...
        'SingleYear',true, ...
        'CSV',o.StableIntervalCSV, ...
        'Plot',false);
    stableInterval.autoSelected=true;
    if o.StableIntervalRequireThresholds && ~stableInterval.passesThresholds
        error('runParmaTrafficTomography:NoStableInterval', ...
            ['No contiguous interval of at least %d days met the configured data-quality/stability thresholds. ', ...
             'Inspect %s, relax StableInterval* thresholds, or set StableIntervalRequireThresholds=false.'], ...
            o.StableIntervalMinDays,o.StableIntervalCSV);
    end
    o.DateRange=[stableInterval.startDate stableInterval.endDate];
    fprintf('Using automatically selected stable interval: %s to %s (%d days).\n', ...
        datestr(stableInterval.startDate,'yyyy-mm-dd'),datestr(stableInterval.endDate,'yyyy-mm-dd'),stableInterval.days);
elseif ~isempty(requestedDateRange)
    stableInterval=struct('autoSelected',false,'startDate',requestedDateRange(1), ...
        'endDate',requestedDateRange(2),'days',days(requestedDateRange(2)-requestedDateRange(1))+1, ...
        'passesThresholds',NaN,'score',NaN);
end

% Display one point/edge per merged measurement.
if o.Plot
    plotParmaMapping(G,match,selectedSections,double(match.Lat),double(match.Lon),string(match.SectionName));
end

%% Divide observations by day type, time span, and calendar date span
o.CVDictionary = double(sig.Aall(sectionEstimationMask,:)); o.CVHits = double(sig.hits(:));
[strata,nativeInfo] = buildPoissonStrata(compactFile,selectedGroups,o,sectionEstimationMask);
if isempty(strata)
    error('runParmaTrafficTomography:NoPoissonStrata', ...
        ['No weekday/time/date/year stratum contains usable observations. ', ...
         'Check sensor coverage, DateRange, and TimeWindowHours.']);
end

fprintf('\nDetected %d usable weekday/time/date/year Poisson strata.\n',numel(strata));
for q=1:numel(strata)
    fprintf('  %2d  %-14s %-8s %s-%s | %s-%s | N=%d, retained estimation sensors=%d/%d\n',q,strata(q).DayType,strata(q).YearLabel, ...
        minuteLabel(strata(q).StartMinute),minuteLabel(strata(q).EndMinute), ...
        strata(q).SeasonStart,strata(q).SeasonEnd,size(strata(q).Y,1), ...
        nnz(strata(q).SensorMask),mEst);
end

% Pooled full-width data provide per-row means, but path-class screening below
% uses ONLY rows assigned to the estimation subset. Individual likelihoods also
% use only slice-complete estimation sensors. Held-out section flows are reserved
% for validation and never enter NNLS screening, projection design, or fitting.
YpoolFull = vertcat(strata.YFull);
% Class screening and the sensor diagnostics use non-test rows only.
poolLearn = true(size(YpoolFull,1),1);
if isfield(nativeInfo,'outerTestDates') && isfield(strata,'Dates')
    poolDates = vertcat(strata.Dates);
    if numel(poolDates)==size(YpoolFull,1), poolLearn = ~ismember(poolDates(:),nativeInfo.outerTestDates); end
end
meanMeasured = mean(YpoolFull(poolLearn,:),1,'omitnan').';
% Sensors with no usable data in any stratum (all removed by the cleaning and
% plausibility checks) are not used at all: filling in a mean from the coarse
% selection data would reintroduce exactly the faulty values the checks removed.
noStrataData = ~isfinite(meanMeasured);
if any(noStrataData & sectionEstimationMask)
    dropped = find(noStrataData & sectionEstimationMask);
    fprintf('Sensors without usable data in any stratum, not used: %d\n',numel(dropped));
    % Why: days removed by the cleaning, and coverage within the strata
    covAll = nan(numel(strata),m);
    for qq = 1:numel(strata)
        if isfield(strata,'SensorCoverage') && numel(strata(qq).SensorCoverage)==m
            covAll(qq,:) = strata(qq).SensorCoverage(:).';
        end
    end
    fc = struct('offDays',zeros(m,1),'stuckDays',zeros(m,1),'abnormalDays',zeros(m,1),'zeroRunIntervals',zeros(m,1));
    if isfield(nativeInfo,'faultCleaning'), fc = nativeInfo.faultCleaning; end
    if ~isfield(fc,'abnormalDays'), fc.abnormalDays = zeros(m,1); end
    for k = dropped(:).'
        cv = covAll(:,k); cv = cv(isfinite(cv));
        if isempty(cv), covText = 'no data in the strata';
        else, covText = sprintf('coverage in strata: median %.0f%%, max %.0f%% (needs %.0f%%)', ...
                100*median(cv),100*max(cv),100*o.SliceSensorMinCoverage); end
        fprintf('  %-45s cleaning removed %d off, %d stuck, %d abnormal days; %s\n', ...
            char(strtrim(string(match.SectionName(selectedSections(k))))), ...
            fc.offDays(k),fc.stuckDays(k),fc.abnormalDays(k),covText);
    end
    sectionEstimationMask(dropped) = false;
end
meanMeasured(noStrataData) = 0;     % never used: these rows are outside the estimation mask

% Dispersion is assessed only on sensors actually used within each stationary
% slice. Store full-width vectors with NaN for excluded sensors for auditing.
stratumFano = cell(numel(strata),1);
allFano = zeros(0,1);
for q=1:numel(strata)
    mask=logical(strata(q).SensorMask(:)) & sectionEstimationMask(:);
    Yf=strata(q).YFull(:,mask);
    mq=mean(Yf,1).';
    fUsed=var(Yf,0,1).' ./ max(mq,eps);
    fUsed(mq<o.FanoMinMeanCount)=NaN;
    fq=nan(m,1); fq(mask)=fUsed;
    stratumFano{q}=fq;
    allFano=[allFano;fUsed(isfinite(fUsed))]; %#ok<AGROW>
end
medianFano = median(allFano,'omitnan');
stratumMedianFanoPre=zeros(numel(strata),1);
for q=1:numel(strata)
    mask=logical(strata(q).SensorMask(:)) & sectionEstimationMask(:);
    stratumMedianFanoPre(q)=medianFanoLocal(strata(q).YFull(:,mask),o.FanoMinMeanCount);
end
nHighFano=nnz(isfinite(stratumMedianFanoPre) & stratumMedianFanoPre>o.MaxMedianFanoPerStratum);
fprintf('Median within-stratum measured-edge Fano factor: %.4g (Poisson benchmark = 1).\n',medianFano);
fprintf('Strata above target median Fano %.2f: %d/%d.\n', ...
    o.MaxMedianFanoPerStratum,nHighFano,numel(strata));
if isfinite(medianFano) && (medianFano > 1.5 || medianFano < 0.67)
    warning('runParmaTrafficTomography:PoissonDispersion', ...
        ['Counts still depart substantially from the Poisson variance benchmark after ', ...
         'day-type/time/season stratification and Fano refinement. Interpret model-based uncertainty cautiously.']);
end

%% Screen the path-signature dictionary using ESTIMATION sections only
% This prevents held-out section flows from leaking into NNLS class screening.
[Aest,classSelection,pilotAll,pilotFit] = selectPathClasses( ...
    sig.Aall(sectionEstimationMask,:),meanMeasured(sectionEstimationMask),o);
if strcmpi(o.ScreeningMethod,'frequency')
    % Same dictionary size as NNLS (before coverage), but the most frequently
    % sampled signatures: tests whether choosing the classes from the data matters.
    nNNLS = nnz(pilotAll > o.MinPathClassPilot);
    if isfinite(o.MaxPathClasses), nNNLS = min(nNNLS,o.MaxPathClasses); end
    [~,ordHits] = sort(double(sig.hits(:)),'descend');
    [Aest,classSelection] = selectPathClassesGiven(sig.Aall(sectionEstimationMask,:), ...
        meanMeasured(sectionEstimationMask),ordHits(1:min(nNNLS,numel(ordHits))),pilotAll,o);
    fprintf('Frequency screening: %d most sampled classes (NNLS would keep %d), %d after coverage.\n', ...
        min(nNNLS,numel(ordHits)),nNNLS,numel(classSelection));
end
selectedClass = classSelection(:);
% Per-sensor residual of the dictionary-wide NNLS fit of the pooled means
% (estimation sensors fitted; validation sensors predicted).
predGlobalNNLS = double(sig.Aall)*pilotAll;
relResGlobalNNLS = (predGlobalNNLS-meanMeasured)./max(meanMeasured,eps);
% ---- Indistinguishable sensors: identical rows of the path dictionary ----------
% Every sampled path through one sensor of such a group crosses the others too,
% so the model can only give them one common flow. Groups of opposite directions
% indicate a matching or network problem (same directed edge, or a road present
% in one direction only); groups along one road indicate missing turning routes.
indistinguishable=table();
try
    estI=find(sectionEstimationMask);
    [~,~,gid]=unique(double(sig.Aall(estI,:)>0),'rows');
    cnt=accumarray(gid,1);
    grpRows=find(cnt(gid)>1);
    if ~isempty(grpRows)
        k=estI(grpRows); g=gid(grpRows);
        nmI=strtrim(string(match.SectionName(selectedSections(k))));
        eds=cellfun(@(e)mat2str(sort(round(double(e(:)')))),measuredEdgeSets(k),'UniformOutput',false);
        [g,ordG]=sort(g); k=k(ordG); nmI=nmI(ordG); eds=eds(ordG);
        indistinguishable=table(g,k,nmI,meanMeasured(k),predGlobalNNLS(k),string(eds(:)), ...
            'VariableNames',{'Group','Measurement','SectionName','MeanCount','PredictedNNLS','MatchedEdges'});
        ug=unique(g);
        fprintf('\nIndistinguishable sensor groups (identical path signatures): %d groups, %d sensors\n',numel(ug),numel(k));
        for gg=ug(:).'
            r=find(g==gg);
            mc=meanMeasured(k(r)); ratio=max(mc)/max(min(mc),0.1);
            sameEdges=numel(unique(eds(r)))<numel(r);
            fprintf('  group %3d: counts %s (max/min %.1f)%s\n    %s\n',gg,mat2str(round(mc(:).',1)),ratio, ...
                ternaryLocal(sameEdges,'  <- some share the same matched edges',''),strjoin(cellstr(nmI(r)),' | '));
        end
        if ~isempty(o.IndistinguishableCSV), writetable(indistinguishable,char(o.IndistinguishableCSV)); end
    end
catch MEind
    warning('runParmaTrafficTomography:Indistinguishable','Indistinguishable-sensor report skipped: %s',MEind.message);
end

% Sensors whose deletion most improves the NNLS fit of the others (report only)
sensorSelection=struct();
if o.SensorSelection
    try
        estIdx=find(sectionEstimationMask);
        selNames=cellstr(strtrim(string(match.SectionName(selectedSections(estIdx)))));
        % Speed: restrict the dictionary to the classes used by the full NNLS fit
        % plus the 1000 most sampled ones (deleting a sensor rarely needs others).
        [~,ordH]=sort(double(sig.hits(:)),'descend');
        cols=false(size(sig.Aall,2),1); cols(ordH(1:min(1000,numel(ordH))))=true; cols=cols | pilotAll(:)>0;
        sensorSelection=selectConsistentSensors(sig.Aall(estIdx,cols),meanMeasured(estIdx),selNames);
        sensorSelection.recommendedMeasurements=estIdx(sensorSelection.recommended);
        fprintf('Sensor selection: NNLS NRMSE %.3f with all estimation sensors, %.3f without the %d recommended:\n', ...
            sensorSelection.baseNRMSE,sensorSelection.finalNRMSE,numel(sensorSelection.recommended));
        st=sensorSelection.steps;
        for k=1:height(st)
            fprintf('  %2d. %-45s relative residual %+6.2f -> NRMSE %.3f\n',k,st.Name{k},st.RelResidual(k),st.NRMSEAfter(k));
        end
        if ~isempty(o.SensorSelectionCSV), writetable(st,char(o.SensorSelectionCSV)); end
    catch MEsel
        warning('runParmaTrafficTomography:SensorSelection','Sensor selection skipped: %s',MEsel.message);
    end
end
% Faulty sensors and wrong directions (leave-one-out, reversal, pair swap, data quality)
sensorDiagnostics = table();
if o.DiagnoseSensors
    try
        fprintf('Sensor diagnostics on the %d sensors with the largest NNLS residuals...\n',o.DiagnoseSensorsTopK);
        sensorDiagnostics = diagnoseParmaSensors(G,measuredEdgeSets,sig,meanMeasured, ...
            sectionEstimationMask,match.SegmentKey(selectedSections),match.SectionName(selectedSections), ...
            YpoolFull,'TopK',o.DiagnoseSensorsTopK);
        % Sensors outside the estimation mask are not validation sensors here:
        % label them correctly and blank their (meaningless) residuals.
        notUsed = ~sectionEstimationMask(:) & ~sectionValidationMask(:);
        sensorDiagnostics.Role(notUsed) = "not used (no stratum data)";
        sensorDiagnostics.RelResidual(notUsed) = NaN;
        sensorDiagnostics.MeanCount(notUsed) = NaN;
        flagged = find(strlength(sensorDiagnostics.Flag)>0);
        fprintf('Flagged sensors: %d (dictionary NNLS misfit %.3f)\n',numel(flagged), ...
            sensorDiagnostics.Properties.UserData.baseMisfit);
        for k = flagged(:).'
            fprintf('  %4d %-40s %-10s %s: %s\n',k,char(sensorDiagnostics.SectionName(k)), ...
                char(sensorDiagnostics.Role(k)),char(sensorDiagnostics.Flag(k)),char(sensorDiagnostics.Advice(k)));
        end
        if ~isempty(o.SensorDiagnosticsCSV)
            writetable(sensorDiagnostics,char(o.SensorDiagnosticsCSV));
            fprintf('  sensor diagnostics written to %s\n',char(o.SensorDiagnosticsCSV));
        end
    catch MEdiag
        warning('runParmaTrafficTomography:SensorDiagnostics','Sensor diagnostics skipped: %s',MEdiag.message);
    end
end
A = double(sig.Aall(:,selectedClass));       % full rows retained only for prediction/validation
Aest = double(Aest);                         % estimation-row view used to define identifiable classes
p = size(A,2);
if p==0
    error('runParmaTrafficTomography:NoEstimationSignatures', ...
        'No identifiable path-signature classes remain after the section hold-out split.');
end
fprintf('Retained %d latent path-signature classes using estimation sections only.\n',p);
fprintf('Estimation-section NNLS screening relative residual: %.4g.\n',pilotFit);
if pilotFit > 0.25
    warning('runParmaTrafficTomography:PathDictionaryFit', ...
        ['The generated path-signature dictionary explains the pooled ESTIMATION-section mean only ', ...
         'approximately (relative residual %.3f). Increase NumPathSamples/PathMaxTries ', ...
         'or MaxPathClasses.'],pilotFit);
end

repEdgePaths = sig.representativeEdgePaths(selectedClass);
repNodePaths = sig.representativeNodePaths(selectedClass);
classEdgePaths = sig.memberEdgePaths(selectedClass);
classPathHits = sig.memberPathHits(selectedClass);
classHits = sig.hits(selectedClass);
classPilot = pilotAll(selectedClass);

%% Estimate one Poisson intensity vector for every day-type/time/date stratum
if o.EstimatorFastMode
    estDefaults = { ...
        'MaxCombinationSize',min(3,mEst), ...
        'MaxCandidates',max(mEst,o.EstimatorFastMaxCandidates), ...
        'NumProjections',[], ...
        'RequireGCDOne',true, ...
        'TrimToIdentifiable',true, ...
        'ComputeCandidateFisherRanks',false, ...
        'MaxIter',60, ...
        'MaxSubproblemIter',40, ...
        'TolGradient',1e-5, ...
        'CandidateSeed',o.Seed+101, ...
        'Display','off'};
else
    estDefaults = { ...
        'MaxCombinationSize',min(3,mEst), ...
        'MaxCandidates',max(2000,20*mEst), ...
        'NumProjections',[], ...
        'RequireGCDOne',true, ...
        'TrimToIdentifiable',false, ...
        'CandidateSeed',o.Seed+101, ...
        'Display','off'};
end
% Real counts are not exactly consistent with Y=AX; a small contamination
% weight keeps every observed projected count finite (see estimator help).
estDefaults = [estDefaults,{'Contamination',o.Contamination,'NumStarts',o.NumStarts,'MultiStartSeed',o.Seed+303, ...
    'MeanAnchorWeight',o.MeanAnchorWeight,'IncludeSingleSensors',logical(o.IncludeSingleSensors), ...
    'UseParallel',logical(o.UseParallel)}];
estOptions = mergeNameValueOptions(estDefaults,o.EstimatorOptions);

% With slice-specific missing-sensor exclusion, different strata may use
% different row subsets of A. Projection matrices therefore cannot blindly be
% shared across every slice. Reuse the expensive D-optimal design only among
% strata with the same sensor mask.
commonProjectionMatrix=[];          % populated only when every usable slice has one mask
commonProjectionDesign=struct();    % legacy compatibility for the one-mask case
commonProjectionDesignSeconds=0;
projectionMatrixByStratum=cell(numel(strata),1);
projectionDesignGroups=struct('SensorMask',{},'SensorRows',{},'Strata',{},'C',{},'Design',{},'Seconds',{},'Identifiable',{},'Reason',{},'ClassReduction',{});
identifiableByStratum=false(numel(strata),1); % reduced-model feasibility
originalClassesIdentifiableByStratum=false(numel(strata),1);
classReductionByStratum=cell(numel(strata),1);
reducedClassCountByStratum=zeros(numel(strata),1);
unobservedClassByStratum=false(p,numel(strata));
mergedOriginalClassCountByStratum=zeros(numel(strata),1);
identifiabilityReason=strings(numel(strata),1);
maskKeys=strings(numel(strata),1);
for q=1:numel(strata)
    mask=logical(strata(q).SensorMask(:)) & sectionEstimationMask(:);
    maskKeys(q)=string(char('0'+double(mask.')));
    if strcmp(o.SliceClassHandling,'merge')
        reduction=reduceSliceClassesLocal(A(mask,:));
    else
        % Identity mapping preserves original column order and original checks.
        reduction=struct('A',A(mask,:),'Expansion',speye(p), ...
            'GroupIndex',(1:p).','GroupSize',ones(p,1), ...
            'Unobserved',~any(A(mask,:)~=0,1).', ...
            'AllocationPolicy',"Strict: original classes; unobserved or duplicate signatures skip the stratum");
    end
    classReductionByStratum{q}=reduction;
    reducedClassCountByStratum(q)=size(reduction.A,2);
    unobservedClassByStratum(:,q)=reduction.Unobserved;
    mergedOriginalClassCountByStratum(q)=nnz(reduction.GroupSize>1);
    [originalClassesIdentifiableByStratum(q),~]=sliceMeasurementIdentifiableLocal(A(mask,:));
    if nnz(mask)<o.MinObservedSections
        identifiabilityReason(q)=sprintf('only %d complete sensors (< MinObservedSections=%d)',nnz(mask),o.MinObservedSections);
        continue
    end
    [identifiableByStratum(q),identifiabilityReason(q)]=sliceMeasurementIdentifiableLocal(reduction.A);
end

if o.ReuseProjectionDesign
    good=find(identifiableByStratum);
    uniqueKeys=unique(maskKeys(good),'stable');
    fprintf('Designing optimized projections for %d distinct complete-sensor mask(s) ...\n',numel(uniqueKeys));
    for ig=1:numel(uniqueKeys)
        idx=good(maskKeys(good)==uniqueKeys(ig));
        mask=logical(strata(idx(1)).SensorMask(:)) & sectionEstimationMask(:);
        reduction=classReductionByStratum{idx(1)};
        Aq=reduction.A;
        YdesignParts=cell(numel(idx),1);
        for kk=1:numel(idx)
            Yk=strata(idx(kk)).YFull(:,mask);
            if isfield(nativeInfo,'outerTestDates') && numel(strata(idx(kk)).Dates)==size(Yk,1)
                Yk=Yk(~ismember(strata(idx(kk)).Dates(:),nativeInfo.outerTestDates),:);   % non-test rows
            end
            YdesignParts{kk}=Yk;
        end
        Ydesign=vertcat(YdesignParts{:});
        pooledPilot=max(lsqnonneg(Aq,mean(Ydesign,1).',optimset('Display','off')),0.05);
        designOpts=mergeNameValueOptions(estOptions,{ ...
            'PilotMu',pooledPilot,'Mu0',pooledPilot,'DesignOnly',true, ...
            'CandidateSeed',o.Seed+101+ig});
        fprintf('  mask %d/%d: %d sensors, %d strata ... ',ig,numel(uniqueKeys),nnz(mask),numel(idx));
        tDesign=tic;
        [~,dsgn]=estimateMuOptimizedProjections(Aq,Ydesign,designOpts{:});
        sec=toc(tDesign); commonProjectionDesignSeconds=commonProjectionDesignSeconds+sec;
        C=dsgn.C;
        for qq=reshape(idx,1,[]), projectionMatrixByStratum{qq}=C; end
        projectionDesignGroups(end+1)=struct('SensorMask',mask,'SensorRows',find(mask), ...
            'Strata',idx(:),'C',C,'Design',dsgn,'Seconds',sec,'Identifiable',true,'Reason',"", ...
            'ClassReduction',reduction); %#ok<AGROW>
        fprintf('%d projections in %.2f s.\n',size(C,1),sec);
    end
    if numel(projectionDesignGroups)==1 && nnz(identifiableByStratum)==numel(strata)
        commonProjectionMatrix=projectionDesignGroups(1).C;
        commonProjectionDesign=projectionDesignGroups(1).Design;
    end
end

Snum=numel(strata);
muHatByStratum=nan(p,Snum);
muRatePerHourByStratum=nan(p,Snum);
infoByStratum=cell(Snum,1);
trainMeanRelErr=nan(Snum,1);
testMeanRelErr=nan(Snum,1);
allMeanRelErr=nan(Snum,1);
firstOrderOptimality=nan(Snum,1);
estimationTime=nan(Snum,1);
nTrainVec=zeros(Snum,1);
nTestVec=zeros(Snum,1);
testDaysByStratum=zeros(Snum,1);     % distinct test days in the stratum
medianFanoByStratum=nan(Snum,1);

% Measured/predicted mean-flow diagnostics.  Each column is one
% day-type/time/date stratum and each row is one measured road section.
predMeanByStratum=nan(m,Snum);
trainMeasuredMeanByStratum=nan(m,Snum);
testMeasuredMeanByStratum=nan(m,Snum);
allMeasuredMeanByStratum=nan(m,Snum);
sectionValidationMeasuredMeanByStratum=nan(m,Snum);
sectionValidationPredictedMeanByStratum=nan(m,Snum);
sectionValidationObservationCountByStratum=zeros(m,Snum);
% Diagnostics added for the paper (baseline, Poisson checks, identifiability).
sectionValidationPredictedMeanNNLSByStratum=nan(m,Snum); % NNLS pilot (means only)
lag1AutocorrByStratum=nan(Snum,1);       % median lag-one autocorrelation over sensors
projectionRankByStratum=nan(Snum,1);     % jump-class rank of the projections used
rankConditionByStratum=false(Snum,1);    % rank == number of reduced classes
reducedClassCountFitted=nan(Snum,1);     % number of reduced classes in the fit
sensitivityConditionByStratum=nan(Snum,1); % condition number of the criterion Hessian
medianRelativeSEByStratum=nan(Snum,1);   % median sandwich SE / estimate (interior classes)
trainMisfitByStratum=nan(Snum,1);        % ||A mu - ybar||/||ybar|| of the estimate
nnlsMisfitByStratum=nan(Snum,1);         % same for the NNLS fit (best the classes allow)
exitflagByStratum=nan(Snum,1);           % trust-region exit flag (0: iteration limit)
downweightedFractionByStratum=nan(Snum,1); % share of projected counts downweighted
multistartAtBestByStratum=nan(Snum,1);   % share of starts reaching the best maximum
multistartDistinctByStratum=nan(Snum,1); % number of distinct local maxima found
multistartMaxRelDistByStratum=nan(Snum,1); % largest ||mu_k-muHat||/||muHat||
multistartGapByStratum=nan(Snum,1);      % criterion gap per observation, best - worst start
reliableByStratum=false(Snum,1);         % fitted and passing all reliability checks
reliabilityReason=strings(Snum,1);       % why a stratum is not reliable
muHatUnreliableByStratum=nan(p,Snum);    % estimates of rejected fits (audit only)

trainMAE=nan(Snum,1); trainRMSE=nan(Snum,1); trainNRMSE=nan(Snum,1);
trainBias=nan(Snum,1); trainBiasPercent=nan(Snum,1);
testMAE=nan(Snum,1); testRMSE=nan(Snum,1); testNRMSE=nan(Snum,1);
testBias=nan(Snum,1); testBiasPercent=nan(Snum,1);
allMAE=nan(Snum,1); allRMSE=nan(Snum,1); allNRMSE=nan(Snum,1);
allBias=nan(Snum,1); allBiasPercent=nan(Snum,1);
sectionValidationMAE=nan(Snum,1); sectionValidationRMSE=nan(Snum,1); sectionValidationNRMSE=nan(Snum,1);
sectionValidationBias=nan(Snum,1); sectionValidationBiasPercent=nan(Snum,1);

fprintf('\nEstimating weekday/time/date/year-specific Poisson intensities...\n');
sensorCompleteByStratum=false(m,Snum);
sensorUsedByStratum=false(m,Snum);
predMeanNNLSByStratum=nan(m,Snum);       % NNLS prediction of all sensors (baseline)
nnlsRelResidualByStratum=nan(m,Snum);    % (predicted-observed)/observed of the stratum NNLS fit
% Strata are independent: fit them with parfor when UseParallel is set (the
% estimator then runs serially inside each worker). Results are identical.
ctx=struct();
ctx.A=A;
ctx.Aq=Aq;
ctx.classReductionByStratum=classReductionByStratum;
ctx.estOptions=estOptions;
ctx.fq=fq;
ctx.i=i;
ctx.identifiableByStratum=identifiableByStratum;
ctx.m=m;
ctx.mEst=mEst;
ctx.mask=mask;
ctx.nativeInfo=nativeInfo;
ctx.o=o;
ctx.p=p;
ctx.projectionMatrixByStratum=projectionMatrixByStratum;
ctx.reduction=reduction;
ctx.sectionEstimationMask=sectionEstimationMask;
ctx.sectionValidationMask=sectionValidationMask;
parWorkers=0; if o.UseParallel, parWorkers=Inf; end
if o.UseParallel
    ctx.estOptions=mergeNameValueOptions(estOptions,{'UseParallel',false});
end
stratumOut=cell(Snum,1);
parfor (q=1:Snum, parWorkers)
    stratumOut{q}=fitStratumLocal(q,strata(q),ctx);
end
for q=1:Snum
    so=stratumOut{q};
    muHatByStratum(:,q)=so.muHatByStratum(:,1);
    muRatePerHourByStratum(:,q)=so.muRatePerHourByStratum(:,1);
    infoByStratum{q}=so.infoByStratum{1};
    trainMeanRelErr(q)=so.trainMeanRelErr(1);
    testMeanRelErr(q)=so.testMeanRelErr(1);
    allMeanRelErr(q)=so.allMeanRelErr(1);
    firstOrderOptimality(q)=so.firstOrderOptimality(1);
    estimationTime(q)=so.estimationTime(1);
    nTrainVec(q)=so.nTrainVec(1);
    nTestVec(q)=so.nTestVec(1);
    testDaysByStratum(q)=so.testDaysByStratum(1);
    medianFanoByStratum(q)=so.medianFanoByStratum(1);
    predMeanByStratum(:,q)=so.predMeanByStratum(:,1);
    trainMeasuredMeanByStratum(:,q)=so.trainMeasuredMeanByStratum(:,1);
    testMeasuredMeanByStratum(:,q)=so.testMeasuredMeanByStratum(:,1);
    allMeasuredMeanByStratum(:,q)=so.allMeasuredMeanByStratum(:,1);
    sectionValidationMeasuredMeanByStratum(:,q)=so.sectionValidationMeasuredMeanByStratum(:,1);
    sectionValidationPredictedMeanByStratum(:,q)=so.sectionValidationPredictedMeanByStratum(:,1);
    sectionValidationObservationCountByStratum(:,q)=so.sectionValidationObservationCountByStratum(:,1);
    sectionValidationPredictedMeanNNLSByStratum(:,q)=so.sectionValidationPredictedMeanNNLSByStratum(:,1);
    lag1AutocorrByStratum(q)=so.lag1AutocorrByStratum(1);
    projectionRankByStratum(q)=so.projectionRankByStratum(1);
    rankConditionByStratum(q)=so.rankConditionByStratum(1);
    reducedClassCountFitted(q)=so.reducedClassCountFitted(1);
    sensitivityConditionByStratum(q)=so.sensitivityConditionByStratum(1);
    medianRelativeSEByStratum(q)=so.medianRelativeSEByStratum(1);
    trainMisfitByStratum(q)=so.trainMisfitByStratum(1);
    nnlsMisfitByStratum(q)=so.nnlsMisfitByStratum(1);
    exitflagByStratum(q)=so.exitflagByStratum(1);
    downweightedFractionByStratum(q)=so.downweightedFractionByStratum(1);
    multistartAtBestByStratum(q)=so.multistartAtBestByStratum(1);
    multistartDistinctByStratum(q)=so.multistartDistinctByStratum(1);
    multistartMaxRelDistByStratum(q)=so.multistartMaxRelDistByStratum(1);
    multistartGapByStratum(q)=so.multistartGapByStratum(1);
    reliableByStratum(q)=so.reliableByStratum(1);
    reliabilityReason(q)=so.reliabilityReason(1);
    muHatUnreliableByStratum(:,q)=so.muHatUnreliableByStratum(:,1);
    trainMAE(q)=so.trainMAE(1);
    trainRMSE(q)=so.trainRMSE(1);
    trainNRMSE(q)=so.trainNRMSE(1);
    trainBias(q)=so.trainBias(1);
    trainBiasPercent(q)=so.trainBiasPercent(1);
    testMAE(q)=so.testMAE(1);
    testRMSE(q)=so.testRMSE(1);
    testNRMSE(q)=so.testNRMSE(1);
    testBias(q)=so.testBias(1);
    testBiasPercent(q)=so.testBiasPercent(1);
    allMAE(q)=so.allMAE(1);
    allRMSE(q)=so.allRMSE(1);
    allNRMSE(q)=so.allNRMSE(1);
    allBias(q)=so.allBias(1);
    allBiasPercent(q)=so.allBiasPercent(1);
    sectionValidationMAE(q)=so.sectionValidationMAE(1);
    sectionValidationRMSE(q)=so.sectionValidationRMSE(1);
    sectionValidationNRMSE(q)=so.sectionValidationNRMSE(1);
    sectionValidationBias(q)=so.sectionValidationBias(1);
    sectionValidationBiasPercent(q)=so.sectionValidationBiasPercent(1);
    sensorCompleteByStratum(:,q)=so.sensorCompleteByStratum(:,1);
    sensorUsedByStratum(:,q)=so.sensorUsedByStratum(:,1);
    predMeanNNLSByStratum(:,q)=so.predMeanNNLSByStratum(:,1);
    nnlsRelResidualByStratum(:,q)=so.nnlsRelResidualByStratum(:,1);
end

% ---- Section audit: where measurement sections are lost --------------------
sectionAudit=struct();
sectionAudit.preselection=table(statusCats(:), ...
    arrayfun(@(c)nnz(matchExport.PreselectionStatus==c),statusCats(:)), ...
    'VariableNames',{'Status','Sections'});
sectionAudit.analysed=m; sectionAudit.estimation=mEst; sectionAudit.validation=mVal;
completeEst=sum(sensorCompleteByStratum(sectionEstimationMask,:),1);
usedEst=sum(sensorUsedByStratum,1);
sectionAudit.completePerStratum=completeEst(:);
sectionAudit.usedPerStratum=usedEst(:);
fracComplete=mean(sensorCompleteByStratum,2);
fracUsed=mean(sensorUsedByStratum,2);
roleNames=repmat("estimation",m,1); roleNames(sectionValidationMask)="validation";
sectionAudit.perSection=table((1:m).',selectedSections(:),roleNames,fracComplete,fracUsed, ...
    'VariableNames',{'Measurement','MappingRow','Role','ShareOfStrataComplete','ShareOfStrataUsed'});
fprintf('\nSection audit\n');
for is=1:height(sectionAudit.preselection)
    fprintf('  preselection %-30s %d\n',char(sectionAudit.preselection.Status(is)), ...
        sectionAudit.preselection.Sections(is));
end
fprintf('  analysed measurements: %d (estimation %d, validation %d)\n',m,mEst,mVal);
fprintf(['  estimation sensors complete in a stratum (coverage >= %.0f%%): ' ...
    'median %g, range %g-%g of %d\n'],100*o.SliceSensorMinCoverage, ...
    median(completeEst),min(completeEst),max(completeEst),mEst);
fprintf('  estimation sensors never complete in any stratum: %d\n', ...
    nnz(fracComplete(sectionEstimationMask)==0));
if ~isempty(o.SectionAuditCSV)
    try
        writetable(sectionAudit.perSection,char(o.SectionAuditCSV));
        fprintf('  per-section audit written to %s (join MappingRow with %s)\n', ...
            char(o.SectionAuditCSV),char(o.MappingCSV));
    catch MEaudit
        warning('runParmaTrafficTomography:SectionAudit','Could not write the section audit: %s',MEaudit.message);
    end
end

% ---- Largest test errors, classified ---------------------------------------------
validationOutliers=table();
try
    vOK=all(isfinite(muHatByStratum),1);
    Mtr=trainMeasuredMeanByStratum; Mte=testMeasuredMeanByStratum; Pr=predMeanByStratum;
    use=sensorUsedByStratum & repmat(vOK,m,1) & isfinite(Mtr) & isfinite(Mte) & isfinite(Pr);
    [ii,qq]=find(use);
    if ~isempty(ii)
        lin=sub2ind(size(Mtr),ii,qq);
        perHour=60./[strata(qq).ObservationMinutes].';
        errH=(Pr(lin)-Mte(lin)).*perHour;
        shift=abs(Mte(lin)-Mtr(lin))./max(Mtr(lin),1);     % test vs training data
        misfit=abs(Pr(lin)-Mtr(lin))./max(Mtr(lin),1);     % prediction vs training data
        cause=repmat("other",numel(lin),1);
        cause(shift>0.5 & misfit<=0.2)="data change on test days";
        cause(misfit>0.5)="fit error at the sensor";
        [~,ord]=sort(abs(errH),'descend'); top=ord(1:min(50,numel(ord)));
        names=strtrim(string(match.SectionName(selectedSections(ii(top)))));
        tDays=testDaysByStratum(qq);
        validationOutliers=table(ii(top),names,qq(top),string({strata(qq(top)).DayType}).',tDays(top), ...
            Mtr(lin(top)).*perHour(top),Mte(lin(top)).*perHour(top),Pr(lin(top)).*perHour(top),errH(top),cause(top), ...
            'VariableNames',{'Measurement','SectionName','Stratum','DayType','TestDays','TrainMeanPerHour', ...
            'TestMeanPerHour','PredictedPerHour','ErrorPerHour','LikelyCause'});
        k10=ord(1:max(1,ceil(0.01*numel(ord))));             % largest 1% of the errors
        fprintf(['\nLargest 1%% of test errors (%d sensor-stratum pairs): %.0f%% data change on test days, ' ...
            '%.0f%% fit error at the sensor, %.0f%% other.\n'],numel(k10), ...
            100*mean(cause(k10)=="data change on test days"),100*mean(cause(k10)=="fit error at the sensor"), ...
            100*mean(cause(k10)=="other"));
        fprintf(['  Among them, %.0f%% come from strata with a single test day (%.0f%% of all pairs do); ' ...
            'median test days per stratum: %g.\n'],100*mean(tDays(k10)==1),100*mean(tDays==1),median(testDaysByStratum(testDaysByStratum>0)));
        for k=1:min(10,numel(top))
            r_=validationOutliers(k,:);
            fprintf('  %-40s %-15s train %7.0f test %7.0f pred %7.0f veh/h  -> %s\n',char(r_.SectionName), ...
                char(r_.DayType),r_.TrainMeanPerHour,r_.TestMeanPerHour,r_.PredictedPerHour,char(r_.LikelyCause));
        end
        if ~isempty(o.ValidationOutlierCSV), writetable(validationOutliers,char(o.ValidationOutlierCSV)); end
    end
catch MEout
    warning('runParmaTrafficTomography:ValidationOutliers','Validation outlier report skipped: %s',MEout.message);
end

% ---- Validation bias summary -----------------------------------------------------
if any(sectionValidationMask)
    vm=sectionValidationMeasuredMeanByStratum(sectionValidationMask,:);
    vp_=sectionValidationPredictedMeanByStratum(sectionValidationMask,:);
    okv=isfinite(vm)&isfinite(vp_)&vm>0&repmat(all(isfinite(muHatByStratum),1),size(vm,1),1);
    if any(okv(:))
        rel=(vp_(okv)-vm(okv))./vm(okv);
        big=vm(okv)>=median(vm(okv));
        fprintf(['\nValidation (held-out sensors, reliable strata): median signed local error %+.0f%%, ' ...
            'under-predicted pairs %.0f%%, predicted ~0 (< 5%% of observed) %.0f%%; ' ...
            'median |error| %.0f%% on the busier half of the pairs.\n'], ...
            100*median(rel),100*mean(rel<0),100*mean(vp_(okv)<0.05*vm(okv)),100*median(abs(rel(big))));
    end
end

% ---- Why strata have no estimate (counts by cause) ----------------------------
causes={'rank','jump-class rank'; 'conditioning','Hessian condition'; ...
        'mean misfit','misfit'; 'too few intervals per class','training intervals <'; ...
        'too few intervals','MinSamplesPerStratum'; 'too few sensors','complete estimation sensors'; ...
        'not identifiable','not identifiable'; 'estimator error','estimator error'};
causeCount=zeros(size(causes,1),1);
for q=1:Snum
    for c=1:size(causes,1)
        causeCount(c)=causeCount(c)+contains(reliabilityReason(q),causes{c,2});
    end
end
fprintf('\nStrata with a reliable estimate: %d of %d. Causes for the others (a stratum can have several):\n', ...
    nnz(reliableByStratum),Snum);
for c=1:size(causes,1)
    if causeCount(c)>0, fprintf('  %-30s %d\n',causes{c,1},causeCount(c)); end
end
fittedCounts=reducedClassCountFitted(isfinite(reducedClassCountFitted));
if ~isempty(fittedCounts)
    fprintf('  classes fitted per stratum: median %g, range %g-%g\n',median(fittedCounts),min(fittedCounts),max(fittedCounts));
end
sectionAudit.rejectionCauses=table(causes(:,1),causeCount,'VariableNames',{'Cause','Strata'});

% ---- Per-sensor NNLS residual report -------------------------------------------
roleRes=repmat("estimation",m,1); roleRes(sectionValidationMask)="validation";
latS=nan(m,1); lonS=nan(m,1); nameS=strings(m,1);
try, latS=double(match.Lat(selectedSections)); lonS=double(match.Lon(selectedSections)); catch, end
try, nameS=string(match.SectionName(selectedSections)); catch, end
medStrata=median(nnlsRelResidualByStratum,2,'omitnan');
sensorResiduals=table((1:m).',selectedSections(:),nameS,roleRes,latS,lonS,meanMeasured(:), ...
    predGlobalNNLS(:),relResGlobalNNLS(:),medStrata,mean(sensorUsedByStratum,2), ...
    'VariableNames',{'Measurement','MappingRow','SectionName','Role','Lat','Lon','MeanCount', ...
    'PredictedDictionaryNNLS','RelResidualDictionaryNNLS','MedianRelResidualStrataNNLS','ShareOfStrataUsed'});
[~,ordRes]=sort(abs(sensorResiduals.RelResidualDictionaryNNLS),'descend','MissingPlacement','last');
fprintf('\nSensors with the largest relative residual of the dictionary-wide NNLS fit (predicted-observed)/observed:\n');
for k=ordRes(1:min(10,m)).'
    fprintf('  %4d %-40s %-10s mean %8.1f  residual %+6.2f  (strata median %+6.2f)\n',k, ...
        char(sensorResiduals.SectionName(k)),char(sensorResiduals.Role(k)),sensorResiduals.MeanCount(k), ...
        sensorResiduals.RelResidualDictionaryNNLS(k),sensorResiduals.MedianRelResidualStrataNNLS(k));
end
if ~isempty(o.SensorResidualCSV)
    try
        writetable(sensorResiduals,char(o.SensorResidualCSV));
        fprintf('  per-sensor residuals written to %s\n',char(o.SensorResidualCSV));
    catch MEres
        warning('runParmaTrafficTomography:SensorResiduals','Could not write the residual table: %s',MEres.message);
    end
end

validFit=all(isfinite(muHatByStratum),1).';
if ~any(validFit)
    lines=strings(Snum,1);
    for q=1:Snum
        r=reliabilityReason(q); if strlength(r)==0, r="no estimate"; end
        lines(q)=sprintf('  stratum %d (%s %s-%s): %s',q,strata(q).DayType, ...
            minuteLabel(strata(q).StartMinute),minuteLabel(strata(q).EndMinute),r);
    end
    nm=nnlsMisfitByStratum(isfinite(nnlsMisfitByStratum));
    hint='';
    if ~isempty(nm)
        hint=sprintf(['\nMedian training-mean misfit of the NNLS fit over the fitted strata: %.3g ' ...
            '(dictionary-wide NNLS residual: %.3g). If this is large, the retained candidate ' ...
            'paths cannot reproduce the sensor means: check the sensor-to-edge matching and ' ...
            'directions, and increase NumPathSamples or MaxPathClasses.'],median(nm),pilotFit);
    end
    error('runParmaTrafficTomography:AllStrataFailed', ...
        ['No stratum has a reliable estimate:\n%s%s\n' ...
         'Relax the checks if appropriate, e.g. ''MaxSensitivityCondition'', ' ...
         '''MaxExcessMisfit'', ''MinSamplesPerClass'', or ''RequireRankCondition'',false.'], ...
        strjoin(lines,newline),hint);
end

% Comparable long-run summary in vehicles/hour. Weight by number of observations
% in each successfully estimated regime, not by regime duration.
w=double(arrayfun(@(x)size(x.Y,1),strata));
w=w(:);                    % force a column vector, independent of strata orientation
w(~validFit(:))=0;
wValid=w(validFit(:));
% Orientation-safe weighted average across successfully fitted strata:
%   (p x Svalid) * (Svalid x 1) -> (p x 1)
meanHourlyMu=(muRatePerHourByStratum(:,validFit)*wValid) ./ max(sum(wValid),eps);

%% Tables
signatureText = strings(p,1);
originLat = nan(p,1); originLon = nan(p,1); destLat = nan(p,1); destLon = nan(p,1);
for j=1:p
    signatureText(j)=char('0'+A(:,j).');
    pn=repNodePaths{j};
    if ~isempty(pn)
        originLat(j)=G.Nodes.Lat(pn(1)); originLon(j)=G.Nodes.Lon(pn(1));
        destLat(j)=G.Nodes.Lat(pn(end)); destLon(j)=G.Nodes.Lon(pn(end));
    end
end
classDefinitionTable=table((1:p).',meanHourlyMu,classPilot(:),classHits(:),signatureText, ...
    originLat,originLon,destLat,destLon, ...
    'VariableNames',{'Class','MeanEstimatedVehiclesPerHour','ScreeningNNLSPilot','SampledPathHits', ...
                     'MeasuredSignature','RepresentativeOriginLat','RepresentativeOriginLon', ...
                     'RepresentativeDestLat','RepresentativeDestLon'});

% Explicit stratum coordinates: day type x time span x calendar date span.
Stratum=(1:Snum).';
DayTypeCode=reshape([strata.DayTypeCode],[],1);
DayType=string({strata.DayType}).'; DayType=DayType(:);
Year=reshape([strata.Year],[],1);
YearLabel=string({strata.YearLabel}).'; YearLabel=YearLabel(:);
StartMinute=reshape([strata.StartMinute],[],1);
EndMinute=reshape([strata.EndMinute],[],1);
DurationMinutes=reshape([strata.DurationMinutes],[],1);
ObservationMinutes=reshape([strata.ObservationMinutes],[],1);
SeasonStartDOY=reshape([strata.SeasonStartDOY],[],1);
SeasonEndDOY=reshape([strata.SeasonEndDOY],[],1);
SeasonStart=string({strata.SeasonStart}).'; SeasonStart=SeasonStart(:);
SeasonEnd=string({strata.SeasonEnd}).'; SeasonEnd=SeasonEnd(:);
SeasonLabel=string({strata.SeasonLabel}).'; SeasonLabel=SeasonLabel(:);
CandidateDates=reshape([strata.CandidateDates],[],1);
CompleteDates=reshape([strata.CompleteDates],[],1);
CandidateObservations=reshape([strata.CandidateObservations],[],1);
CompleteObservations=reshape([strata.CompleteObservations],[],1);
RefinementDepth=reshape([strata.RefinementDepth],[],1);
RefinedForFano=reshape([strata.RefinedForFano],[],1);
OutlierValuesRemoved=reshape([strata.OutlierValuesRemoved],[],1);
OutlierObservationRowsAffected=reshape([strata.OutlierObservationRowsAffected],[],1);
OutlierFractionOfCandidateValues=reshape([strata.OutlierFractionOfCandidateValues],[],1);
CompleteSensors=arrayfun(@(x)nnz(x.SensorMask),strata).'; CompleteSensors=CompleteSensors(:);
SensorsExcludedMissing=mEst-CompleteSensors;
RowRetentionFraction=reshape([strata.RowRetentionFraction],[],1);
RowsRetainedAfterMissingFilter=reshape([strata.RowsRetainedAfterMissingFilter],[],1);
MedianEstimationSensorCoverage=nan(Snum,1);
for qCov=1:Snum
    covq=double(strata(qCov).SensorCoverage(:));
    covq=covq(sectionEstimationMask(:));
    MedianEstimationSensorCoverage(qCov)=median(covq,'omitnan');
end
SensorsUsed=sum(sensorUsedByStratum,1).';
CompleteSensorRows=strings(Snum,1); SensorRowsUsed=strings(Snum,1);
for qMask=1:Snum
    CompleteSensorRows(qMask)=strjoin(string(find(strata(qMask).SensorMask(:)).'),',');
    SensorRowsUsed(qMask)=strjoin(string(find(sensorUsedByStratum(:,qMask)).'),',');
end
AllSensorCompleteObservations=reshape([strata.AllSensorCompleteObservations],[],1);
N=arrayfun(@(x)size(x.Y,1),strata); N=N(:);

StartTime=strings(Snum,1);
EndTime=strings(Snum,1);
for qTable=1:Snum
    StartTime(qTable)=string(minuteLabel(StartMinute(qTable)));
    EndTime(qTable)=string(minuteLabel(EndMinute(qTable)));
end

nTrainVec=nTrainVec(:); nTestVec=nTestVec(:);
medianFanoByStratum=medianFanoByStratum(:);
trainMeanRelErr=trainMeanRelErr(:); testMeanRelErr=testMeanRelErr(:); allMeanRelErr=allMeanRelErr(:);
firstOrderOptimality=firstOrderOptimality(:); estimationTime=estimationTime(:);
EstimationSensorsUsed=sum(sensorUsedByStratum,1).';
ValidationSensorsAvailable=sum(sectionValidationObservationCountByStratum>0,1).';
sectionValidationMAE=sectionValidationMAE(:); sectionValidationRMSE=sectionValidationRMSE(:); sectionValidationNRMSE=sectionValidationNRMSE(:);
sectionValidationBias=sectionValidationBias(:); sectionValidationBiasPercent=sectionValidationBiasPercent(:);
trainMAE=trainMAE(:); trainRMSE=trainRMSE(:); trainNRMSE=trainNRMSE(:);
trainBias=trainBias(:); trainBiasPercent=trainBiasPercent(:);
testMAE=testMAE(:); testRMSE=testRMSE(:); testNRMSE=testNRMSE(:);
testBias=testBias(:); testBiasPercent=testBiasPercent(:);
allMAE=allMAE(:); allRMSE=allRMSE(:); allNRMSE=allNRMSE(:);
allBias=allBias(:); allBiasPercent=allBiasPercent(:);
Estimated=logical(validFit(:));

strataTable=table(Stratum,DayTypeCode,DayType,Year,YearLabel,StartMinute,EndMinute,StartTime,EndTime, ...
    SeasonStartDOY,SeasonEndDOY,SeasonStart,SeasonEnd,SeasonLabel,DurationMinutes,ObservationMinutes, ...
    CandidateDates,CompleteDates,CandidateObservations,CompleteObservations,AllSensorCompleteObservations,CompleteSensors,SensorsExcludedMissing,MedianEstimationSensorCoverage,RowsRetainedAfterMissingFilter,RowRetentionFraction,CompleteSensorRows,SensorsUsed,SensorRowsUsed,RefinementDepth,RefinedForFano, ...
    OutlierValuesRemoved,OutlierObservationRowsAffected,OutlierFractionOfCandidateValues,N,nTrainVec,nTestVec,medianFanoByStratum, ...
    trainMeanRelErr,testMeanRelErr,allMeanRelErr, ...
    trainMAE,trainRMSE,trainNRMSE,trainBias,trainBiasPercent, ...
    testMAE,testRMSE,testNRMSE,testBias,testBiasPercent, ...
    allMAE,allRMSE,allNRMSE,allBias,allBiasPercent, ...
    EstimationSensorsUsed,ValidationSensorsAvailable,sectionValidationMAE,sectionValidationRMSE,sectionValidationNRMSE,sectionValidationBias,sectionValidationBiasPercent, ...
    firstOrderOptimality,estimationTime,Estimated, ...
    'VariableNames',{'Stratum','DayTypeCode','DayType','Year','YearLabel','StartMinute','EndMinute','StartTime','EndTime', ...
    'SeasonStartDOY','SeasonEndDOY','SeasonStart','SeasonEnd','SeasonLabel','DurationMinutes','ObservationMinutes', ...
    'CandidateDates','CompleteDates','CandidateObservations','CompleteObservations','AllSensorCompleteObservations','CompleteSensors','SensorsExcludedMissing','MedianEstimationSensorCoverage','RowsRetainedAfterMissingFilter','RowRetentionFraction','CompleteSensorRows','SensorsUsed','SensorRowsUsed','RefinementDepth','RefinedForFano', ...
    'OutlierValuesRemoved','OutlierObservationRowsAffected','OutlierFractionOfCandidateValues','N','NTrain','NTest','MedianFano', ...
    'TrainMeanRelativeError','TestMeanRelativeError','AllMeanRelativeError', ...
    'TrainMAE','TrainRMSE','TrainNRMSE','TrainBias','TrainBiasPercent', ...
    'TestMAE','TestRMSE','TestNRMSE','TestBias','TestBiasPercent', ...
    'AllMAE','AllRMSE','AllNRMSE','AllBias','AllBiasPercent', ...
    'EstimationSensorsUsed','ValidationSensorsAvailable','SectionValidationMAE','SectionValidationRMSE','SectionValidationNRMSE','SectionValidationBias','SectionValidationBiasPercent', ...
    'FirstOrderOptimality','EstimationTimeSeconds','Estimated'});
strataTable.ReducedClassCount=reducedClassCountByStratum;
strataTable.UnobservedClassCount=sum(unobservedClassByStratum,1).';
strataTable.MergedOriginalClassCount=mergedOriginalClassCountByStratum;
strataTable.OriginalClassesIdentifiable=originalClassesIdentifiableByStratum;
writetable(strataTable,o.StrataCSV);
fprintf('Poisson strata table written to: %s\n',o.StrataCSV);

% Long-form class x stratum table.
rows=Snum*p;
Class=repmat((1:p).',Snum,1);
StratumLong=kron((1:Snum).',ones(p,1));
DayTypeLong=repelem(DayType,p,1);
YearLong=repelem(Year,p,1);
YearLabelLong=repelem(YearLabel,p,1);
StartTimeLong=repelem(StartTime,p,1);
EndTimeLong=repelem(EndTime,p,1);
SeasonStartLong=repelem(SeasonStart,p,1);
SeasonEndLong=repelem(SeasonEnd,p,1);
MuHat=reshape(muHatByStratum,p*Snum,1);
MuPerHour=reshape(muRatePerHourByStratum,p*Snum,1);
Signature=repmat(signatureText,Snum,1);
classTable=table(Class,StratumLong,DayTypeLong,YearLong,YearLabelLong,StartTimeLong,EndTimeLong, ...
    SeasonStartLong,SeasonEndLong,MuHat,MuPerHour,Signature, ...
    'VariableNames',{'Class','Stratum','DayType','Year','YearLabel','StartTime','EndTime', ...
    'SeasonStart','SeasonEnd','MuHat','MuPerHour','MeasuredSignature'});
% Report which original-class values are allocated or unobserved conventions.
ReducedClassLong=zeros(p,Snum); GroupSizeLong=zeros(p,Snum);
for qRed=1:Snum
    ReducedClassLong(:,qRed)=classReductionByStratum{qRed}.GroupIndex;
    GroupSizeLong(:,qRed)=classReductionByStratum{qRed}.GroupSize;
end
classTable.ReducedClass=ReducedClassLong(:);
classTable.MergedGroupSize=GroupSizeLong(:);
classTable.UnobservedInStratum=unobservedClassByStratum(:);
classTable.EqualShareAllocation=GroupSizeLong(:)>1;
assert(height(classTable)==rows);
writetable(classTable,o.ClassCSV);
fprintf('Path-signature intensity table written to: %s\n',o.ClassCSV);

%% Measured-versus-predicted flow errors
% Long-form sensor x stratum table. There are two distinct validation layers:
% (1) temporal hold-out within ESTIMATION sections (YTrain/YTest), and
% (2) true SECTION hold-out rows that never enter the estimator at all.
nFlowRows=m*Snum;
MeasuredRow=repmat((1:m).',Snum,1);
SectionIndexLong=repmat(selectedSections(:),Snum,1);
SectionIDLong=repmat(match.SectionID(selectedSections),Snum,1);
SectionNameLong=repmat(match.SectionName(selectedSections),Snum,1);
OSMRoadLong=repmat(match.RoadName(selectedSections),Snum,1);
LaneCountLong=repmat(double(match.LaneCount(selectedSections)),Snum,1);
LaneMembersLong=repmat(string(match.LaneMemberSectionIDs(selectedSections)),Snum,1);
StratumFlow=kron((1:Snum).',ones(m,1));
DayTypeFlow=repelem(DayType,m,1);
YearFlow=repelem(Year,m,1);
YearLabelFlow=repelem(YearLabel,m,1);
StartTimeFlow=repelem(StartTime,m,1);
EndTimeFlow=repelem(EndTime,m,1);
SeasonStartFlow=repelem(SeasonStart,m,1);
SeasonEndFlow=repelem(SeasonEnd,m,1);
DurationFlow=repelem(DurationMinutes,m,1);
ObservationFlow=repelem(ObservationMinutes,m,1);
NTrainFlow=repelem(nTrainVec,m,1);
NTestFlow=repelem(nTestVec,m,1);
SectionRole=repmat(strings(m,1),Snum,1);
baseRole=repmat("Estimation",m,1); baseRole(sectionValidationMask)="Validation";
SectionRole=repmat(baseRole,Snum,1);
NSectionValidationSamples=reshape(sectionValidationObservationCountByStratum,m*Snum,1);
MeasuredTrainMean=reshape(trainMeasuredMeanByStratum,m*Snum,1);
MeasuredTestMean=reshape(testMeasuredMeanByStratum,m*Snum,1);
MeasuredAllMean=reshape(allMeasuredMeanByStratum,m*Snum,1);
PredictedMean=reshape(predMeanByStratum,m*Snum,1);
SectionValidationMeasuredMean=reshape(sectionValidationMeasuredMeanByStratum,m*Snum,1);
SectionValidationPredictedMean=reshape(sectionValidationPredictedMeanByStratum,m*Snum,1);
SectionValidationError=SectionValidationPredictedMean-SectionValidationMeasuredMean;
SectionValidationAbsoluteError=abs(SectionValidationError);
SectionValidationPercentError=100*SectionValidationError./max(abs(SectionValidationMeasuredMean),eps);
UsedInEstimation=reshape(sensorUsedByStratum,m*Snum,1);
TrainError=PredictedMean-MeasuredTrainMean;
TestError=PredictedMean-MeasuredTestMean;
AllError=PredictedMean-MeasuredAllMean;
TestAbsoluteError=abs(TestError);
TestPercentError=100*TestError./max(abs(MeasuredTestMean),eps);
MeasuredTestPerHour=MeasuredTestMean./max(ObservationFlow/60,eps);
PredictedPerHour=PredictedMean./max(ObservationFlow/60,eps);
TestErrorPerHour=PredictedPerHour-MeasuredTestPerHour;
flowErrorTable=table(MeasuredRow,SectionIndexLong,SectionIDLong,SectionNameLong,LaneCountLong,LaneMembersLong,OSMRoadLong, ...
    StratumFlow,DayTypeFlow,YearFlow,YearLabelFlow,StartTimeFlow,EndTimeFlow,SeasonStartFlow,SeasonEndFlow,DurationFlow,ObservationFlow,NTrainFlow,NTestFlow,SectionRole,NSectionValidationSamples, ...
    MeasuredTrainMean,MeasuredTestMean,MeasuredAllMean,PredictedMean,UsedInEstimation, ...
    SectionValidationMeasuredMean,SectionValidationPredictedMean,SectionValidationError,SectionValidationAbsoluteError,SectionValidationPercentError, ...
    TrainError,TestError,AllError,TestAbsoluteError,TestPercentError, ...
    MeasuredTestPerHour,PredictedPerHour,TestErrorPerHour, ...
    'VariableNames',{'MeasuredRow','SectionIndex','SectionID','SectionName','LaneCount','LaneMemberSectionIDs','OSMRoad', ...
    'Stratum','DayType','Year','YearLabel','StartTime','EndTime','SeasonStart','SeasonEnd','DurationMinutes','ObservationMinutes','NTrain','NTest','SectionRole','NSectionValidationSamples', ...
    'MeasuredTrainMean','MeasuredTestMean','MeasuredAllMean','PredictedMean','UsedInEstimation', ...
    'SectionValidationMeasuredMean','SectionValidationPredictedMean','SectionValidationError','SectionValidationAbsoluteError','SectionValidationPercentError', ...
    'TrainError','TestError','AllError','TestAbsoluteError','TestPercentError', ...
    'MeasuredTestPerHour','PredictedPerHour','TestErrorPerHour'});
assert(height(flowErrorTable)==nFlowRows);
writetable(flowErrorTable,o.FlowErrorCSV);
fprintf('Measured/predicted flow errors written to: %s\n',o.FlowErrorCSV);
sectionValidationTable=flowErrorTable(flowErrorTable.SectionRole=="Validation", ...
    {'MeasuredRow','SectionIndex','SectionID','SectionName','LaneCount','LaneMemberSectionIDs','OSMRoad', ...
     'Stratum','DayType','Year','YearLabel','StartTime','EndTime','SeasonStart','SeasonEnd', ...
     'NSectionValidationSamples','SectionValidationMeasuredMean','SectionValidationPredictedMean', ...
     'SectionValidationError','SectionValidationAbsoluteError','SectionValidationPercentError'});
writetable(sectionValidationTable,o.SectionValidationCSV);
fprintf('Section hold-out validation table written to: %s\n',o.SectionValidationCSV);

% Sensor-level hold-out summary, weighted by the number of test repetitions in
% each successfully fitted day-type/time/date regime.
sensorMeasuredMean=nan(m,1); sensorPredictedMean=nan(m,1);
sensorMAE=nan(m,1); sensorRMSE=nan(m,1); sensorNRMSE=nan(m,1);
sensorBias=nan(m,1); sensorBiasPercent=nan(m,1); sensorMaxAbsError=nan(m,1);
for iSensor=1:m
    ok=validFit(:).' & isfinite(testMeasuredMeanByStratum(iSensor,:)) & ...
        isfinite(predMeanByStratum(iSensor,:));
    if ~any(ok), continue; end
    ws=nTestVec(ok); ws=double(ws(:));
    obs=testMeasuredMeanByStratum(iSensor,ok).';
    pred=predMeanByStratum(iSensor,ok).';
    sw=max(sum(ws),eps);
    err=pred-obs;
    sensorMeasuredMean(iSensor)=sum(ws.*obs)/sw;
    sensorPredictedMean(iSensor)=sum(ws.*pred)/sw;
    sensorMAE(iSensor)=sum(ws.*abs(err))/sw;
    sensorRMSE(iSensor)=sqrt(sum(ws.*err.^2)/sw);
    rmsObs=sqrt(sum(ws.*obs.^2)/sw);
    sensorNRMSE(iSensor)=sensorRMSE(iSensor)/max(rmsObs,eps);
    sensorBias(iSensor)=sum(ws.*err)/sw;
    sensorBiasPercent(iSensor)=100*sensorBias(iSensor)/max(abs(sensorMeasuredMean(iSensor)),eps);
    sensorMaxAbsError(iSensor)=max(abs(err));
end

% Aggregate true section-holdout validation by held-out sensor. Weights are the
% number of actually observed native intervals in each successfully fitted slice.
sectionValMeasuredMean=nan(m,1); sectionValPredictedMean=nan(m,1);
sectionValMAE=nan(m,1); sectionValRMSE=nan(m,1); sectionValNRMSE=nan(m,1);
sectionValBias=nan(m,1); sectionValBiasPercent=nan(m,1); sectionValMaxAbsError=nan(m,1);
for iSensor=reshape(validationRows,1,[])
    obs=sectionValidationMeasuredMeanByStratum(iSensor,:).';
    pred=sectionValidationPredictedMeanByStratum(iSensor,:).';
    ws=double(sectionValidationObservationCountByStratum(iSensor,:).');
    ok=validFit(:) & isfinite(obs) & isfinite(pred) & ws>0;
    if ~any(ok), continue; end
    obs=obs(ok); pred=pred(ok); ws=ws(ok); sw=max(sum(ws),eps); err=pred-obs;
    sectionValMeasuredMean(iSensor)=sum(ws.*obs)/sw;
    sectionValPredictedMean(iSensor)=sum(ws.*pred)/sw;
    sectionValMAE(iSensor)=sum(ws.*abs(err))/sw;
    sectionValRMSE(iSensor)=sqrt(sum(ws.*err.^2)/sw);
    sectionValNRMSE(iSensor)=sectionValRMSE(iSensor)/max(sqrt(sum(ws.*obs.^2)/sw),eps);
    sectionValBias(iSensor)=sum(ws.*err)/sw;
    sectionValBiasPercent(iSensor)=100*sectionValBias(iSensor)/max(abs(sectionValMeasuredMean(iSensor)),eps);
    sectionValMaxAbsError(iSensor)=max(abs(err));
end
sensorErrorTable=table((1:m).',selectedSections(:),match.SectionID(selectedSections), ...
    match.SectionName(selectedSections),double(match.LaneCount(selectedSections)), ...
    string(match.LaneMemberSectionIDs(selectedSections)),match.RoadName(selectedSections),baseRole, ...
    sensorMeasuredMean,sensorPredictedMean,sensorMAE,sensorRMSE,sensorNRMSE, ...
    sensorBias,sensorBiasPercent,sensorMaxAbsError, ...
    sectionValMeasuredMean,sectionValPredictedMean,sectionValMAE,sectionValRMSE,sectionValNRMSE, ...
    sectionValBias,sectionValBiasPercent,sectionValMaxAbsError, ...
    'VariableNames',{'MeasuredRow','SectionIndex','SectionID','SectionName','LaneCount','LaneMemberSectionIDs','OSMRoad','SectionRole', ...
    'WeightedMeasuredTestMean','WeightedPredictedMean','MAE','RMSE','NRMSE', ...
    'Bias','BiasPercent','MaxAbsoluteError', ...
    'SectionValidationMeasuredMean','SectionValidationPredictedMean','SectionValidationMAE','SectionValidationRMSE','SectionValidationNRMSE', ...
    'SectionValidationBias','SectionValidationBiasPercent','SectionValidationMaxAbsoluteError'});
writetable(sensorErrorTable,o.SensorErrorCSV);
fprintf('Sensor error summary written to: %s\n',o.SensorErrorCSV);

selectedSectionTable = table((1:m).',selectedSections(:), ...
    match.SectionID(selectedSections),match.SectionName(selectedSections), ...
    double(match.LaneCount(selectedSections)),string(match.LaneMemberSectionIDs(selectedSections)), ...
    match.RoadName(selectedSections),match.OSMWayID(selectedSections), ...
    match.DistanceM(selectedSections),match.NameScore(selectedSections), ...
    match.MatchScore(selectedSections),match.DataCoverage(selectedSections),baseRole, ...
    'VariableNames',{'MeasuredRow','SectionIndex','SectionID','SectionName','LaneCount','LaneMemberSectionIDs','OSMRoad', ...
                     'OSMWayID','DistanceM','NameScore','MatchScore','DataCoverage','SectionRole'});
sectionSplitTable=selectedSectionTable;
sectionSplitTable.UsedForEstimation=sectionEstimationMask(:);
sectionSplitTable.HeldOutForValidation=sectionValidationMask(:);
writetable(sectionSplitTable,o.SectionSplitCSV);
fprintf('Section estimation/validation split written to: %s\n',o.SectionSplitCSV);

%% Summary and plots
fprintf('\n================ STRATIFIED PARMA RESULTS ================\n');
fprintf('Observed measurements after lane merging  : %d\n',m);
fprintf('Selected raw lane sections represented    : %d\n',sum(double(match.LaneCount(selectedSections))));
fprintf('Estimation / validation measurements      : %d / %d\n',mEst,mVal);
fprintf('Common path-signature classes             : %d\n',p);
fprintf('Usable weekday/time/date/year strata      : %d\n',Snum);
fprintf('Merge weekdays / merge years              : %d / %d\n',o.MergeWeekdays,o.MergeYears);
if ~isempty(o.DateRange)
    fprintf('Effective estimation date range            : %s to %s (%d days)\n', ...
        datestr(o.DateRange(1),'yyyy-mm-dd'),datestr(o.DateRange(2),'yyyy-mm-dd'),days(o.DateRange(2)-o.DateRange(1))+1);
end
fprintf('Successfully estimated strata             : %d/%d\n',nnz(validFit),Snum);
fprintf('Median retained estimation sensors/slice  : %.1f/%d\n',median(CompleteSensors,'omitnan'),mEst);
fprintf('Median row retention after missing filter : %.1f%%%%\n',100*median(RowRetentionFraction,'omitnan'));
fprintf('Slices without a feasible reduced model : %d\n',nnz(~identifiableByStratum));
fprintf('Common projection design time [s]          : %.3f\n',commonProjectionDesignSeconds);
fprintf('Total stratum estimator time [s]           : %.3f\n',sum(estimationTime(validFit),'omitnan'));
fprintf('Outlier sensor counts removed             : %d\n',nativeInfo.outlierRemoval.valuesRemoved);
if isfield(nativeInfo,'faultCleaning')
    fc=nativeInfo.faultCleaning;
    bad=fc.offDays+fc.stuckDays;
    [~,ordF]=sort(bad,'descend');
    fprintf('Sensors with most off/stuck days set missing:\n');
    for k=ordF(1:min(10,nnz(bad>0))).'
        fprintf('  %4d %-40s off %3d days, stuck %3d days, zero-run intervals %d\n',k, ...
            char(string(match.SectionName(selectedSections(k)))),fc.offDays(k),fc.stuckDays(k),fc.zeroRunIntervals(k));
    end
end
fprintf('Native observation rows affected          : %d\n',nativeInfo.outlierRemoval.observationRowsAffected);
fprintf('Median within-stratum Fano factor         : %.4g\n',medianFano);
fprintf('Median test mean-flow relative error       : %.4g\n',median(testMeanRelErr(validFit),'omitnan'));
fprintf('Median test MAE [vehicles/regime]          : %.4g\n',median(testMAE(validFit),'omitnan'));
fprintf('Median test RMSE [vehicles/regime]         : %.4g\n',median(testRMSE(validFit),'omitnan'));
fprintf('Median test NRMSE                          : %.4g\n',median(testNRMSE(validFit),'omitnan'));
fprintf('Median absolute temporal-test bias [%%]    : %.4g\n',median(abs(testBiasPercent(validFit)),'omitnan'));
if mVal>0
    vv=validFit(:) & isfinite(sectionValidationNRMSE);
    fprintf('Median SECTION-validation MAE              : %.4g\n',median(sectionValidationMAE(vv),'omitnan'));
    fprintf('Median SECTION-validation RMSE             : %.4g\n',median(sectionValidationRMSE(vv),'omitnan'));
    fprintf('Median SECTION-validation NRMSE            : %.4g\n',median(sectionValidationNRMSE(vv),'omitnan'));
    fprintf('Median absolute SECTION-validation bias [%%]: %.4g\n',median(abs(sectionValidationBiasPercent(vv)),'omitnan'));
end
fprintf('==========================================================\n');

if o.Plot
    plotEstimatedClasses(G,repEdgePaths,meanHourlyMu,measuredEdgeSets,selectedSections,string(match.SectionName));

    figure('Name','Parma day-type/time/date Poisson regimes','Color','w');
    tiledlayout('flow');
    groupKey=YearLabel + " | " + DayType;
    groupList=unique(groupKey,'stable');
    for iw=1:numel(groupList)
        nexttile; hold on;
        ix=find(groupKey==groupList(iw));
        totalRate=sum(muRatePerHourByStratum(:,ix),1,'omitnan');
        for j=1:numel(ix)
            q=ix(j);
            rectangle('Position',[StartMinute(q)/60,0,DurationMinutes(q)/60,max(totalRate(j),eps)], ...
                'FaceColor',[0.85 0.85 0.85],'EdgeColor',[0.3 0.3 0.3]);
            text((StartMinute(q)+EndMinute(q))/120,max(totalRate(j),eps),sprintf('S%d',q), ...
                'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',7);
        end
        xlabel('Hour of day'); ylabel('Estimated total latent vehicles/hour');
        title(sprintf('%s (seasonal strata overlaid)',char(groupList(iw)))); grid on; xlim(o.TimeWindowHours);
    end

    % True SECTION hold-out measured versus predicted means.
    obsPlot=SectionValidationMeasuredMean(isfinite(SectionValidationMeasuredMean) & isfinite(SectionValidationPredictedMean));
    predPlot=SectionValidationPredictedMean(isfinite(SectionValidationMeasuredMean) & isfinite(SectionValidationPredictedMean));
    if ~isempty(obsPlot)
        figure('Name','Parma section-holdout validation: measured vs predicted','Color','w');
        scatter(obsPlot,predPlot,26,'filled'); hold on;
        lim=max([obsPlot;predPlot]);
        if ~isfinite(lim) || lim<=0, lim=1; end
        plot([0 lim],[0 lim],'k--','LineWidth',1.2);
        xlabel('Measured hold-out mean flow'); ylabel('Predicted mean flow A\mu');
        title('Parma: measured vs predicted flow by sensor and Poisson stratum');
        axis equal; xlim([0 lim]); ylim([0 lim]); grid on; hold off;
    end

    figure('Name','Parma sensor prediction errors','Color','w');
    if any(sectionValidationMask)
        bar(1:m,sectionValNRMSE);
        ylabel('Section-holdout NRMSE');
        title('Parma: validation error on sensors excluded from estimation');
    else
        bar(1:m,sensorNRMSE);
        ylabel('Temporal hold-out NRMSE');
        title('Parma: temporal prediction error by estimation sensor');
    end
    xlabel('Measured sensor index'); grid on;
end

%% Results
results=struct();
results.codeVersion=codeVersion;
results.animationReady=true;
results.settings=o;
results.requestedDateRange=requestedDateRange;
results.effectiveDateRange=o.DateRange;
results.stableInterval=stableInterval;
results.compactFile=compactFile;
results.osmCacheFile=osmCacheFile;
results.bbox=bbox;
results.G=G;
results.officialRoadMatching=officialRoadSummaryLocal(officialRoads);
results.sectionMatch=match;
results.selectedSections=selectedSections;
results.selectedSectionGroups=selectedGroups;
results.laneMergeTable=laneMergeTable;
results.selectedSectionTable=selectedSectionTable;
results.sectionSplitTable=sectionSplitTable;
results.sectionEstimationMask=sectionEstimationMask;
results.sectionValidationMask=sectionValidationMask;
results.estimationSectionRows=estimationRows;
results.validationSectionRows=validationRows;
results.estimationSections=selectedSections(sectionEstimationMask);
results.validationSections=selectedSections(sectionValidationMask);
results.sectionSplitInfo=sectionSplitInfo;
results.selectionDates=datesSelection;
results.Yselection=Yselection;
results.A=A;
results.AEstimation=A(sectionEstimationMask,:);
results.AValidation=A(sectionValidationMask,:);
results.allCandidateSignatures=sig.Aall;
results.selectedSignatureIndices=selectedClass;
results.representativeEdgePaths=repEdgePaths;
results.representativeNodePaths=repNodePaths;
% Every retained measurement-signature class keeps all DISTINCT physical
% shortest paths sampled for that signature.  The animator divides the
% inferred class rate across these member paths instead of duplicating the
% full class rate on each alternative. Sampling multiplicities are saved
% separately so visualization can weight frequently sampled alternatives.
results.classEdgePaths=classEdgePaths;
% Sampling multiplicity of each distinct physical member path. These weights
% allow the animator to represent all considered member paths without
% pretending that rarely and frequently sampled alternatives are identical.
results.classPathHits=classPathHits;
results.pathClassHits=classHits;
results.pathClassPilot=classPilot;
results.strata=strata;
results.strataTable=strataTable;
results.muHatByStratum=muHatByStratum;
results.muRatePerHourByStratum=muRatePerHourByStratum;
results.meanHourlyMu=meanHourlyMu;
results.infoByStratum=infoByStratum;
results.classDefinitionTable=classDefinitionTable;
results.classTable=classTable;
results.flowErrorTable=flowErrorTable;
results.sectionValidationTable=sectionValidationTable;
results.sensorErrorTable=sensorErrorTable;
results.predictedMeanByStratum=predMeanByStratum;
results.predictedMeanNNLSByStratum=predMeanNNLSByStratum;
results.trainMeasuredMeanByStratum=trainMeasuredMeanByStratum;
results.testMeasuredMeanByStratum=testMeasuredMeanByStratum;
results.allMeasuredMeanByStratum=allMeasuredMeanByStratum;
results.sectionValidationMeasuredMeanByStratum=sectionValidationMeasuredMeanByStratum;
results.sectionValidationPredictedMeanByStratum=sectionValidationPredictedMeanByStratum;
results.sectionValidationPredictedMeanNNLSByStratum=sectionValidationPredictedMeanNNLSByStratum;
results.lag1AutocorrByStratum=lag1AutocorrByStratum;
results.medianFanoByStratum=medianFanoByStratum;
results.projectionRankByStratum=projectionRankByStratum;
results.rankConditionByStratum=rankConditionByStratum;
results.reducedClassCountFitted=reducedClassCountFitted;
results.estimationTimeByStratum=estimationTime;
results.nTrainByStratum=nTrainVec(:);
if exist('machineSummary','file')==2, results.machine=machineSummary(); end
results.validFitByStratum=validFit;
results.reliableByStratum=reliableByStratum;
results.reliabilityReason=reliabilityReason;
results.sensitivityConditionByStratum=sensitivityConditionByStratum;
results.medianRelativeSEByStratum=medianRelativeSEByStratum;
results.trainMisfitByStratum=trainMisfitByStratum;
results.sectionAudit=sectionAudit;
results.sensorResiduals=sensorResiduals;
if isfield(nativeInfo,'faultCleaning'), results.faultCleaning=nativeInfo.faultCleaning; end
results.sensorDiagnostics=sensorDiagnostics;
results.sensorSelection=sensorSelection;
results.indistinguishableSensors=indistinguishable;
results.validationOutliers=validationOutliers;
results.testDaysByStratum=testDaysByStratum;
results.nnlsRelResidualByStratum=nnlsRelResidualByStratum;
results.nnlsMisfitByStratum=nnlsMisfitByStratum;
results.exitflagByStratum=exitflagByStratum;
results.downweightedFractionByStratum=downweightedFractionByStratum;
results.multistartAtBestByStratum=multistartAtBestByStratum;
results.multistartDistinctByStratum=multistartDistinctByStratum;
results.multistartMaxRelDistByStratum=multistartMaxRelDistByStratum;
results.multistartGapByStratum=multistartGapByStratum;
results.muHatUnreliableByStratum=muHatUnreliableByStratum;
results.sectionValidationObservationCountByStratum=sectionValidationObservationCountByStratum;
results.pooledScreeningMean=meanMeasured;
results.pooledScreeningMeanEstimation=meanMeasured(sectionEstimationMask);
results.stratumFano=stratumFano;
results.commonProjectionMatrix=commonProjectionMatrix;
results.commonProjectionDesign=commonProjectionDesign;
results.commonProjectionDesignSeconds=commonProjectionDesignSeconds;
results.projectionDesignGroups=projectionDesignGroups;
results.sensorCompleteByStratum=sensorCompleteByStratum;
results.sensorUsedByStratum=sensorUsedByStratum;
results.identifiableByStratum=identifiableByStratum;
results.originalClassesIdentifiableByStratum=originalClassesIdentifiableByStratum;
results.classReductionByStratum=classReductionByStratum;
results.reducedClassCountByStratum=reducedClassCountByStratum;
results.unobservedClassByStratum=unobservedClassByStratum;
results.sliceClassHandling=o.SliceClassHandling;
if strcmp(o.SliceClassHandling,'merge')
    results.classAllocationPolicy="Equal shares within nonzero reduced classes; zero for unobserved classes (not an estimate)";
else
    results.classAllocationPolicy="Strict: original classes; skip strata with unobserved or indistinguishable classes";
end
results.identifiabilityReason=identifiabilityReason;
results.diagnostics=struct('medianWithinStratumFano',medianFano, ...
    'estimatorFastMode',o.EstimatorFastMode, ...
    'reuseProjectionDesign',o.ReuseProjectionDesign, ...
    'sliceMissingDataProtection',o.ExcludeMissingSensorsPerSlice, ...
    'sliceSensorMinCoverage',o.SliceSensorMinCoverage, ...
    'autoStableDateRange',o.AutoStableDateRange, ...
    'stableInterval',stableInterval, ...
    'sectionHoldoutValidation',o.SectionHoldoutValidation, ...
    'sectionEstimationCount',mEst,'sectionValidationCount',mVal, ...
    'sectionValidationTargetCount',sectionSplitInfo.TargetValidationCount, ...
    'medianSectionValidationMAE',median(sectionValidationMAE(validFit(:) & isfinite(sectionValidationMAE(:))),'omitnan'), ...
    'medianSectionValidationRMSE',median(sectionValidationRMSE(validFit(:) & isfinite(sectionValidationRMSE(:))),'omitnan'), ...
    'medianSectionValidationNRMSE',median(sectionValidationNRMSE(validFit(:) & isfinite(sectionValidationNRMSE(:))),'omitnan'), ...
    'projectionDesignMaskGroups',numel(projectionDesignGroups), ...
    'unidentifiableSlicesAfterSensorExclusion',nnz(~identifiableByStratum), ...
    'commonProjectionDesignSeconds',commonProjectionDesignSeconds, ...
    'totalStratumEstimatorSeconds',sum(estimationTime(validFit),'omitnan'), ...
    'outlierValuesRemoved',nativeInfo.outlierRemoval.valuesRemoved, ...
    'outlierObservationRowsAffected',nativeInfo.outlierRemoval.observationRowsAffected, ...
    'strataAboveFanoTarget',nHighFano, ...
    'fanoTarget',o.MaxMedianFanoPerStratum, ...
    'medianTestMAE',median(testMAE(validFit),'omitnan'), ...
    'medianTestRMSE',median(testRMSE(validFit),'omitnan'), ...
    'medianTestNRMSE',median(testNRMSE(validFit),'omitnan'), ...
    'medianAbsoluteTestBiasPercent',median(abs(testBiasPercent(validFit)),'omitnan'), ...
    'pathDictionaryPilotRelativeResidual',pilotFit, ...
    'sourceAggregationIntervalMinutes',intervalMinutes, ...
    'selectionAggregation',dataInfo,'nativeStratification',nativeInfo);

% Dedicated Fano decomposition.  This is diagnostic only: it uses the final
% cleaned native-interval observations already stored in STRATA and never
% changes the tomography estimates or the temporal partition.
if o.FanoDecomposition
    if exist('diagnoseParmaFanoDecomposition','file')==2
        try
            fanoDecomposition=diagnoseParmaFanoDecomposition(results, ...
                'Plot',false, ...
                'MinGroupSamples',o.FanoDecompositionMinGroupSamples, ...
                'MinMeanCount',o.FanoMinMeanCount, ...
                'SummaryCSV',o.FanoDecompositionCSV, ...
                'SensorCSV',o.FanoDecompositionSensorCSV);
            results.fanoDecomposition=fanoDecomposition;
            results.diagnostics.fanoDecompositionGlobal=fanoDecomposition.global;
        catch MEfanoDecomp
            warning('runParmaTrafficTomography:FanoDecompositionFailed', ...
                'Fano decomposition diagnostic failed: %s',MEfanoDecomp.message);
        end
    else
        warning('runParmaTrafficTomography:FanoDecompositionMissing', ...
            ['FanoDecomposition=true but diagnoseParmaFanoDecomposition.m is not on the MATLAB path. ', ...
             'Estimation results are still valid.']);
    end
end

if ~isempty(o.SaveFile)
    save(o.SaveFile,'results','-v7.3');
    fprintf('Results saved to: %s\n',o.SaveFile);
    fprintf('Animation can be launched later without recomputation using:\n');
    fprintf('  animateParmaTrafficFlows(''%s'',''DayType'',''%s'',''SeasonDate'',''auto'',''PoissonSlice'',''auto'');\n', ...
        o.SaveFile,o.AnimationDayType);
end

% Optional moving-car visualization.  This is a visualization of estimated
% latent path flow rates, not a claim that individual vehicles were tracked.
if o.AnimateTraffic
    if exist('animateParmaTrafficFlows','file')==2
        try
            animateParmaTrafficFlows(results, ...
                'DayType',o.AnimationDayType, ...
                'SeasonDate',o.AnimationSeasonDate, ...
                'Year',o.AnimationYear, ...
                'PoissonSlice',o.AnimationPoissonSlice, ...
                'DurationSeconds',o.AnimationDurationSeconds, ...
                'FrameRate',o.AnimationFrameRate, ...
                'MaxCars',o.AnimationMaxCars, ...
                'VehicleSpeedKmH',o.AnimationVehicleSpeedKmH, ...
                'MinEdgeTraversalsPerSecond',o.AnimationMinEdgeTraversalsPerSecond, ...
                'VideoFile',o.AnimationVideoFile, ...
                'Seed',o.Seed+909);
        catch MEanimation
            warning('runParmaTrafficTomography:AnimationFailed', ...
                'Traffic animation could not be displayed: %s',MEanimation.message);
        end
    else
        warning('runParmaTrafficTomography:AnimationMissing', ...
            'AnimateTraffic=true but animateParmaTrafficFlows.m is not on the MATLAB path.');
    end
end
end


% =========================================================================
function r = stratumLag1AutocorrLocal(Ys,stratum,meanCounts,minMean)
% Median, over sensors with mean count >= minMean, of the correlation between
% counts in consecutive native intervals of the same day within the stratum.
r=NaN;
if ~isfield(stratum,'Slots') || isempty(stratum.Slots) || size(Ys,1)<3, return; end
d=stratum.Dates(:); sl=double(stratum.Slots(:));
if numel(d)~=size(Ys,1), return; end
pair=find(d(1:end-1)==d(2:end) & sl(2:end)-sl(1:end-1)==1);
if numel(pair)<3, return; end
keep=find(meanCounts(:)>=minMean);
rs=nan(numel(keep),1);
for t=1:numel(keep)
    a=Ys(pair,keep(t)); b=Ys(pair+1,keep(t));
    if std(a)>0 && std(b)>0
        cc=corrcoef(a,b);
        rs(t)=cc(1,2);
    end
end
r=median(rs,'omitnan');
end

function tf = isScalarLogicalLike(x)
tf = (islogical(x) || isnumeric(x)) && isscalar(x) && any(x==[0 1]);
end

function tf = isPositiveInteger(x)
tf = isnumeric(x) && isscalar(x) && isfinite(x) && x>=1 && x==floor(x);
end


% =========================================================================
function bbox = sensorBoundingBox(lat,lon,marginKm)
lat=lat(:); lon=lon(:);
lat0=median(lat);
kmPerDegLat=111.32;
kmPerDegLon=111.32*max(cosd(lat0),0.1);
dLat=marginKm/kmPerDegLat;
dLon=marginKm/kmPerDegLon;
bbox=[min(lat)-dLat,min(lon)-dLon,max(lat)+dLat,max(lon)+dLon];
end


% =========================================================================
function [osm,cacheFile] = getOSMRoadData(bbox,o)
% Download the road graph robustly. Large Overpass queries frequently time out
% on public instances, so the sensor bounding box is split into smaller tiles.
% Each successful tile is cached independently and duplicate OSM elements are
% removed after merging. If one endpoint fails, public mirrors are tried.

if ~exist(o.CacheDir,'dir'), mkdir(o.CacheDir); end
classes=strjoin(o.RoadClasses,'|');

% Overall merged cache. Include tiling in the signature because a different
% tile size can change which boundary-crossing ways are returned/cached.
signature=sprintf('%.5f_%.5f_%.5f_%.5f_%s_T%.3f', ...
    bbox(1),bbox(2),bbox(3),bbox(4),classes,o.OverpassTileKm);
h=simpleHash(signature);
cacheFile=fullfile(o.CacheDir,sprintf('parma_roads_%08x.mat',h));
if o.UseCache && exist(cacheFile,'file')==2
    S=load(cacheFile,'osm');
    osm=S.osm;
    fprintf('Loaded merged OSM road data from cache: %s\n',cacheFile);
    return
end

% Primary endpoint first, followed by mirrors. UNIQUE(...,'stable') prevents
% retrying the same endpoint when the user supplied one of the mirrors.
endpoints=[{o.OverpassURL}; o.OverpassMirrors(:)];
endpoints=cellstr(unique(string(endpoints),'stable'));

tiles=splitBBoxKm(bbox,o.OverpassTileKm);
nTiles=size(tiles,1);
fprintf(['Downloading Parma roads from Overpass in %d tiles ', ...
    '(target %.1f km per tile).\n'],nTiles,o.OverpassTileKm);

parts=cell(nTiles,1);
for q=1:nTiles
    tile=tiles(q,:);
    tileSignature=sprintf('%.7f_%.7f_%.7f_%.7f_%s', ...
        tile(1),tile(2),tile(3),tile(4),classes);
    tileHash=simpleHash(tileSignature);
    tileCache=fullfile(o.CacheDir,sprintf('parma_roads_tile_%08x.mat',tileHash));

    if o.UseCache && exist(tileCache,'file')==2
        S=load(tileCache,'part');
        if isfield(S,'part') && isfield(S.part,'elements') && ~isempty(S.part.elements)
            parts{q}=S.part;
            fprintf('  tile %d/%d: cache\n',q,nTiles);
            continue
        end
    end

    fprintf('  tile %d/%d: [%.5f %.5f %.5f %.5f]\n', ...
        q,nTiles,tile(1),tile(2),tile(3),tile(4));
    [part,query,endpointUsed]=downloadOverpassTile(tile,classes,endpoints,o);
    parts{q}=part;
    save(tileCache,'part','tile','query','endpointUsed','-v7.3');
end

osm=mergeOSMResponses(parts);
if ~isfield(osm,'elements') || isempty(osm.elements)
    error('runParmaTrafficTomography:OverpassEmpty', ...
        'The tiled Overpass download returned no road data.');
end

save(cacheFile,'osm','bbox','tiles','endpoints','-v7.3');
fprintf('Saved merged OSM cache: %s\n',cacheFile);
end


% =========================================================================
function tiles = splitBBoxKm(bbox,targetKm)
% Split [south west north east] into approximately targetKm x targetKm tiles.
south=bbox(1); west=bbox(2); north=bbox(3); east=bbox(4);
lat0=(south+north)/2;
heightKm=max((north-south)*111.32,0);
widthKm=max((east-west)*111.32*max(cosd(lat0),0.1),0);
nLat=max(1,ceil(heightKm/targetKm));
nLon=max(1,ceil(widthKm/targetKm));
latBreaks=linspace(south,north,nLat+1);
lonBreaks=linspace(west,east,nLon+1);
tiles=zeros(nLat*nLon,4);
q=0;
for i=1:nLat
    for j=1:nLon
        q=q+1;
        tiles(q,:)=[latBreaks(i),lonBreaks(j),latBreaks(i+1),lonBreaks(j+1)];
    end
end
end


% =========================================================================
function [part,query,endpointUsed] = downloadOverpassTile(tile,classes,endpoints,o)
south=tile(1); west=tile(2); north=tile(3); east=tile(4);
query=sprintf([ ...
    '[out:json][timeout:%d];' ...
    '(way["highway"~"^(%s)$"](%.7f,%.7f,%.7f,%.7f););' ...
    'out body;>;out skel qt;'], ...
    round(o.OverpassTimeout),classes,south,west,north,east);
opts=weboptions('Timeout',o.OverpassTimeout+30,'ContentType','json', ...
    'UserAgent',o.UserAgent);

lastMessage="";
for ie=1:numel(endpoints)
    endpoint=endpoints{ie};
    for attempt=1:o.OverpassRetries
        try
            part=webread(endpoint,'data',query,opts);
            if isfield(part,'elements') && ~isempty(part.elements)
                endpointUsed=endpoint;
                fprintf('    OK: %s\n',endpoint);
                return
            end
            lastMessage="empty response";
        catch ME
            lastMessage=string(ME.message);
        end

        fprintf('    attempt %d/%d failed on %s: %s\n', ...
            attempt,o.OverpassRetries,endpoint,char(lastMessage));
        if o.OverpassRetryPause>0
            pause(o.OverpassRetryPause*attempt);
        end
    end
end

error('runParmaTrafficTomography:OverpassTileFailed', ...
    ['All Overpass endpoints failed for tile [%.6f %.6f %.6f %.6f].\n' ...
     'Last error: %s\n' ...
     'Successful earlier tiles remain cached, so rerunning will resume.'], ...
    south,west,north,east,char(lastMessage));
end


% =========================================================================
function osm = mergeOSMResponses(parts)
% Merge node/way elements by (type,id). A way intersecting adjacent tiles can
% appear multiple times; retaining one copy is sufficient because Overpass
% returns the complete way body plus referenced nodes for each tile query.
seen=containers.Map('KeyType','char','ValueType','logical');
merged=cell(0,1);
for q=1:numel(parts)
    P=parts{q};
    if isempty(P) || ~isfield(P,'elements'), continue; end
    E=P.elements;
    for k=1:elementCount(E)
        e=getElement(E,k);
        if ~isstruct(e) || ~isfield(e,'type') || ~isfield(e,'id')
            continue
        end
        key=sprintf('%s_%s',char(string(e.type)),idKey(e.id));
        if ~isKey(seen,key)
            seen(key)=true;
            merged{end+1,1}=e; %#ok<AGROW>
        end
    end
end
osm=struct();
osm.elements=merged;
end


% =========================================================================
function G = buildDirectedRoadGraph(osm)
E=osm.elements;
numElements=elementCount(E);
nodeIDs=zeros(0,1); nodeLat=zeros(0,1); nodeLon=zeros(0,1);
for k=1:numElements
    e=getElement(E,k);
    if isfield(e,'type') && strcmpi(string(e.type),'node') && ...
            isfield(e,'id') && isfield(e,'lat') && isfield(e,'lon')
        nodeIDs(end+1,1)=double(e.id); %#ok<AGROW>
        nodeLat(end+1,1)=double(e.lat); %#ok<AGROW>
        nodeLon(end+1,1)=double(e.lon); %#ok<AGROW>
    end
end
if isempty(nodeIDs)
    error('runParmaTrafficTomography:NoOSMNodes','No OSM nodes found.');
end
[nodeIDs,ia]=unique(nodeIDs,'stable');
nodeLat=nodeLat(ia); nodeLon=nodeLon(ia);
nNodes=numel(nodeIDs);
nodeMap=containers.Map('KeyType','char','ValueType','double');
for k=1:nNodes, nodeMap(idKey(nodeIDs(k)))=k; end

s=zeros(0,1); t=zeros(0,1); w=zeros(0,1);
highway=strings(0,1); roadName=strings(0,1); wayID=zeros(0,1);
for k=1:numElements
    e=getElement(E,k);
    if ~isfield(e,'type') || ~strcmpi(string(e.type),'way') || ~isfield(e,'nodes'), continue; end
    tags=struct(); if isfield(e,'tags'), tags=e.tags; end
    hw=getTag(tags,'highway','');
    if strlength(hw)==0, continue; end
    access=lower(getTag(tags,'access',''));
    if any(access==["no","private"]), continue; end
    name=getTag(tags,'name','');
    ids=double(e.nodes(:));
    if numel(ids)<2, continue; end
    mode=roadDirection(tags,hw);
    for j=1:numel(ids)-1
        k1=idKey(ids(j)); k2=idKey(ids(j+1));
        if ~isKey(nodeMap,k1) || ~isKey(nodeMap,k2), continue; end
        u=nodeMap(k1); v=nodeMap(k2);
        d=haversineKm(nodeLat(u),nodeLon(u),nodeLat(v),nodeLon(v));
        if ~isfinite(d) || d<=0, continue; end
        if mode>=0
            s(end+1,1)=u; t(end+1,1)=v; w(end+1,1)=d; %#ok<AGROW>
            highway(end+1,1)=hw; roadName(end+1,1)=name; wayID(end+1,1)=double(e.id); %#ok<AGROW>
        end
        if mode<=0
            s(end+1,1)=v; t(end+1,1)=u; w(end+1,1)=d; %#ok<AGROW>
            highway(end+1,1)=hw; roadName(end+1,1)=name; wayID(end+1,1)=double(e.id); %#ok<AGROW>
        end
    end
end
if isempty(s)
    error('runParmaTrafficTomography:NoOSMEdges','No directed road edges were constructed.');
end

T=table(s,t,w,highway,roadName,wayID, ...
    'VariableNames',{'S','T','Weight','Highway','RoadName','OSMWayID'});
T=sortrows(T,{'S','T','Weight'});
[~,keep]=unique([T.S T.T],'rows','stable');
T=T(keep,:);
G=digraph(T.S,T.T,T.Weight,nNodes);
edgeIDs=findedge(G,T.S,T.T);
H=strings(numedges(G),1); R=strings(numedges(G),1); WID=nan(numedges(G),1);
H(edgeIDs)=T.Highway; R(edgeIDs)=T.RoadName; WID(edgeIDs)=T.OSMWayID;
G.Edges.Highway=H; G.Edges.RoadName=R; G.Edges.OSMWayID=WID;
G.Nodes.OSMID=nodeIDs; G.Nodes.Lat=nodeLat; G.Nodes.Lon=nodeLon;

% Keep largest weak component.
bin=conncomp(G,'Type','weak');
counts=accumarray(bin(:),1);
[~,largest]=max(counts);
G=subgraph(G,find(bin==largest));
iso=indegree(G)+outdegree(G)==0;
if any(iso), G=rmnode(G,find(iso)); end
end


% =========================================================================
function mode = roadDirection(tags,highway)
oneway=lower(getTag(tags,'oneway',''));
if any(oneway==["no","0","false"])
    mode=0; return
elseif any(oneway==["yes","1","true"])
    mode=1; return
elseif oneway=="-1"
    mode=-1; return
end
junction=lower(getTag(tags,'junction',''));
if any(junction==["roundabout","circular"]) || lower(highway)=="motorway"
    mode=1;
else
    mode=0;
end
end


% =========================================================================
function match = matchSectionsToRoads(G,sectionId,sectionName,lat,lon,validGPS,cityCenter,o,officialRoads)
% Match traffic sections to OSM physical segments, optionally assisted by the
% official Comune di Parma street-centerline layer and documented detector
% direction rules.  The official layer identifies the road/arm; OSM remains
% the directed graph used by tomography.
if nargin<9 || isempty(officialRoads)
    officialRoads=emptyOfficialRoadLayerLocal();
end
[sEdge,tEdge]=edgeEndpointIndices(G,(1:numedges(G)).');
pair=sort([sEdge tEdge],2);
[groupPair,~,grp]=unique(pair,'rows','stable');
nGroups=size(groupPair,1);
groupEdges=accumarray(grp,(1:numedges(G)).',[],@(x){x});

% Representative group metadata.
groupName=strings(nGroups,1); groupWay=nan(nGroups,1);
for g=1:nGroups
    ids=groupEdges{g};
    names=G.Edges.RoadName(ids);
    q=find(strlength(strtrim(names))>0,1);
    if ~isempty(q), groupName(g)=names(q); else, groupName(g)=""; end
    wid=G.Edges.OSMWayID(ids);
    q=find(isfinite(wid),1); if ~isempty(q), groupWay(g)=wid(q); end
end
u=groupPair(:,1); v=groupPair(:,2);
lat1=double(G.Nodes.Lat(u)); lon1=double(G.Nodes.Lon(u));
lat2=double(G.Nodes.Lat(v)); lon2=double(G.Nodes.Lon(v));
groupBearing=bearingDegLocal(lat1,lon1,lat2,lon2);

n=numel(sectionId);
RoadName=strings(n,1); OSMWayID=nan(n,1); DistanceM=nan(n,1);
NameScore=nan(n,1); MatchScore=nan(n,1); SecondScore=nan(n,1);
ScoreMargin=nan(n,1); DirectionRule=strings(n,1); SegmentKey=strings(n,1);
DirectedEdgeIDs=cell(n,1); Accepted=false(n,1);
OfficialRoadName=strings(n,1); OfficialRoadDistanceM=nan(n,1);
OfficialRoadNameScore=nan(n,1); OfficialOSMNameScore=nan(n,1);
OfficialRoadUsed=false(n,1); DirectionEvidence=strings(n,1);
AssignmentSource=strings(n,1); MatchConfidence=strings(n,1); Ambiguous=false(n,1);

for i=1:n
    if ~validGPS(i), continue; end

    off=matchOfficialRoadForSensorLocal(officialRoads,sectionName(i),lat(i),lon(i),o);
    OfficialRoadName(i)=off.name;
    OfficialRoadDistanceM(i)=off.distanceM;
    OfficialRoadNameScore(i)=off.sensorNameScore;
    OfficialRoadUsed(i)=off.used;

    d=pointToSegmentsMeters(lat(i),lon(i),lat1,lon1,lat2,lon2);
    near=find(d<=o.MaxSensorDistanceM);
    if isempty(near)
        [~,ord]=sort(d,'ascend'); near=ord(1:min(12,numel(ord)));
    else
        [~,ord]=sort(d(near),'ascend'); near=near(ord(1:min(50,numel(ord))));
    end
    scores=zeros(numel(near),1); ns=zeros(numel(near),1);
    ons=zeros(numel(near),1); bs=zeros(numel(near),1);
    for q=1:numel(near)
        g=near(q);
        ns(q)=roadNameSimilarity(sectionName(i),groupName(g));
        ds=exp(-0.5*(d(g)/o.DistanceScaleM)^2);
        if off.used
            ons(q)=roadNameSimilarity(off.name,groupName(g));
            if isfinite(off.bearingDeg)
                bs(q)=undirectedBearingCompatibilityLocal(groupBearing(g),off.bearingDeg);
            else
                bs(q)=0.5;
            end
            % Official geometry/name is deliberately given substantial weight:
            % it is the authoritative municipal street axis, while OSM supplies
            % the directed network topology.
            if strlength(strtrim(groupName(g)))==0
                scores(q)=0.52*ds+0.18*ns(q)+0.10*bs(q)+0.20*off.sensorNameScore;
            else
                scores(q)=0.38*ds+0.20*ns(q)+0.32*ons(q)+0.10*bs(q);
            end
        else
            if strlength(strtrim(groupName(g)))==0
                scores(q)=0.75*ds;
            else
                scores(q)=0.65*ds+0.35*ns(q);
            end
        end
    end
    [scores,ord]=sort(scores,'descend'); near=near(ord); ns=ns(ord); ons=ons(ord);
    best=near(1);
    RoadName(i)=groupName(best);
    OSMWayID(i)=groupWay(best);
    DistanceM(i)=d(best);
    NameScore(i)=ns(1);
    MatchScore(i)=scores(1);
    OfficialOSMNameScore(i)=ons(1);
    if numel(scores)>=2, SecondScore(i)=scores(2); else, SecondScore(i)=0; end
    ScoreMargin(i)=MatchScore(i)-SecondScore(i);
    Ambiguous(i)=numel(scores)>=2 && ScoreMargin(i)<o.MinSensorMatchMargin;
    SegmentKey(i)=sprintf('%d-%d',groupPair(best,1),groupPair(best,2));

    hint=officialDirectionHintLocal(sectionId(i),sectionName(i));
    [chosen,rule,evidence]=chooseDirectedEdges(G,groupEdges{best},sectionName(i),cityCenter,hint);
    DirectedEdgeIDs{i}=chosen(:).';
    DirectionRule(i)=rule;
    DirectionEvidence(i)=evidence;

    baseAccept=DistanceM(i)<=o.MaxSensorDistanceM && ...
        MatchScore(i)>=o.MinSensorMatchScore && ~isempty(chosen);
    Accepted(i)=baseAccept && (~o.RejectAmbiguousMatches || ~Ambiguous(i));

    if off.used && hint.used
        AssignmentSource(i)="Comune centerline + municipal detector direction + OSM";
    elseif off.used
        AssignmentSource(i)="Comune centerline + OSM";
    elseif hint.used
        AssignmentSource(i)="Municipal detector direction + OSM GPS/name";
    else
        AssignmentSource(i)="OSM GPS/name fallback";
    end
    if Accepted(i) && MatchScore(i)>=0.70 && ScoreMargin(i)>=max(0.06,o.MinSensorMatchMargin) && ...
            ((off.used && OfficialOSMNameScore(i)>=0.5) || hint.used)
        MatchConfidence(i)="HIGH";
    elseif Accepted(i) && ScoreMargin(i)>=o.MinSensorMatchMargin
        MatchConfidence(i)="MEDIUM";
    else
        MatchConfidence(i)="LOW";
    end
end

match=table(sectionId,sectionName,lat(:),lon(:),RoadName,OSMWayID,DistanceM,NameScore, ...
    MatchScore,SecondScore,ScoreMargin,OfficialRoadName,OfficialRoadDistanceM, ...
    OfficialRoadNameScore,OfficialOSMNameScore,OfficialRoadUsed,DirectionRule, ...
    DirectionEvidence,AssignmentSource,MatchConfidence,Ambiguous,SegmentKey,DirectedEdgeIDs,Accepted, ...
    'VariableNames',{'SectionID','SectionName','Lat','Lon','RoadName','OSMWayID', ...
    'DistanceM','NameScore','MatchScore','SecondScore','ScoreMargin','OfficialRoadName', ...
    'OfficialRoadDistanceM','OfficialRoadNameScore','OfficialOSMNameScore','OfficialRoadUsed', ...
    'DirectionRule','DirectionEvidence','AssignmentSource','MatchConfidence','Ambiguous', ...
    'SegmentKey','DirectedEdgeIDs','Accepted'});
match.DataCoverage=nan(n,1);
end


% =========================================================================
function [chosen,rule,evidence] = chooseDirectedEdges(G,edgeIDs,sensorName,center,hint)
edgeIDs=edgeIDs(:);
if nargin<5 || isempty(hint), hint=emptyDirectionHintLocal(); end
if isempty(edgeIDs)
    chosen=edgeIDs; rule="unmatched"; evidence="no candidate directed OSM edge"; return
end
if numel(edgeIDs)<=1
    chosen=edgeIDs;
    if hint.used
        rule=hint.rule; evidence=hint.evidence+"; OSM segment is one-way";
    else
        rule="oneway"; evidence="OSM one-way topology";
    end
    return
end
name=lower(normalizeText(sensorName));
[s,t]=edgeEndpointIndices(G,edgeIDs);
latS=double(G.Nodes.Lat(s)); lonS=double(G.Nodes.Lon(s));
latT=double(G.Nodes.Lat(t)); lonT=double(G.Nodes.Lon(t));

% Municipal engineering evidence has priority when a documented flow token is
% present in the section label.
if hint.used
    [chosen,ok]=chooseEdgesByRuleLocal(edgeIDs,latS,lonS,latT,lonT,center,hint.rule);
    if ok
        rule=hint.rule; evidence=hint.evidence; return
    end
end

% Direction toward/away from Parma center from API section-name cues.
if contains(name,"ingresso") || contains(name,"verso centro") || ...
        contains(name,"dir centro") || contains(name,"direzione centro")
    [chosen,~]=chooseEdgesByRuleLocal(edgeIDs,latS,lonS,latT,lonT,center,"toward_center");
    rule="toward_center"; evidence="traffic-section name cue"; return
elseif contains(name,"uscita")
    [chosen,~]=chooseEdgesByRuleLocal(edgeIDs,latS,lonS,latT,lonT,center,"away_from_center");
    rule="away_from_center"; evidence="traffic-section name cue"; return
end

% Cardinal words, only when present as separate tokens in the section name.
words=split(name); words=words(words~="");
rule="both_directions";
if any(words=="nord"), rule="north";
elseif any(words=="sud"), rule="south";
elseif any(words=="ovest"), rule="west";
elseif any(words=="est"), rule="east";
end
if rule~="both_directions"
    [chosen,~]=chooseEdgesByRuleLocal(edgeIDs,latS,lonS,latT,lonT,center,rule);
    evidence="traffic-section cardinal cue";
else
    % Without a direction cue, retain both directions of the same physical
    % segment. This treats the section as a cross-section count.
    chosen=edgeIDs; evidence="no directional evidence; both OSM arcs retained";
end
end


function [chosen,ok]=chooseEdgesByRuleLocal(edgeIDs,latS,lonS,latT,lonT,center,rule)
rule=string(rule); chosen=edgeIDs; ok=true;
switch rule
    case "toward_center"
        ds=haversineKm(latS,lonS,center(1),center(2));
        dt=haversineKm(latT,lonT,center(1),center(2));
        [~,q]=max(ds-dt); chosen=edgeIDs(q);
    case "away_from_center"
        ds=haversineKm(latS,lonS,center(1),center(2));
        dt=haversineKm(latT,lonT,center(1),center(2));
        [~,q]=max(dt-ds); chosen=edgeIDs(q);
    case {"north","south","east","west"}
        dy=latT-latS; dx=(lonT-lonS).*cosd((latT+latS)/2);
        nrm=max(hypot(dx,dy),eps);
        switch rule
            case "north", target=[1 0];
            case "south", target=[-1 0];
            case "east", target=[0 1];
            otherwise, target=[0 -1];
        end
        score=(dy./nrm)*target(1)+(dx./nrm)*target(2);
        [~,q]=max(score); chosen=edgeIDs(q);
    otherwise
        ok=false;
end
end


function hint=officialDirectionHintLocal(sectionID,sectionName)
% Documented flow labels from Comune di Parma traffic-detector engineering
% documents.  The rule is applied only when both the site/road and its flow
% token are present in the API section label; this avoids guessing directions.
hint=emptyDirectionHintLocal();
n=lower(normalizeText(sectionName));
w=split(n); w=w(w~="");
containsToken=@(tok) any(w==lower(string(tok)));

if contains(n,"emilia") && containsToken("pr")
    hint=makeDirectionHintLocal("toward_center", ...
        "Comune detector plan: Emilia Ovest-Pini, flow Emilia->PR");
elseif contains(n,"emilia") && containsToken("pc")
    hint=makeDirectionHintLocal("away_from_center", ...
        "Comune detector plan: Emilia Ovest-Pini, flow Emilia->PC");
elseif contains(n,"cremonese") && containsToken("pr")
    hint=makeDirectionHintLocal("toward_center", ...
        "Comune detector plan: Cremonese-Battibue, flow Cremonese->PR");
elseif contains(n,"cremonese") && containsToken("cr")
    hint=makeDirectionHintLocal("away_from_center", ...
        "Comune detector plan: Cremonese-Battibue, flow Cremonese->CR");
elseif contains(n,"gramsci") && containsToken("pr")
    hint=makeDirectionHintLocal("toward_center", ...
        "Comune detector plan: Gramsci-Osacca, flow Gramsci->PR");
elseif contains(n,"gramsci") && containsToken("pc")
    hint=makeDirectionHintLocal("away_from_center", ...
        "Comune detector plan: Gramsci-Osacca, flow Gramsci->PC");
elseif contains(n,"ponte italia") && (containsToken("ovest") || containsToken("west"))
    hint=makeDirectionHintLocal("west", ...
        "Comune detector plan: P.le Marsala, flow Ponte Italia->Ovest");
elseif contains(n,"pastrengo") && containsToken("nord")
    hint=makeDirectionHintLocal("north", ...
        "Comune detector plan: Montebello-Bizzozzero, flow Pastrengo->Nord");
elseif contains(n,"pastrengo") && containsToken("sud")
    hint=makeDirectionHintLocal("south", ...
        "Comune detector plan: Montebello-Bizzozzero, flow Pastrengo->Sud");
elseif contains(n,"montebello") && containsToken("ovest")
    hint=makeDirectionHintLocal("west", ...
        "Comune detector plan: Montebello-Parco del Dono, flow Montebello->Ovest");
end
% Keep SECTIONID in the signature for future exact-ID extensions and to make
% this helper deliberately compatible with numeric or textual section IDs.
if nargin<1 || isempty(sectionID), return; end %#ok<NASGU>
end

function h=emptyDirectionHintLocal()
h=struct('used',false,'rule',"",'evidence',"");
end

function h=makeDirectionHintLocal(rule,evidence)
h=struct('used',true,'rule',string(rule),'evidence',string(evidence));
end


function off=matchOfficialRoadForSensorLocal(Roads,sensorName,lat,lon,o)
off=struct('used',false,'name',"",'distanceM',NaN,'sensorNameScore',NaN,'bearingDeg',NaN);
if ~Roads.available || isempty(Roads.lat1), return; end
d=pointToSegmentsMeters(lat,lon,Roads.lat1,Roads.lon1,Roads.lat2,Roads.lon2);
near=find(d<=o.OfficialRoadSearchRadiusM);
if isempty(near)
    [~,ord]=sort(d,'ascend'); near=ord(1:min(20,numel(ord)));
end
scores=zeros(numel(near),1); ns=zeros(numel(near),1);
for k=1:numel(near)
    q=near(k);
    ns(k)=roadNameSimilarity(sensorName,Roads.name(q));
    ds=exp(-0.5*(d(q)/o.OfficialRoadDistanceScaleM)^2);
    if strlength(strtrim(Roads.name(q)))==0
        scores(k)=0.8*ds;
    else
        scores(k)=0.58*ds+0.42*ns(k);
    end
end
[bestScore,k]=max(scores); q=near(k);
% Require a credible geometry/name combination before allowing the official
% layer to influence the OSM selection.  Very close unnamed axes remain useful.
use=bestScore>=0.38 && d(q)<=max(o.OfficialRoadSearchRadiusM,2*o.OfficialRoadDistanceScaleM);
if use
    off.used=true; off.name=Roads.name(q); off.distanceM=d(q);
    off.sensorNameScore=ns(k); off.bearingDeg=Roads.bearingDeg(q);
end
end


function score=undirectedBearingCompatibilityLocal(a,b)
if ~(isfinite(a)&&isfinite(b)), score=0.5; return; end
d=mod(abs(a-b),180); d=min(d,180-d);
score=max(0,cosd(d))^2;
end

function b=bearingDegLocal(lat1,lon1,lat2,lon2)
dy=double(lat2)-double(lat1);
dx=(double(lon2)-double(lon1)).*cosd((double(lat1)+double(lat2))/2);
b=mod(atan2d(dx,dy),360);
end


function Roads=emptyOfficialRoadLayerLocal()
Roads=struct('available',false,'source',"",'shapefile',"",'nameField',"", ...
    'lat1',zeros(0,1),'lon1',zeros(0,1),'lat2',zeros(0,1),'lon2',zeros(0,1), ...
    'name',strings(0,1),'bearingDeg',zeros(0,1),'message',"");
end

function Roads=loadParmaOfficialRoadLayerLocal(o)
% Try to load the official Comune di Parma 2025 street-centerline shapefile.
% This enhancement is optional: failure leaves the original OSM matcher active.
Roads=emptyOfficialRoadLayerLocal();
if ~o.OfficialDataAssistedMatching, Roads.message="disabled"; return; end
Roads.source="Comune di Parma - Livello cartografico del reticolo stradale di Parma 2025";
shp=string(o.OfficialRoadShapefile);
cacheDir=char(string(o.OfficialDataCacheDir));
if strlength(shp)==0
    if exist(cacheDir,'dir')==7
        d=dir(fullfile(cacheDir,'**','Grafo_stradale_PARMA_2025_RDN2008.shp'));
        if ~isempty(d), shp=string(fullfile(d(1).folder,d(1).name)); end
    end
end
if strlength(shp)==0 && o.OfficialRoadAutoDownload
    try
        if exist(cacheDir,'dir')~=7, mkdir(cacheDir); end
        opts=weboptions('Timeout',45,'UserAgent',o.UserAgent,'ContentType','text');
        html=webread(o.OfficialRoadResourcePage,opts);
        toks=regexp(html,'href\s*=\s*["'']([^"'']+)["'']','tokens','ignorecase');
        links=strings(0,1);
        for k=1:numel(toks)
            z=string(toks{k}{1});
            z=replace(z,'&amp;','&');
            if contains(lower(z),'zip') || contains(lower(z),'download')
                links(end+1,1)=z; %#ok<AGROW>
            end
        end
        hit=find(contains(lower(links),'rdn2008') & contains(lower(links),'zip'),1);
        if isempty(hit), hit=find(contains(lower(links),'zip'),1); end
        if ~isempty(hit)
            url=links(hit);
            if startsWith(url,"//"), url="https:"+url;
            elseif startsWith(url,"/"), url="https://opendata.comune.parma.it"+url; end
            zipFile=fullfile(cacheDir,'Grafo_stradale_PARMA_2025_RDN2008.zip');
            websave(zipFile,char(url),weboptions('Timeout',90,'UserAgent',o.UserAgent));
            unzip(zipFile,cacheDir);
            d=dir(fullfile(cacheDir,'**','Grafo_stradale_PARMA_2025_RDN2008.shp'));
            if ~isempty(d), shp=string(fullfile(d(1).folder,d(1).name)); end
        end
    catch ME
        warning('runParmaTrafficTomography:OfficialRoadDownload', ...
            'Could not auto-download the Comune road layer; using OSM fallback (%s).',ME.message);
    end
end
if strlength(shp)==0 || exist(char(shp),'file')~=2
    Roads.message="official shapefile unavailable; OSM fallback";
    fprintf('Official Parma road layer unavailable; using OSM GPS/name matcher plus municipal direction rules.\n');
    return
end
if ~(exist('shaperead','file')==2 && exist('projcrs','file')==2 && exist('projinv','file')==2)
    Roads.message="Mapping Toolbox shapefile/projection functions unavailable; OSM fallback";
    warning('runParmaTrafficTomography:OfficialRoadToolbox', ...
        ['Official Parma road SHP found but shaperead/projcrs/projinv are unavailable. ', ...
         'Install Mapping Toolbox or supply matching without the official centerline layer.']);
    return
end
try
    S=shaperead(char(shp));
    if isempty(S), error('empty shapefile'); end
    nameField=detectOfficialRoadNameFieldLocal(S,o.OfficialRoadNameField);
    crs=projcrs(7791);
    lat1=zeros(0,1); lon1=zeros(0,1); lat2=zeros(0,1); lon2=zeros(0,1); names=strings(0,1);
    for i=1:numel(S)
        X=double(S(i).X(:)); Y=double(S(i).Y(:));
        good=isfinite(X)&isfinite(Y);
        if nnz(good)<2, continue; end
        % Convert each uninterrupted polyline run to elementary line segments.
        breaks=[0;find(~good);numel(X)+1];
        nm=string(S(i).(nameField));
        for r=1:numel(breaks)-1
            idx=(breaks(r)+1):(breaks(r+1)-1);
            if numel(idx)<2, continue; end
            [la,lo]=projinv(crs,X(idx),Y(idx));
            keep=isfinite(la(1:end-1))&isfinite(lo(1:end-1))&isfinite(la(2:end))&isfinite(lo(2:end));
            if ~any(keep), continue; end
            lat1=[lat1;la(find(keep))]; %#ok<FNDSB,AGROW>
            lon1=[lon1;lo(find(keep))]; %#ok<FNDSB,AGROW>
            k2=find(keep)+1;
            lat2=[lat2;la(k2)]; lon2=[lon2;lo(k2)]; %#ok<AGROW>
            names=[names;repmat(nm,nnz(keep),1)]; %#ok<AGROW>
        end
    end
    Roads.available=~isempty(lat1); Roads.shapefile=shp; Roads.nameField=string(nameField);
    Roads.lat1=lat1; Roads.lon1=lon1; Roads.lat2=lat2; Roads.lon2=lon2; Roads.name=names;
    Roads.bearingDeg=bearingDegLocal(lat1,lon1,lat2,lon2);
    Roads.message=sprintf('%d official centerline segments loaded',numel(lat1));
    fprintf('Official Parma road layer: %d line segments loaded, name field "%s".\n',numel(lat1),nameField);
catch ME
    Roads=emptyOfficialRoadLayerLocal(); Roads.source="Comune di Parma road centerline 2025";
    Roads.message="load failed: "+string(ME.message);
    warning('runParmaTrafficTomography:OfficialRoadLoad', ...
        'Could not read official Parma road layer; using OSM fallback (%s).',ME.message);
end
end

function field=detectOfficialRoadNameFieldLocal(S,userField)
fn=fieldnames(S);
if ~strcmpi(string(userField),"auto")
    hit=find(strcmpi(fn,char(string(userField))),1);
    if isempty(hit), error('Requested OfficialRoadNameField "%s" not found.',string(userField)); end
    field=fn{hit}; return
end
exclude=lower(string({'Geometry','BoundingBox','X','Y'}));
best=-Inf; field='';
for k=1:numel(fn)
    f=fn{k}; if any(lower(string(f))==exclude), continue; end
    key=lower(string(f)); lex=0;
    if contains(key,'denom'), lex=lex+8; end
    if contains(key,'nome'), lex=lex+7; end
    if contains(key,'strad'), lex=lex+6; end
    if contains(key,'topon'), lex=lex+6; end
    if contains(key,'name'), lex=lex+5; end
    vals=strings(min(numel(S),200),1);
    ok=true;
    for i=1:numel(vals)
        try, vals(i)=string(S(i).(f)); catch, ok=false; break; end
    end
    if ~ok, continue; end
    nonempty=strlength(strtrim(vals))>0 & vals~="<missing>";
    if ~any(nonempty), continue; end
    diversity=numel(unique(vals(nonempty)))/nnz(nonempty);
    score=lex+2*mean(nonempty)+diversity;
    if score>best, best=score; field=f; end
end
if isempty(field), error('Could not identify a street-name attribute in the official shapefile.'); end
end

function summary=officialRoadSummaryLocal(Roads)
summary=struct('available',Roads.available,'source',Roads.source,'shapefile',Roads.shapefile, ...
    'nameField',Roads.nameField,'nSegments',numel(Roads.lat1),'message',Roads.message);
end


% =========================================================================
function d = pointToSegmentsMeters(lat0,lon0,lat1,lon1,lat2,lon2)
% Local tangent-plane point-to-segment distance; adequate at city scale.
ky=111320;
kx=111320*cosd(lat0);
x1=(lon1-lon0)*kx; y1=(lat1-lat0)*ky;
x2=(lon2-lon0)*kx; y2=(lat2-lat0)*ky;
vx=x2-x1; vy=y2-y1;
den=vx.^2+vy.^2;
t=-(x1.*vx+y1.*vy)./max(den,eps);
t=max(0,min(1,t));
xc=x1+t.*vx; yc=y1+t.*vy;
d=sqrt(xc.^2+yc.^2);
end


% =========================================================================
function score = roadNameSimilarity(sensorName,roadName)
a=normalizeText(sensorName); b=normalizeText(roadName);
if strlength(a)==0 || strlength(b)==0, score=0; return; end
if contains(a,b) || contains(b,a)
    substringScore=1;
else
    substringScore=0;
end
ta=roadTokens(a); tb=roadTokens(b);
if isempty(ta) || isempty(tb)
    tokenScore=0;
else
    common=numel(intersect(ta,tb));
    tokenScore=common/max(1,min(numel(ta),numel(tb)));
end
score=max(substringScore,tokenScore);
end

function v=ternaryLocal(c,a,b)
if c, v=a; else, v=b; end
end

function c=cleanSectionNamesLocal(names)
% Section names without invisible padding: null, control and non-breaking-space
% characters become spaces, runs of spaces collapse, ends are trimmed.
c=string(names(:));
c=regexprep(c,'[\x00-\x1F\x7F\xA0]',' ');
c=regexprep(c,'\s+',' ');
c=strtrim(c);
end

function [match,laneTable] = mergeCorsiaMatchesLocal(match,sectionId,sectionName,o)
%MERGECORSIAMATCHESLOCAL Collapse colocated lane-specific detector sections.
% A group is formed only when the section name explicitly contains "corsia",
% the lane qualifier can be removed to obtain the same normalized base name,
% and the GPS locations belong to the same local cluster.
n=height(match);
originalAccepted=logical(match.Accepted);
match.LaneCount=ones(n,1);
match.LaneMergedInto=(1:n).';
match.IsLaneRepresentative=true(n,1);
match.LaneGroupKey=strings(n,1);
match.LaneMemberSections=cell(n,1);
match.LaneMemberSectionIDs=strings(n,1);
for i=1:n
    match.LaneMemberSections{i}=i;
    match.LaneMemberSectionIDs(i)=string(sectionId(i));
    match.LaneGroupKey(i)=sprintf('section_%d',i);
end

if o.MergeCorsie
    key=strings(n,1); label=strings(n,1); isLane=false(n,1);
    for i=1:n
        [key(i),label(i),isLane(i)]=laneBaseNameLocal(sectionName(i));
    end
    laneRows=find(isLane & isfinite(double(match.Lat)) & isfinite(double(match.Lon)) & key~="");
    ukey=unique(key(laneRows),'stable');
    for ku=1:numel(ukey)
        pool=laneRows(key(laneRows)==ukey(ku));
        if numel(pool)<2, continue; end
        remaining=pool(:);
        while ~isempty(remaining)
            % Cluster around one physical GPS location. Use distance to the seed,
            % not transitive chaining, so repeated sensors along the same road
            % cannot be accidentally joined through a sequence of 100 m hops.
            seed=remaining(1); remaining(1)=[];
            if isempty(remaining)
                cluster=seed;
            else
                dd=1000*haversineKm(double(match.Lat(remaining)),double(match.Lon(remaining)), ...
                    double(match.Lat(seed)),double(match.Lon(seed)));
                take=dd<=o.LaneMergeMaxDistanceM;
                cluster=[seed;remaining(take)]; %#ok<AGROW>
                remaining=remaining(~take);
            end
            members=sort(unique(cluster(:)));
            if numel(members)<2, continue; end

            acceptedMembers=members(originalAccepted(members));
            if isempty(acceptedMembers), candidates=members; else, candidates=acceptedMembers; end
            sc=double(match.MatchScore(candidates)); sc(~isfinite(sc))=-Inf;
            [~,ib]=max(sc); rep=candidates(ib);
            if isempty(rep) || ~isfinite(rep), rep=members(1); end

            groupKey=sprintf('corsia_%s_%d',char(ukey(ku)),rep);
            memberIDs=strjoin(string(sectionId(members)),'|');
            for mm=reshape(members,1,[])
                match.LaneCount(mm)=numel(members);
                match.LaneMergedInto(mm)=rep;
                match.IsLaneRepresentative(mm)=(mm==rep);
                match.LaneGroupKey(mm)=groupKey;
                match.LaneMemberSections{mm}=members(:).';
                match.LaneMemberSectionIDs(mm)=memberIDs;
                if mm~=rep, match.Accepted(mm)=false; end
            end
            match.Accepted(rep)=any(originalAccepted(members));
            goodGPS=isfinite(double(match.Lat(members))) & isfinite(double(match.Lon(members)));
            if any(goodGPS)
                match.Lat(rep)=mean(double(match.Lat(members(goodGPS))),'omitnan');
                match.Lon(rep)=mean(double(match.Lon(members(goodGPS))),'omitnan');
            end
            if strlength(label(rep))>0, match.SectionName(rep)=label(rep); end
        end
    end
end

reps=find(match.IsLaneRepresentative);
RepresentativeSection=reps(:);
RepresentativeSectionID=string(match.SectionID(reps));
SectionName=string(match.SectionName(reps));
LaneCount=double(match.LaneCount(reps));
MemberSectionIndices=strings(numel(reps),1);
MemberSectionIDs=strings(numel(reps),1);
Accepted=logical(match.Accepted(reps));
Latitude=double(match.Lat(reps)); Longitude=double(match.Lon(reps));
OSMRoad=string(match.RoadName(reps)); OSMWayID=double(match.OSMWayID(reps));
Direction=string(match.DirectionRule(reps));
for k=1:numel(reps)
    rr=reps(k);
    MemberSectionIndices(k)=strjoin(string(match.LaneMemberSections{rr}),',');
    MemberSectionIDs(k)=match.LaneMemberSectionIDs(rr);
end
laneTable=table(RepresentativeSection,RepresentativeSectionID,SectionName,LaneCount, ...
    MemberSectionIndices,MemberSectionIDs,Accepted,Latitude,Longitude,OSMRoad,OSMWayID,Direction);
end

function [key,label,isLane]=laneBaseNameLocal(name)
raw=string(name);
% The section code before '|' (e.g. 'PRM16-M4 | ') differs between lanes of the
% same road, so it must not enter the grouping key.
rawNoCode=regexprep(raw,'^\s*[^|]*\|\s*','');
normRaw=normalizeText(raw);
isLane=contains(normRaw,"corsia");
if ~isLane
    key=normRaw; label=raw; return
end
% Remove the lane word and a conventional lane designator immediately after it.
% Word boundaries are written as lookarounds: in MATLAB regular expressions '\b'
% is a backspace, not a word boundary.
% Covers e.g. "corsia 1", "corsia n 2", "corsia II", "corsia dx/sx".
pat='(?i)(?<![A-Za-z0-9])corsia\s*(?:n(?:r|umero)?\.?\s*)?(?:[0-9]+|[ivxlcdm]+|[a-z]|dx|sx|destra|sinistra|centrale|centro)?(?![A-Za-z0-9])';
patBefore='(?i)(?<![A-Za-z0-9])(?:[0-9]+|[ivxlcdm]+)\s*[\^°ºª]?\s*corsia(?![A-Za-z0-9])';
label=regexprep(raw,pat,' ');
label=regexprep(label,patBefore,' ');
label=regexprep(label,'\s+',' '); label=strtrim(label);
if strlength(label)==0, label=raw; end
keyText=regexprep(rawNoCode,pat,' ');
keyText=regexprep(keyText,patBefore,' ');
key=normalizeText(keyText);
end

function groups=normalizeSectionGroupsLocal(sectionGroups)
if isnumeric(sectionGroups)
    x=double(sectionGroups(:));
    groups=arrayfun(@(v)reshape(v,1,1),x,'UniformOutput',false);
elseif iscell(sectionGroups)
    groups=sectionGroups(:);
    for k=1:numel(groups)
        groups{k}=unique(round(double(groups{k}(:).')),'stable');
        groups{k}=groups{k}(isfinite(groups{k}) & groups{k}>=1);
        if isempty(groups{k}), error('runParmaTrafficTomography:EmptyMeasurementGroup','Measurement group %d is empty.',k); end
    end
else
    error('runParmaTrafficTomography:MeasurementGroups','Section groups must be numeric indices or a cell array of index vectors.');
end
end

function Y=aggregateLaneGroups2DLocal(Yraw,rawSections,groups,requireAll)
nd=size(Yraw,1); ng=numel(groups); Y=nan(nd,ng);
for g=1:ng
    [tf,pos]=ismember(groups{g},rawSections);
    pos=pos(tf & pos>0);
    if isempty(pos), continue; end
    B=Yraw(:,pos);
    if numel(pos)==1
        Y(:,g)=B;
    elseif requireAll
        good=all(isfinite(B),2);
        tmp=sum(B,2,'omitnan'); tmp(~good)=NaN; Y(:,g)=tmp;
    else
        good=any(isfinite(B),2);
        tmp=sum(B,2,'omitnan'); tmp(~good)=NaN; Y(:,g)=tmp;
    end
end
end

function [X,info]=cleanFaultySensorDataLocal(X,intervalMinutes,o,dates,isTestDay)
% Mark as missing (NaN) the data of sensors that are off, stuck or abnormal.
% X is dates x slots x sensors. See the CleanFaultySensorData option.
[D,T,m]=size(X);
info=struct('offDays',zeros(m,1),'stuckDays',zeros(m,1),'zeroRunIntervals',zeros(m,1),'abnormalDays',zeros(m,1));
if nargin<4 || isempty(dates), dates=[]; end
if nargin<5 || isempty(isTestDay), isTestDay=false(D,1); end
learn=~isTestDay(:);                      % thresholds are learned on these days only
stuckLen=max(2,ceil(o.StuckRunMinutes/intervalMinutes));
zeroLen=max(2,ceil(o.ZeroRunMinutes/intervalMinutes));
for j=1:m
    Xj=X(:,:,j);
    fin=isfinite(Xj);
    Xz=Xj; Xz(~fin)=0;
    dayTot=sum(Xz,2);
    full=sum(fin,2)>=T/2;
    if ~any(full), continue; end
    active=full & dayTot>0 & learn;      % typical level from non-test days when the sensor works
    if ~any(active), continue; end
    typical=median(dayTot(active));
    slotTypical=median(Xj(learn,:),1,'omitnan');
    for d=1:D
        if ~full(d), continue; end
        row=Xj(d,:);
        % Off day: almost no traffic recorded on a normally busy sensor
        if typical>=o.OffMinDailyCount && dayTot(d)<=o.OffDayFraction*typical
            Xj(d,:)=NaN; info.offDays(j)=info.offDays(j)+1; continue
        end
        % Stuck day: long run of identical, non-trivial counts
        if longestRunLocal(row,o.StuckMinValue)>=stuckLen
            Xj(d,:)=NaN; info.stuckDays(j)=info.stuckDays(j)+1; continue
        end
        % Partial outage: long run of zeros where the sensor is normally busy
        z=row==0 & slotTypical>=o.ZeroRunMinTypical;
        [st,len]=runsLocal(z);
        for r=find(len>=zeroLen)
            idx=st(r):st(r)+len(r)-1;
            Xj(d,idx)=NaN; info.zeroRunIntervals(j)=info.zeroRunIntervals(j)+numel(idx);
        end
    end
    X(:,:,j)=Xj;
end
% Abnormal days: total far from the sensor's median for the same day type.
% Holidays are not known here and count as their weekday; with a factor of 3
% only gross departures are removed.
if isfield(o,'AbnormalDayFactor') && isfinite(o.AbnormalDayFactor) && ~isempty(dates)
    wd=weekday(dates(:)); dtype=ones(D,1); dtype(wd==7)=2; dtype(wd==1)=3;
    for j=1:m
        Xj=X(:,:,j); fin=isfinite(Xj);
        Xz=Xj; Xz(~fin)=0; tot=sum(Xz,2); full=sum(fin,2)>=T/2;
        for c=1:3
            days=find(dtype==c & full & tot>0);
            ref=days(learn(days));               % median from non-test days
            if numel(ref)<4, continue; end
            med=median(tot(ref));
            if med<o.OffMinDailyCount, continue; end
            bad=days(tot(days)<med/o.AbnormalDayFactor | tot(days)>med*o.AbnormalDayFactor);
            Xj(bad,:)=NaN; info.abnormalDays(j)=info.abnormalDays(j)+numel(bad);
        end
        X(:,:,j)=Xj;
    end
end
end

function r=longestRunLocal(x,minValue)
% Longest run of identical consecutive finite counts that are >= minValue.
r=0; cur=1;
for k=2:numel(x)
    if isfinite(x(k)) && isfinite(x(k-1)) && x(k)==x(k-1) && x(k)>=minValue
        cur=cur+1;
    else
        cur=1;
    end
    r=max(r,cur);
end
end

function [st,len]=runsLocal(z)
% Start indices and lengths of runs of true values in logical vector z.
z=logical(z(:)).'; d=diff([0 z 0]);
st=find(d==1); en=find(d==-1)-1; len=en-st+1;
end

function [B,bad]=markImplausibleDetectorsLocal(B,dimSensor,maxPerUnit)
% Treat as missing every raw detector whose median positive count per unit
% (interval or day) exceeds maxPerUnit: physically impossible traffic.
if dimSensor==2
    m=size(B,2); bad=false(1,m);
    for j=1:m, v=B(:,j); v=v(isfinite(v)&v>0); if ~isempty(v), bad(j)=median(v)>maxPerUnit; end, end
    B(:,bad)=NaN;
else
    m=size(B,3); bad=false(1,m);
    for j=1:m, v=B(:,:,j); v=v(isfinite(v)&v>0); if ~isempty(v), bad(j)=median(v)>maxPerUnit; end, end
    B(:,:,bad)=NaN;
end
end

function [B,dead]=markDeadDetectorsLocal(B,dimSensor,maxZeroShare)
% Treat as missing every raw detector whose finite counts are zero at least
% maxZeroShare of the time (dimSensor: 2 for days x sensors, 3 for days x slots x sensors).
if dimSensor==2
    fin=isfinite(B); nz=sum(fin,1); z=sum(fin & B==0,1);
    dead=nz>0 & z./max(nz,1)>=maxZeroShare;
    B(:,dead)=NaN;
else
    fin=isfinite(B); nz=squeeze(sum(sum(fin,1),2)); z=squeeze(sum(sum(fin & B==0,1),2));
    dead=(nz>0 & z./max(nz,1)>=maxZeroShare).';
    B(:,:,dead)=NaN;
end
end

function X=aggregateLaneGroups3DLocal(Xraw,rawSections,groups,requireAll)
nd=size(Xraw,1); T=size(Xraw,2); ng=numel(groups); X=nan(nd,T,ng);
for g=1:ng
    [tf,pos]=ismember(groups{g},rawSections);
    pos=pos(tf & pos>0);
    if isempty(pos), continue; end
    B=Xraw(:,:,pos);
    if numel(pos)==1
        X(:,:,g)=B(:,:,1);
    elseif requireAll
        good=all(isfinite(B),3);
        tmp=sum(B,3,'omitnan'); tmp(~good)=NaN; X(:,:,g)=tmp;
    else
        good=any(isfinite(B),3);
        tmp=sum(B,3,'omitnan'); tmp(~good)=NaN; X(:,:,g)=tmp;
    end
end
end

function s = normalizeText(x)
s=lower(string(x));
replFrom=["à","á","è","é","ì","í","ò","ó","ù","ú","'","-","/",".",",","(",")"];
replTo  =["a","a","e","e","i","i","o","o","u","u"," "," "," "," "," "," "," "];
for k=1:numel(replFrom), s=replace(s,replFrom(k),replTo(k)); end
s=regexprep(s,'[^a-z0-9 ]',' ');
s=regexprep(s,'\s+',' ');
s=strtrim(s);
end

function t = roadTokens(s)
raw=split(string(s)); raw=raw(raw~="");
stop=["via","viale","v","vl","strada","str","piazza","piazzale","p.le", ...
      "borgo","ponte","rotatoria","rotonda","parma","localita","loc", ...
      "direzione","dir","verso","ingresso","uscita","centro", ...
      "nord","sud","est","ovest","km","svincolo","corsia"];
t=unique(raw(~ismember(raw,stop) & strlength(raw)>=2));
end


% =========================================================================
function [strata,info] = buildPoissonStrata(fileName,sectionGroups,o,modelSensorMask)
sectionGroups=normalizeSectionGroupsLocal(sectionGroups);
rawSections=unique([sectionGroups{:}]);
% Construct exhaustive Poisson strata of the form
%   (time-of-day span) x (weekday/day-type group) x (calendar month/day span)
%   x (year group).
%
% MERGEWEEKDAYS=true:
%   1 Workday        = Monday-Friday excluding holidays
%   2 Saturday       = non-holiday Saturdays
%   3 Sunday/Holiday = Sundays and holidays (holiday status overrides weekday)
%
% MERGEWEEKDAYS=false:
%   Monday, Tuesday, Wednesday, Thursday, Friday, Saturday, Sunday/Holiday.
%   Holidays always enter Sunday/Holiday regardless of their nominal weekday.
%
% MERGEYEARS=true pools all available years on a canonical month/day axis.
% MERGEYEARS=false builds independent time/season strata within each year.
if nargin<4 || isempty(modelSensorMask), modelSensorMask=true(numel(sectionGroups),1); end
modelSensorMask=logical(modelSensorMask(:));
if numel(modelSensorMask)~=numel(sectionGroups)
    error('runParmaTrafficTomography:ModelSensorMask','modelSensorMask must match the number of measurement groups.');
end
S=load(fileName,'typeNames','dateStart','dateEnd','nRecords');
typeIdx=find(string(S.typeNames(:))==o.TrafficType,1);
if isempty(typeIdx)
    error('runParmaTrafficTomography:TrafficType','Traffic type "%s" is unavailable.',o.TrafficType);
end
M=matfile(fileName);
nRecords=double(S.nRecords);
recType=double(M.recType(1:nRecords,1));
recDay=double(M.recDay(1:nRecords,1));
recSection=double(M.recSection(1:nRecords,1));
recInterval=double(M.recInterval(1:nRecords,1));
recStart=double(M.recSignalStart(1:nRecords,1));
recLength=double(M.recSignalLength(1:nRecords,1));
base=recType==typeIdx & ismember(recSection,rawSections) & recLength>0 & ...
    isfinite(recInterval) & recInterval>0;
if ~any(base), error('runParmaTrafficTomography:NoTrafficRecords','No eligible records.'); end

lastAvailable=max(recDay(base)); firstAvailable=min(recDay(base));
if isempty(o.DateRange)
    dayEnd=lastAvailable; dayStart=firstAvailable;
else
    dr=o.DateRange(:);
    dayStart=max(1,round((dr(1)-S.dateStart)/days(1))+1);
    dayEnd=min(round((S.dateEnd-S.dateStart)/days(1))+1,round((dr(2)-S.dateStart)/days(1))+1);
end
allDays=(dayStart:dayEnd).';
allDates=S.dateStart+days(allDays-1);
% Test days are drawn first, so that no learning step uses them: cleaning
% thresholds, segmentation, outlier scores, class screening and projection
% design all use the other days only (common to every stratum and criterion).
% Final test days, common to every stratum and to both stratification criteria
rsTest=RandStream('mt19937ar','Seed',o.Seed+4242);
nTestDays=max(1,round((1-o.TrainFraction)*numel(allDates)));
if strcmpi(o.TestDaySampling,'systematic')
    outerTestDates=systematicTestDates(allDates,o,rsTest);
else
    outerTestDates=allDates(sort(randperm(rsTest,numel(allDates),nTestDays)));
end
isOuterDay=ismember(allDates,outerTestDates);
if isempty(allDays), error('runParmaTrafficTomography:NoEligibleDays','No eligible days.'); end
if o.WeekdaysOnly
    warning('runParmaTrafficTomography:WeekdaysOnlyIgnored', ...
        'WeekdaysOnly is ignored because the configured calendar grouping uses all data.');
end

eligible=base & ismember(recDay,allDays);
intervalMode=mode(recInterval(eligible));
eligible=eligible & abs(recInterval-intervalMode)<1e-6;
maxL=max(recLength(eligible));
allSlotMinutes=(0:maxL-1).'*intervalMode;
windowStart=o.TimeWindowHours(1)*60; windowEnd=o.TimeWindowHours(2)*60;
slotKeep=allSlotMinutes>=windowStart & allSlotMinutes<windowEnd;
slotIndex=find(slotKeep); slotMinutes=allSlotMinutes(slotKeep);
if isempty(slotIndex), error('runParmaTrafficTomography:TimeWindow','No native samples lie in TimeWindowHours.'); end

nd=numel(allDays); m=numel(sectionGroups); T=numel(slotIndex); nRaw=numel(rawSections);
Xraw=nan(nd,T,nRaw);
[~,dayPos]=ismember(recDay,allDays);
[~,secPos]=ismember(recSection,rawSections);
rows=find(eligible & dayPos>0 & secPos>0);
signalChunks=buildSignalReadChunks(rows,recStart,recLength,o.SignalReadChunkMB);
fprintf('Native-strata signal I/O: %d records -> %d chunk reads.\n',numel(rows),numel(signalChunks));
for cRead=1:numel(signalChunks)
    cFirst=signalChunks(cRead).first;
    cLast=signalChunks(cRead).last;
    signalChunk=M.signalData(cFirst:cLast,1);
    rr=signalChunks(cRead).rows;
    for z=1:numel(rr)
        r=rr(z); L=recLength(r);
        src=slotIndex(slotIndex<=L);
        if isempty(src), continue; end
        localFirst=recStart(r)-cFirst+1;
        x=double(signalChunk(localFirst:localFirst+L-1));
        vals=x(src); vals(~isfinite(vals) | vals<0)=NaN;
        d=dayPos(r); j=secPos(r); tgt=1:numel(src);
        old=squeeze(Xraw(d,tgt,j)); old=old(:);
        fill=isnan(old); old(fill)=vals(fill);
        both=~fill & isfinite(vals);
        old(both)=0.5*(old(both)+vals(both));
        Xraw(d,tgt,j)=reshape(old,1,[]);
    end
end

% Collapse raw lane-specific detector arrays into one measurement per physical
% cross-section. Counts are summed across corsie at the native interval level.
[Xraw,deadRaw]=markDeadDetectorsLocal(Xraw,3,o.MaxSensorZeroShare);
[Xraw,implRaw]=markImplausibleDetectorsLocal(Xraw,3,o.MaxFlowPerHour*intervalMode/60);
if any(implRaw), fprintf('Implausible raw detectors (median count > %g veh/h): %d\n',o.MaxFlowPerHour,nnz(implRaw)); end
if any(deadRaw), fprintf('Dead raw detectors (zero counts >= %.0f%% of the time): %d\n',100*o.MaxSensorZeroShare,nnz(deadRaw)); end
X=aggregateLaneGroups3DLocal(Xraw,rawSections,sectionGroups,o.LaneMergeRequireAllLanes);
[X,implSec]=markImplausibleDetectorsLocal(X,3,o.MaxSectionFlowPerHour*intervalMode/60);
if any(implSec)
    fprintf('Implausible measurements after lane summation (median count > %g veh/h): %d (measurement %s)\n', ...
        o.MaxSectionFlowPerHour,nnz(implSec),mat2str(find(implSec)));
end
faultInfo=struct('offDays',zeros(size(X,3),1),'stuckDays',zeros(size(X,3),1), ...
    'zeroRunIntervals',zeros(size(X,3),1),'abnormalDays',zeros(size(X,3),1));
if o.CleanFaultySensorData
    [X,faultInfo]=cleanFaultySensorDataLocal(X,intervalMode,o,allDates,isOuterDay);
    fprintf(['Faulty-sensor cleaning: %d off sensor-days, %d stuck sensor-days, %d abnormal sensor-days, ' ...
        '%d intervals in zero runs set missing.\n'], ...
        sum(faultInfo.offDays),sum(faultInfo.stuckDays),sum(faultInfo.abnormalDays),sum(faultInfo.zeroRunIntervals));
end

% Optional coarser time aggregation (after cleaning, which works at the native interval)
if ~isempty(o.AggregationMinutes) && o.AggregationMinutes>intervalMode+1e-9
    k=o.AggregationMinutes/intervalMode;
    if abs(k-round(k))>1e-9
        error('runParmaTrafficTomography:AggregationMinutes', ...
            'AggregationMinutes=%g must be a multiple of the native interval (%g min).',o.AggregationMinutes,intervalMode);
    end
    k=round(k); T2=floor(size(X,2)/k);
    Xa=nan(size(X,1),T2,size(X,3));
    for t=1:T2
        Xa(:,t,:)=sum(X(:,(t-1)*k+(1:k),:),2);       % NaN if any native interval is missing
    end
    X=Xa; slotMinutes=slotMinutes(1:k:k*T2); T=T2; intervalMode=intervalMode*k;
    fprintf('Counts aggregated to %g-minute intervals (%d per day).\n',intervalMode,T2);
end

[groupCode,holidayMask,groupNames]=classifyTrafficDayType(allDates,o);
calendarDOY=canonicalMonthDayIndex(allDates);
yearAxis=year(allDates);
if o.MergeYears
    yearGroups=0; % 0 denotes pooled AllYears
else
    yearGroups=unique(yearAxis,'stable').';
end
nGroups=numel(groupNames); nYears=numel(yearGroups);

strata=repmat(struct( ...
    'DayTypeCode',0,'DayType','', ...
    'Weekday',0,'WeekdayName','', ...
    'Year',0,'YearLabel','', ...
    'StartMinute',0,'EndMinute',0,'DurationMinutes',0,'ObservationMinutes',intervalMode, ...
    'SeasonStartDOY',0,'SeasonEndDOY',0,'SeasonStart','', 'SeasonEnd','', 'SeasonLabel','', ...
    'Y',[],'YFull',[],'SensorMask',false(1,m),'SensorRows',zeros(0,1),'SensorCoverage',nan(1,m),'Dates',NaT(0,1),'NativeSlots',[],'Slots',zeros(0,1),'BIC',NaN,'TimeBIC',NaN,'SeasonBIC',NaN, ...
    'CandidateDates',0,'CompleteDates',0,'CandidateObservations',0,'CompleteObservations',0,'AllSensorCompleteObservations',0, ...
    'RowsRetainedAfterMissingFilter',0,'RowRetentionFraction',NaN, ...
    'RefinementDepth',0,'RefinedForFano',false),0,1);
timeSegmentation=cell(nYears,nGroups);
seasonSegmentation=cell(nYears,nGroups);
minSlots=max(1,ceil(o.MinPoissonSegmentMinutes/intervalMode));
useCV=strcmpi(o.StrataCriterion,'cv');
cvWorkers=0; if isfield(o,'UseParallel') && o.UseParallel, cvWorkers=Inf; end
isValDay=false(size(allDates)); Acv=[];
if useCV
    % Inner validation days among the non-test days, and a compact class matrix:
    % classes positive in the NNLS fit of the pooled means, plus the most sampled.
    rsVal=RandStream('mt19937ar','Seed',o.Seed+4343);
    isValDay=~isOuterDay(:) & rand(rsVal,numel(allDates),1)<o.CVValidationShare;
    Xm=reshape(X(~isOuterDay,:,modelSensorMask),[],nnz(modelSensorMask));
    ybarAll=zeros(size(Xm,2),1);
    for i=1:size(Xm,2), v=Xm(:,i); v=v(isfinite(v)); if ~isempty(v), ybarAll(i)=mean(v); end, end
    muAll=lsqnonneg(o.CVDictionary,ybarAll,optimset('Display','off'));
    % Classes used by the NNLS fit of the pooled means, plus, for every sensor
    % they leave uncovered, its most sampled crossing class.
    cvCols=muAll(:)>0;
    for i=find(~any(o.CVDictionary(:,cvCols)>0,2)).'
        cand=find(o.CVDictionary(i,:)>0);
        if isempty(cand), continue; end
        [~,jb]=max(o.CVHits(cand)); cvCols(cand(jb))=true;
    end
    Acv=o.CVDictionary(:,cvCols);
    cvTimeBlock=max(1,round(o.CVTimeBlockMinutes/intervalMode));
    fprintf('  Strata criterion: cross-validated NNLS error (%d classes, %d validation days, %d test days)\n', ...
        size(Acv,2),nnz(isValDay),nnz(isOuterDay));
end

for iy=1:nYears
    yv=yearGroups(iy);
    if yv==0, ymask=true(size(allDates)); else, ymask=yearAxis==yv; end
    for dt=1:nGroups
        dmask=groupCode==dt & ymask;
        if ~any(dmask), continue; end
        Xd=X(dmask & ~isOuterDay(:),:,modelSensorMask);   % segmentation learns on non-test days
        C=zeros(T,1); E=zeros(T,1);
        for t=1:T
            aa=squeeze(Xd(:,t,:));
            good=isfinite(aa) & aa>=0;
            C(t)=sum(aa(good)); E(t)=nnz(good);
        end
        segLimit=o.MaxPoissonSegmentsPerWeekday(min(dt,numel(o.MaxPoissonSegmentsPerWeekday)));
        % A window must hold MinStratumIntervals intervals over all days of this type
        minSlotsDT=max(minSlots,ceil(o.MinStratumIntervals/max(nnz(dmask),1)));
        if useCV
            dIdx=find(dmask);
            tr=dIdx(~isOuterDay(dIdx) & ~isValDay(dIdx)); va=dIdx(isValDay(dIdx));
            [Str,Ntr]=slotStatsLocal(X(tr,:,modelSensorMask),false);
            [Sv,Nv,Qv]=slotStatsLocal(X(va,:,modelSensorMask),true);
            [timeBounds,timeBIC]=cvPiecewiseSegments(Str,Ntr,Sv,Qv,Nv,Acv,minSlotsDT,segLimit,o.CVTolerance,cvWorkers,cvTimeBlock);
            fprintf('    CV time windows, day type %d: %d windows\n',dt,size(timeBounds,1));
        else
            [timeBounds,timeBIC]=poissonPiecewiseSegments(C,E,minSlotsDT,segLimit);
        end
        timeSegmentation{iy,dt}=timeBounds;
        seasonSegmentation{iy,dt}=cell(size(timeBounds,1),1);

        for rt=1:size(timeBounds,1)
            a=timeBounds(rt,1); b=timeBounds(rt,2);
            Cseason=zeros(366,1); Eseason=zeros(366,1);
            didx=find(dmask);
            for kk=1:numel(didx)
                d=didx(kk);
                if isOuterDay(d), continue; end               % test days are not used
                v=squeeze(X(d,a:b,modelSensorMask));
                good=isfinite(v) & v>=0;
                if any(good(:))
                    cd=calendarDOY(d);
                    Cseason(cd)=Cseason(cd)+sum(v(good));
                    Eseason(cd)=Eseason(cd)+nnz(good);
                end
            end
            seasonLimit=o.MaxSeasonSegmentsPerTimeInterval(min(dt,numel(o.MaxSeasonSegmentsPerTimeInterval)));
            % A calendar period must contain enough days of this type for the window
            dayCount=accumarray(calendarDOY(didx(:)),1,[366 1]);
            minDays=ceil(o.MinStratumIntervals/(b-a+1));
            if useCV
                mS=nnz(modelSensorMask);
                Str=zeros(366,mS); Ntr=Str; Sv=Str; Nv=Str; Qv=Str;
                for kk=1:numel(didx)
                    d=didx(kk); if isOuterDay(d), continue; end
                    v=reshape(X(d,a:b,modelSensorMask),b-a+1,mS); good=isfinite(v)&v>=0; v(~good)=0;
                    cd=calendarDOY(d);
                    if isValDay(d)
                        Sv(cd,:)=Sv(cd,:)+sum(v,1); Qv(cd,:)=Qv(cd,:)+sum(v.^2,1); Nv(cd,:)=Nv(cd,:)+sum(good,1);
                    else
                        Str(cd,:)=Str(cd,:)+sum(v,1); Ntr(cd,:)=Ntr(cd,:)+sum(good,1);
                    end
                end
                [seasonBounds,seasonBIC]=cvPiecewiseSegments(Str,Ntr,Sv,Qv,Nv,Acv,o.MinSeasonSegmentDays,seasonLimit,o.CVTolerance,cvWorkers,o.CVCalendarBlockDays,dayCount,minDays);
            else
                [seasonBounds,seasonBIC]=poissonPiecewiseSegments( ...
                    Cseason,Eseason,o.MinSeasonSegmentDays,seasonLimit,dayCount,minDays);
            end
            seasonSegmentation{iy,dt}{rt}=seasonBounds;

            for rs=1:size(seasonBounds,1)
                s0=seasonBounds(rs,1); s1=seasonBounds(rs,2);
                sq=makeCalendarStratumLocal(X,allDates,groupCode,calendarDOY,yearAxis,dt,yv,a,b,s0,s1, ...
                    slotMinutes,intervalMode,groupNames,timeBIC,seasonBIC,0,false,o,modelSensorMask);
                if isempty(sq) || isempty(sq.Y), continue; end
                strata(end+1,1)=sq; %#ok<AGROW>
            end
        end
    end
end

outlierMask=false(size(X));
outlierInfo=struct('enabled',o.OutlierRemoval,'valuesRemoved',0, ...
    'observationRowsAffected',0,'preliminaryStrata',numel(strata), ...
    'threshold',o.OutlierMADThreshold,'minSamples',o.OutlierMinSamples, ...
    'maxFractionPerSensor',o.OutlierMaxFractionPerSensor,'modelSensorMask',modelSensorMask);
if o.OutlierRemoval && ~isempty(strata)
    preliminaryStrata=strata;
    [X,outlierMask,outlierInfo]=removeCalendarOutliersLocal( ...
        X,allDates,groupCode,calendarDOY,yearAxis,preliminaryStrata,o,modelSensorMask,isOuterDay);
    rebuilt=repmat(preliminaryStrata(1),0,1);
    for q0=1:numel(preliminaryStrata)
        parent=preliminaryStrata(q0);
        a=parent.NativeSlots(1); b=parent.NativeSlots(end);
        sq=makeCalendarStratumLocal(X,allDates,groupCode,calendarDOY,yearAxis,parent.DayTypeCode,parent.Year, ...
            a,b,parent.SeasonStartDOY,parent.SeasonEndDOY,slotMinutes,intervalMode,groupNames, ...
            parent.TimeBIC,parent.SeasonBIC,parent.RefinementDepth,parent.RefinedForFano,o,modelSensorMask);
        if ~isempty(sq) && ~isempty(sq.Y), rebuilt(end+1,1)=sq; end %#ok<AGROW>
    end
    strata=rebuilt;
end

refinementInfo=struct('enabled',o.FanoRefinement,'splits',0,'passes',0, ...
    'targetMedianFano',o.MaxMedianFanoPerStratum);
if o.FanoRefinement && ~isempty(strata)
    [strata,refinementInfo]=refineCalendarStrataByFanoLocal( ...
        X,allDates,groupCode,calendarDOY,yearAxis,slotMinutes,intervalMode,groupNames,strata,o,modelSensorMask);
end
strata=annotateOutlierStatsLocal(strata,outlierMask,allDates,groupCode,calendarDOY,yearAxis);

fprintf('Calendar grouping: %d weekday/day-type groups; %d year group(s); holidays=%d.\n', ...
    nGroups,nYears,nnz(holidayMask));
for dt=1:nGroups
    fprintf('  %-16s %d dates\n',char(groupNames(dt)),nnz(groupCode==dt));
end
if o.MergeYears
    fprintf('  Years pooled: %d-%d\n',min(yearAxis),max(yearAxis));
else
    fprintf('  Years separated: %s\n',strjoin(string(yearGroups),', '));
end
info=struct('faultCleaning',faultInfo,'outerTestDates',outerTestDates,'strataCriterion',string(o.StrataCriterion), ...
    'intervalMinutes',intervalMode,'slotMinutes',slotMinutes, ...
    'dayStartIndex',dayStart,'dayEndIndex',dayEnd,'allDates',allDates, ...
    'dayTypeCode',groupCode,'dayTypeNames',groupNames,'holidayMask',holidayMask, ...
    'yearAxis',yearAxis,'yearGroups',yearGroups,'mergeYears',o.MergeYears,'mergeWeekdays',o.MergeWeekdays, ...
    'timeSegmentation',{timeSegmentation},'seasonSegmentation',{seasonSegmentation}, ...
    'timeWindowHours',o.TimeWindowHours,'minPoissonSegmentMinutes',o.MinPoissonSegmentMinutes, ...
    'maxTimeSegmentsPerDayType',o.MaxPoissonSegmentsPerWeekday, ...
    'minSeasonSegmentDays',o.MinSeasonSegmentDays, ...
    'maxSeasonSegmentsPerTimeInterval',o.MaxSeasonSegmentsPerTimeInterval, ...
    'outlierRemoval',outlierInfo,'fanoRefinement',refinementInfo, ...
    'excludeMissingSensorsPerSlice',o.ExcludeMissingSensorsPerSlice,'sliceSensorMinCoverage',o.SliceSensorMinCoverage, ...
    'measurementSectionGroups',{sectionGroups},'laneGroupSizes',cellfun(@numel,sectionGroups(:)), ...
    'modelSensorMask',modelSensorMask, ...
    'laneRequireAll',o.LaneMergeRequireAllLanes);
end

function sq=makeCalendarStratumLocal(X,allDates,groupCode,calendarDOY,yearAxis,dt,yearValue,a,b,s0,s1, ...
    slotMinutes,intervalMode,groupNames,timeBIC,seasonBIC,refDepth,refinedForFano,o,modelSensorMask)
yearMask=(yearValue==0) | (yearAxis==yearValue);
dateIdx=find(groupCode==dt & yearMask & calendarDOY>=s0 & calendarDOY<=s1);
if isempty(dateIdx), sq=[]; return; end
D=numel(dateIdx); K=b-a+1; m=size(X,3);
block=X(dateIdx,a:b,:);
Yfull=reshape(block,D*K,m);
obsDates=repmat(allDates(dateIdx),K,1);
obsSlots=kron((a:b).',ones(D,1));
valid=isfinite(Yfull) & Yfull>=0;
allSensorCompleteOriginal=nnz(all(valid,2));
modelSensorMask=logical(modelSensorMask(:)).';
if numel(modelSensorMask)~=m
    error('runParmaTrafficTomography:ModelSensorMask','Internal modelSensorMask width mismatch.');
end
sensorCoverage=mean(valid,1);

% Missing-data rule. We never pass NaN/negative entries to the estimator.
% The old implementation required 100%% completeness over an entire slice,
% which can reject many otherwise useful sensors because of one isolated gap.
% The new default retains ESTIMATION sensors with high per-slice coverage, then
% removes only observation rows that are incomplete across those retained sensors.
% Validation sensors never influence either the retained-sensor mask or row filter.
if o.ExcludeMissingSensorsPerSlice
    sensorMask=modelSensorMask & sensorCoverage>=o.SliceSensorMinCoverage;
else
    sensorMask=modelSensorMask;
end
% Nearly constant sensors (stuck within this stratum): Fano factor far below 1.
if isfield(o,'MinSensorFanoPerStratum') && o.MinSensorFanoPerStratum>0
    Yv=Yfull; Yv(~valid)=NaN;
    mu_s=mean(Yv,1,'omitnan'); var_s=var(Yv,0,1,'omitnan');
    nearlyConstant=mu_s>=5 & var_s./max(mu_s,eps)<o.MinSensorFanoPerStratum;
    sensorMask=sensorMask & ~nearlyConstant;
end
if ~any(sensorMask), sq=[]; return; end
complete=all(valid(:,sensorMask),2);
rowsRetained=nnz(complete);
if rowsRetained==0, sq=[]; return; end
Yfull=Yfull(complete,:); obsDates=obsDates(complete); obsSlots=obsSlots(complete);
Y=round(Yfull(:,sensorMask));
Yfull=round(Yfull);
rowRetentionFraction=rowsRetained/max(D*K,1);
obsTimes=obsDates+minutes(slotMinutes(obsSlots));
[~,ord]=sort(obsTimes); Y=Y(ord,:); Yfull=Yfull(ord,:); obsDates=obsDates(ord); obsSlots=obsSlots(ord);
dtName=char(groupNames(dt));
if yearValue==0, yearLabel='AllYears'; else, yearLabel=sprintf('%d',yearValue); end
sq=struct( ...
    'DayTypeCode',dt,'DayType',dtName,'Weekday',dt,'WeekdayName',dtName, ...
    'Year',yearValue,'YearLabel',yearLabel, ...
    'StartMinute',slotMinutes(a),'EndMinute',slotMinutes(b)+intervalMode, ...
    'DurationMinutes',(b-a+1)*intervalMode,'ObservationMinutes',intervalMode, ...
    'SeasonStartDOY',s0,'SeasonEndDOY',s1, ...
    'SeasonStart',canonicalMonthDayLabel(s0),'SeasonEnd',canonicalMonthDayLabel(s1),'SeasonLabel','', ...
    'Y',Y,'YFull',Yfull,'SensorMask',logical(sensorMask),'SensorRows',find(sensorMask(:)),'SensorCoverage',sensorCoverage,'Dates',obsDates,'NativeSlots',(a:b).','Slots',obsSlots(:), ...
    'BIC',seasonBIC,'TimeBIC',timeBIC,'SeasonBIC',seasonBIC, ...
    'CandidateDates',D,'CompleteDates',numel(unique(obsDates)), ...
    'CandidateObservations',D*K,'CompleteObservations',size(Y,1), ...
    'AllSensorCompleteObservations',allSensorCompleteOriginal, ...
    'RowsRetainedAfterMissingFilter',rowsRetained,'RowRetentionFraction',rowRetentionFraction, ...
    'RefinementDepth',refDepth,'RefinedForFano',logical(refinedForFano));
sq.SeasonLabel=sprintf('%s-%s',sq.SeasonStart,sq.SeasonEnd);
end

function [Xclean,outlierMask,info]=removeCalendarOutliersLocal( ...
    X,allDates,groupCode,calendarDOY,yearAxis,strata,o,modelSensorMask,isTestDay)
% The median and MAD of each score are learned on non-test days; the rule is
% then applied to all days.
if nargin<9 || isempty(isTestDay), isTestDay=false(size(X,1),1); end
Xclean=X; outlierMask=false(size(X)); rowsAffected=false(size(X,1),size(X,2)); valuesRemoved=0;
for q=1:numel(strata)
    st=strata(q); yearMask=(st.Year==0) | (yearAxis==st.Year);
    dateIdx=find(groupCode==st.DayTypeCode & yearMask & ...
        calendarDOY>=st.SeasonStartDOY & calendarDOY<=st.SeasonEndDOY);
    slots=double(st.NativeSlots(:));
    if isempty(dateIdx) || isempty(slots), continue; end
    block=Xclean(dateIdx,slots,:); D=numel(dateIdx); K=numel(slots); m=size(X,3);
    B=reshape(block,D*K,m); F=false(size(B));
    learnRow=repmat(~isTestDay(dateIdx(:)),K,1);      % rows of B: dates vary fastest
    for j=reshape(find(modelSensorMask),1,[])
        v=B(:,j); good=isfinite(v) & v>=0; ng=nnz(good);
        if ng<o.OutlierMinSamples || nnz(good & learnRow)<o.OutlierMinSamples, continue; end
        zl=log1p(v(good & learnRow)); med=median(zl,'omitnan'); mad0=median(abs(zl-med),'omitnan'); sigma=1.4826*mad0;
        z=log1p(v(good));
        if ~(isfinite(sigma) && sigma>sqrt(eps))
            medCount=median(v(good & learnRow),'omitnan'); sigma=sqrt(max(medCount,1))/max(1+medCount,1);
        end
        if ~(isfinite(sigma) && sigma>sqrt(eps)), continue; end
        rz=abs(z-med)/sigma; flagGood=rz>o.OutlierMADThreshold; nFlag=nnz(flagGood);
        if nFlag==0, continue; end
        maxRemove=floor(o.OutlierMaxFractionPerSensor*ng);
        if o.OutlierMaxFractionPerSensor>0 && maxRemove<1, maxRemove=1; end
        if o.OutlierMaxFractionPerSensor==0
            flagGood(:)=false;
        elseif nFlag>maxRemove
            idxGood=find(good); candidate=idxGood(flagGood); score=rz(flagGood); [~,ord]=sort(score,'descend');
            keep=candidate(ord(1:maxRemove)); flagGood(:)=false; [~,loc]=ismember(keep,idxGood); flagGood(loc)=true;
        end
        idxGood=find(good); F(idxGood(flagGood),j)=true;
    end
    if any(F(:))
        maskBlock=reshape(F,D,K,m); existing=outlierMask(dateIdx,slots,:);
        outlierMask(dateIdx,slots,:)=existing | maskBlock;
        rowMask=reshape(any(F,2),D,K); rowsAffected(dateIdx,slots)=rowsAffected(dateIdx,slots) | rowMask;
        valuesRemoved=valuesRemoved+nnz(F);
    end
end
Xclean(outlierMask)=NaN;
info=struct('enabled',true,'valuesRemoved',valuesRemoved, ...
    'observationRowsAffected',nnz(rowsAffected),'preliminaryStrata',numel(strata), ...
    'threshold',o.OutlierMADThreshold,'minSamples',o.OutlierMinSamples, ...
    'maxFractionPerSensor',o.OutlierMaxFractionPerSensor,'modelSensorMask',modelSensorMask);
fprintf('Outlier removal (estimation sensors only): %d sensor counts removed across %d native observation rows (threshold %.2f log-MAD).\n', ...
    info.valuesRemoved,info.observationRowsAffected,o.OutlierMADThreshold);
end

function strata=annotateOutlierStatsLocal(strata,outlierMask,allDates,groupCode,calendarDOY,yearAxis)
for q=1:numel(strata)
    st=strata(q); yearMask=(st.Year==0) | (yearAxis==st.Year);
    dateIdx=find(groupCode==st.DayTypeCode & yearMask & ...
        calendarDOY>=st.SeasonStartDOY & calendarDOY<=st.SeasonEndDOY);
    slots=double(st.NativeSlots(:));
    if isempty(dateIdx) || isempty(slots) || isempty(outlierMask)
        nVal=0; nRow=0; frac=0;
    else
        B=outlierMask(dateIdx,slots,:); nVal=nnz(B);
        RR=reshape(B,numel(dateIdx)*numel(slots),size(outlierMask,3)); nRow=nnz(any(RR,2));
        denom=max(numel(dateIdx)*numel(slots)*size(outlierMask,3),1); frac=nVal/denom;
    end
    strata(q).OutlierValuesRemoved=nVal;
    strata(q).OutlierObservationRowsAffected=nRow;
    strata(q).OutlierFractionOfCandidateValues=frac;
end
end

function [strata,info]=refineCalendarStrataByFanoLocal( ...
    X,allDates,groupCode,calendarDOY,yearAxis,slotMinutes,intervalMode,groupNames,strata,o,modelSensorMask)
info=struct('enabled',true,'splits',0,'passes',0, ...
    'targetMedianFano',o.MaxMedianFanoPerStratum,'unresolved',0);
minSlots=max(1,ceil(o.MinPoissonSegmentMinutes/intervalMode));
for pass=1:o.FanoRefinementMaxPasses
    didSplit=false; out=repmat(strata(1),0,1);
    for q=1:numel(strata)
        parent=strata(q); f0=medianFanoForModelSensorsLocal(parent,modelSensorMask,o.FanoMinMeanCount);
        if ~(isfinite(f0) && f0>o.MaxMedianFanoPerStratum) || size(parent.Y,1)<2*o.FanoRefinementMinSamples
            out(end+1,1)=parent; continue; %#ok<AGROW>
        end
        a=parent.NativeSlots(1); b=parent.NativeSlots(end); dt=parent.DayTypeCode;
        yv=parent.Year; s0=parent.SeasonStartDOY; s1=parent.SeasonEndDOY;
        yearMask=(yv==0) | (yearAxis==yv);
        bestScore=Inf; bestChildren=[];

        if (b-a+1)>=2*minSlots
            dateIdx=find(groupCode==dt & yearMask & calendarDOY>=s0 & calendarDOY<=s1);
            C=zeros(b-a+1,1); E=zeros(b-a+1,1);
            for u=1:numel(C)
                v=squeeze(X(dateIdx,a+u-1,modelSensorMask)); good=isfinite(v) & v>=0;
                C(u)=sum(v(good)); E(u)=nnz(good);
            end
            cutRel=bestForcedPoissonSplitLocal(C,E,minSlots);
            if cutRel>0
                cut=a+cutRel-1;
                c1=makeCalendarStratumLocal(X,allDates,groupCode,calendarDOY,yearAxis,dt,yv,a,cut,s0,s1, ...
                    slotMinutes,intervalMode,groupNames,parent.TimeBIC,parent.SeasonBIC,parent.RefinementDepth+1,true,o,modelSensorMask);
                c2=makeCalendarStratumLocal(X,allDates,groupCode,calendarDOY,yearAxis,dt,yv,cut+1,b,s0,s1, ...
                    slotMinutes,intervalMode,groupNames,parent.TimeBIC,parent.SeasonBIC,parent.RefinementDepth+1,true,o,modelSensorMask);
                [score,ok]=childFanoScoreLocal(c1,c2,o,modelSensorMask);
                if ok && score<bestScore, bestScore=score; bestChildren=[c1;c2]; end %#ok<AGROW>
            end
        end

        if (s1-s0+1)>=2*o.MinSeasonSegmentDays
            dateIdx=find(groupCode==dt & yearMask & calendarDOY>=s0 & calendarDOY<=s1);
            C=zeros(s1-s0+1,1); E=zeros(s1-s0+1,1);
            for kk=1:numel(dateIdx)
                d=dateIdx(kk); cd=calendarDOY(d)-s0+1; v=squeeze(X(d,a:b,modelSensorMask)); good=isfinite(v) & v>=0;
                C(cd)=C(cd)+sum(v(good)); E(cd)=E(cd)+nnz(good);
            end
            cutRel=bestForcedPoissonSplitLocal(C,E,o.MinSeasonSegmentDays);
            if cutRel>0
                cut=s0+cutRel-1;
                c1=makeCalendarStratumLocal(X,allDates,groupCode,calendarDOY,yearAxis,dt,yv,a,b,s0,cut, ...
                    slotMinutes,intervalMode,groupNames,parent.TimeBIC,parent.SeasonBIC,parent.RefinementDepth+1,true,o,modelSensorMask);
                c2=makeCalendarStratumLocal(X,allDates,groupCode,calendarDOY,yearAxis,dt,yv,a,b,cut+1,s1, ...
                    slotMinutes,intervalMode,groupNames,parent.TimeBIC,parent.SeasonBIC,parent.RefinementDepth+1,true,o,modelSensorMask);
                [score,ok]=childFanoScoreLocal(c1,c2,o,modelSensorMask);
                if ok && score<bestScore, bestScore=score; bestChildren=[c1;c2]; end %#ok<AGROW>
            end
        end

        if ~isempty(bestChildren) && bestScore <= f0*(1-o.FanoRefinementMinImprovement)
            out(end+1,1)=bestChildren(1); out(end+1,1)=bestChildren(2); %#ok<AGROW>
            info.splits=info.splits+1; didSplit=true;
        else
            out(end+1,1)=parent; %#ok<AGROW>
        end
    end
    strata=out; info.passes=pass; if ~didSplit, break; end
end
f=zeros(numel(strata),1); for q=1:numel(strata), f(q)=medianFanoForModelSensorsLocal(strata(q),modelSensorMask,o.FanoMinMeanCount); end
info.unresolved=nnz(isfinite(f) & f>o.MaxMedianFanoPerStratum);
fprintf('Fano refinement: %d additional splits in %d pass(es); %d strata remain above target %.2f.\n', ...
    info.splits,info.passes,info.unresolved,o.MaxMedianFanoPerStratum);
end

function [score,ok]=childFanoScoreLocal(c1,c2,o,modelSensorMask)
ok=~isempty(c1) && ~isempty(c2) && ...
    size(c1.Y,1)>=o.FanoRefinementMinSamples && size(c2.Y,1)>=o.FanoRefinementMinSamples;
if ~ok, score=Inf; return; end
f1=medianFanoForModelSensorsLocal(c1,modelSensorMask,o.FanoMinMeanCount); f2=medianFanoForModelSensorsLocal(c2,modelSensorMask,o.FanoMinMeanCount);
if ~isfinite(f1) || ~isfinite(f2), score=Inf; ok=false; return; end
n1=size(c1.Y,1); n2=size(c2.Y,1); score=(n1*f1+n2*f2)/(n1+n2);
end

function f=medianFanoForModelSensorsLocal(st,modelSensorMask,minMean)
rows=double(st.SensorRows(:));
if isempty(rows), f=NaN; return; end
keep=modelSensorMask(rows);
if ~any(keep), f=NaN; return; end
f=medianFanoLocal(st.Y(:,keep),minMean);
end

function f=medianFanoLocal(Y,minMean)
if isempty(Y) || size(Y,1)<2, f=NaN; return; end
mu=mean(Y,1).'; vv=var(Y,0,1).'; good=isfinite(mu) & isfinite(vv) & mu>=minMean;
if ~any(good), f=NaN; return; end
ff=vv(good)./max(mu(good),eps); f=median(ff,'omitnan');
end

function cut=bestForcedPoissonSplitLocal(C,E,minWidth)
C=C(:); E=E(:); T=numel(C); cut=0; if T<2*minWidth, return; end
cs=[0;cumsum(C)]; es=[0;cumsum(E)]; best=Inf;
for h=minWidth:T-minWidth
    c1=cs(h+1); e1=es(h+1); c2=cs(end)-cs(h+1); e2=es(end)-es(h+1);
    if e1<=0 || e2<=0, continue; end
    v=poissonProfileCostLocal(c1,e1)+poissonProfileCostLocal(c2,e2);
    if v<best, best=v; cut=h; end
end
end

function v=poissonProfileCostLocal(c,e)
if e<=0, v=Inf; elseif c>0, v=c-c*log(c/e); else, v=0; end
end

function [S,N,Q]=slotStatsLocal(Xd,wantSquares)
% Per-slot sums, counts (and sums of squares) of finite counts over days.
% Xd: days x slots x sensors -> slots x sensors.
good=isfinite(Xd)&Xd>=0; Xz=Xd; Xz(~good)=0;
S=reshape(sum(Xz,1),size(Xd,2),size(Xd,3));
N=reshape(sum(good,1),size(Xd,2),size(Xd,3));
Q=[];
if wantSquares, Q=reshape(sum(Xz.^2,1),size(Xd,2),size(Xd,3)); end
end

function [bounds,bestCost] = cvPiecewiseSegments(Str,Ntr,Sv,Qv,Nv,A,minLen,maxK,tol,workers,blockSize,cnt,minCnt)
% Contiguous segmentation of positions (time slots or calendar days) that
% minimizes the validation squared error of sensor means predicted by an NNLS
% fit of the training means of each segment. Inputs are positions x sensors:
% training sums and counts, validation sums, sums of squares and counts. The
% number of segments (<= maxK) is the smallest whose total error is within a
% relative tolerance TOL (default 1%) of the best: extra segments that barely
% change the validation error are not kept (parsimony rule).
% Segments without training or validation data are not allowed. Positions
% outside the data range are attached to the first and last segments.
if nargin<9 || isempty(tol), tol=0.01; end
if nargin<10 || isempty(workers), workers=0; end
if nargin<11 || isempty(blockSize), blockSize=1; end
if nargin<12 || isempty(cnt), cnt=ones(size(Str,1),1); minCnt=0; end
P=size(Str,1); m=size(Str,2);
has=any(Ntr>0,2)|any(Nv>0,2);
first=find(has,1,'first'); last=find(has,1,'last');
if isempty(first), bounds=[1 P]; bestCost=NaN; return; end
r=first:last; Lp=numel(r);
% Boundaries on a grid of blockSize positions: sum the statistics by block
bs=max(1,round(blockSize)); L=ceil(Lp/bs); blk=ceil((1:Lp)'/bs);
agg=@(Z) reshape(accumarrayCols(blk,Z(r,:),L),L,m);
cStr=[zeros(1,m);cumsum(agg(Str),1)]; cNtr=[zeros(1,m);cumsum(agg(Ntr),1)];
cSv=[zeros(1,m);cumsum(agg(Sv),1)];   cQv=[zeros(1,m);cumsum(agg(Qv),1)];
cNv=[zeros(1,m);cumsum(agg(Nv),1)];
cCnt=[0;cumsum(accumarray(blk,cnt(r),[L 1]))];
minLen=max(1,min(ceil(minLen/bs),L)); maxK=max(1,min(maxK,floor(L/minLen)));
opt=optimset('Display','off');
cost=inf(L,L);
parfor (i=1:L, workers)                 % rows are independent (serial when workers=0)
    rc=inf(1,L);
    for j=i+minLen-1:L
        nt=cNtr(j+1,:)-cNtr(i,:); nv=cNv(j+1,:)-cNv(i,:);
        sens=nt>0;
        if ~any(sens) || ~any(nv>0) || cCnt(j+1)-cCnt(i)<minCnt, continue; end
        ybar=(cStr(j+1,sens)-cStr(i,sens))./nt(sens);
        mu=lsqnonneg(A(sens,:),ybar(:),opt);
        pr=(A*mu).'; ev=nv>0;
        sv=cSv(j+1,:)-cSv(i,:); qv=cQv(j+1,:)-cQv(i,:);
        rc(j)=sum(qv(ev)-2*pr(ev).*sv(ev)+pr(ev).^2.*nv(ev));
    end
    cost(i,:)=rc;
end
DP=inf(maxK,L); prev=zeros(maxK,L);
for j=minLen:L, DP(1,j)=cost(1,j); end
for k=2:maxK
    for j=k*minLen:L
        best=Inf; arg=0;
        for h=(k-1)*minLen:j-minLen
            v=DP(k-1,h)+cost(h+1,j);
            if v<best, best=v; arg=h; end
        end
        DP(k,j)=best; prev(k,j)=arg;
    end
end
bestCost=min(DP(:,L));
if ~isfinite(bestCost), bounds=[1 P]; bestCost=NaN; return; end
K=find(DP(:,L)<=(1+tol)*bestCost,1,'first'); bestCost=DP(K,L);
b=zeros(K,2); j=L;
for k=K:-1:1
    if k==1, h=0; else, h=prev(k,j); end
    b(k,:)=[h+1 j]; j=h;
end
% Block indices back to positions
b=[(b(:,1)-1)*bs+1, min(b(:,2)*bs,Lp)];
bounds=b+first-1; bounds(1,1)=1; bounds(end,2)=P;
end

function S=accumarrayCols(idx,Z,n)
% Column-wise sums of the rows of Z grouped by idx (1..n).
S=zeros(n,size(Z,2));
for c=1:size(Z,2), S(:,c)=accumarray(idx(:),Z(:,c),[n 1]); end
end

function [bounds,bestBIC] = poissonPiecewiseSegments(C,E,minSlots,maxK,cnt,minCnt)
% Exact dynamic programming for a piecewise-constant Poisson rate with
% exposure E_t. Segment MLE is lambda=sum(C)/sum(E); BIC selects K.
% Optional cnt (per position) and minCnt: every segment needs sum(cnt)>=minCnt.
C=C(:); E=E(:); T=numel(C);
if nargin<5 || isempty(cnt), cnt=ones(T,1); minCnt=0; end
cc=[0;cumsum(cnt(:))];
minSlots=max(1,min(minSlots,T));
maxK=max(1,min(maxK,floor(T/minSlots)));
cs=[0;cumsum(C)]; es=[0;cumsum(E)];
cost=inf(T,T);
for i=1:T
    for j=i:T
        if j-i+1<minSlots, continue; end
        c=cs(j+1)-cs(i); e=es(j+1)-es(i);
        if e<=0 || cc(j+1)-cc(i)<minCnt, continue; end
        if c>0
            cost(i,j)=c-c*log(c/e);
        else
            cost(i,j)=0;
        end
    end
end
DP=inf(maxK,T); prev=zeros(maxK,T);
for j=minSlots:T, DP(1,j)=cost(1,j); end
for k=2:maxK
    for j=k*minSlots:T
        best=Inf; arg=0;
        lo=(k-1)*minSlots;
        hi=j-minSlots;
        for h=lo:hi
            v=DP(k-1,h)+cost(h+1,j);
            if v<best, best=v; arg=h; end
        end
        DP(k,j)=best; prev(k,j)=arg;
    end
end
nExposure=max(sum(E),2);
bic=inf(maxK,1);
for k=1:maxK
    if isfinite(DP(k,T)), bic(k)=2*DP(k,T)+k*log(nExposure); end
end
[bestBIC,K]=min(bic);
if ~isfinite(bestBIC), bounds=[1 T]; bestBIC=NaN; return; end
bounds=zeros(K,2); j=T;
for k=K:-1:1
    if k==1, h=0; else, h=prev(k,j); end
    bounds(k,:)=[h+1 j]; j=h;
end
end

function [code,isHoliday,names] = classifyTrafficDayType(dates,o)
% Flexible weekday grouping. Holidays always override the nominal weekday.
dates=dates(:); wd=weekday(dates); isHoliday=false(size(dates));
if o.TreatHolidaysAsSunday, isHoliday=italianParmaHolidayMask(dates,o); end
if o.MergeWeekdays
    names=["Workday","Saturday","Sunday/Holiday"];
    code=ones(size(dates)); code(wd==7)=2; code(wd==1 | isHoliday)=3;
else
    names=["Monday","Tuesday","Wednesday","Thursday","Friday","Saturday","Sunday/Holiday"];
    code=zeros(size(dates));
    code(wd==2)=1; code(wd==3)=2; code(wd==4)=3; code(wd==5)=4;
    code(wd==6)=5; code(wd==7)=6; code(wd==1)=7;
    code(isHoliday)=7;
end
end

function tf = italianParmaHolidayMask(dates,o)
dates=dateshift(dates(:),'start','day');
y=year(dates); mo=month(dates); da=day(dates);
tf=false(size(dates));
fixed=[1 1;1 6;4 25;5 1;6 2;8 15;11 1;12 8;12 25;12 26];
for k=1:size(fixed,1)
    tf=tf | (mo==fixed(k,1) & da==fixed(k,2));
end
if o.IncludeParmaPatronHoliday
    tf=tf | (mo==1 & da==13);
end
% San Francesco d'Assisi: national holiday from 2026 onward.
tf=tf | (y>=2026 & mo==10 & da==4);
% Easter Monday varies by year.
yu=unique(y(:)).';
for yy=yu
    easter=easterSundayGregorian(yy);
    tf=tf | dates==(easter+days(1));
end
if ~isempty(o.AdditionalHolidayDates)
    extra=dateshift(o.AdditionalHolidayDates(:),'start','day');
    tf=tf | ismember(dates,extra);
end
end

function d = easterSundayGregorian(y)
% Meeus/Jones/Butcher Gregorian computus.
a=mod(y,19); b=floor(y/100); c=mod(y,100); d0=floor(b/4); e=mod(b,4);
f=floor((b+8)/25); g=floor((b-f+1)/3); h=mod(19*a+b-d0-g+15,30);
i=floor(c/4); k=mod(c,4); l=mod(32+2*e+2*i-h-k,7);
m=floor((a+11*h+22*l)/451);
mon=floor((h+l-7*m+114)/31);
daynum=mod(h+l-7*m+114,31)+1;
d=datetime(y,mon,daynum);
end

function idx = canonicalMonthDayIndex(dates)
% Month/day -> 1..366 on leap year 2000, allowing all years to be merged.
ref=datetime(2000,month(dates),day(dates));
idx=day(ref,'dayofyear');
idx=double(idx(:));
end

function s = canonicalMonthDayLabel(idx)
idx=max(1,min(366,round(double(idx))));
d=datetime(2000,1,1)+days(idx-1);
s=datestr(d,'mmm dd');
end

function s = weekdayName(wd)
names=["Sunday","Monday","Tuesday","Wednesday","Thursday","Friday","Saturday"];
s=char(names(wd));
end

function s = minuteLabel(m)
m=round(m); h=floor(mod(m,1440)/60); mm=mod(m,60);
if m>=1440, h=24; mm=0; end
s=sprintf('%02d:%02d',h,mm);
end


% =========================================================================
function [Y,dates,coverage,intervalMode,info] = aggregateDailyTraffic(fileName,sectionGroups,o)
sectionGroups=normalizeSectionGroupsLocal(sectionGroups);
rawSections=unique([sectionGroups{:}]);
S=load(fileName,'typeNames','dateStart','dateEnd','nRecords');
typeIdx=find(string(S.typeNames(:))==o.TrafficType,1);
if isempty(typeIdx)
    availableTypes = strjoin(string(S.typeNames(:)).', ', ');
    error('runParmaTrafficTomography:TrafficType', ...
        'Traffic type "%s" is not present. Available: %s',o.TrafficType,availableTypes);
end
M=matfile(fileName);
nRecords=double(S.nRecords);
if nRecords<1, error('runParmaTrafficTomography:NoRecords','Compact archive has no records.'); end

recType=double(M.recType(1:nRecords,1));
recDay=double(M.recDay(1:nRecords,1));
recSection=double(M.recSection(1:nRecords,1));
recInterval=double(M.recInterval(1:nRecords,1));
recStart=double(M.recSignalStart(1:nRecords,1));
recLength=double(M.recSignalLength(1:nRecords,1));

isType=recType==typeIdx;
isSec=ismember(recSection,rawSections);
base=isType & isSec & recLength>0 & isfinite(recInterval) & recInterval>0;
if ~any(base)
    error('runParmaTrafficTomography:NoTrafficRecords', ...
        'No %s records were found for the matched sections.',o.TrafficType);
end

% Date range: user supplied, otherwise all available days of this type.
lastAvailable=max(recDay(base));
firstAvailable=min(recDay(base));
if isempty(o.DateRange)
    dayEnd=lastAvailable;
    dayStart=firstAvailable;
else
    dr=o.DateRange(:);
    dayStart=max(1,round((dr(1)-S.dateStart)/days(1))+1);
    dayEnd=min(round((S.dateEnd-S.dateStart)/days(1))+1, ...
        round((dr(2)-S.dateStart)/days(1))+1);
end
if dayEnd<dayStart
    error('runParmaTrafficTomography:DateRange','Requested date range contains no records.');
end
allDays=(dayStart:dayEnd).';
allDates=S.dateStart+days(allDays-1);
if o.WeekdaysOnly
    wd=weekday(allDates);
    keepDay=wd>=2 & wd<=6;
    allDays=allDays(keepDay); allDates=allDates(keepDay);
end
if isempty(allDays), error('runParmaTrafficTomography:DateRange','No eligible days remain.'); end

% Use the modal positive aggregation interval among eligible records.
eligible=base & ismember(recDay,allDays);
vals=recInterval(eligible);
if isempty(vals), error('runParmaTrafficTomography:NoEligibleRecords','No eligible records in date range.'); end
intervalMode=mode(vals);
eligible=eligible & abs(recInterval-intervalMode)<1e-6;

q=numel(sectionGroups); nd=numel(allDays); nRaw=numel(rawSections);
Yraw=nan(nd,nRaw);
[~,dayPos]=ismember(recDay,allDays);
[~,secPos]=ismember(recSection,rawSections);
rows=find(eligible & dayPos>0 & secPos>0);
windowStart=o.TimeWindowHours(1)*60;
windowEnd=o.TimeWindowHours(2)*60;
usedRecords=0;

signalChunks=buildSignalReadChunks(rows,recStart,recLength,o.SignalReadChunkMB);
fprintf('Selection-window signal I/O: %d records -> %d chunk reads.\n',numel(rows),numel(signalChunks));
for cRead=1:numel(signalChunks)
    cFirst=signalChunks(cRead).first;
    cLast=signalChunks(cRead).last;
    signalChunk=M.signalData(cFirst:cLast,1); % one bounded contiguous HDF5 read
    rr=signalChunks(cRead).rows;
    for z=1:numel(rr)
        r=rr(z);
        L=recLength(r);
        if L<=0, continue; end
        localFirst=recStart(r)-cFirst+1;
        x=double(signalChunk(localFirst:localFirst+L-1));
        minutes=(0:L-1).'*intervalMode;
        ix=minutes>=windowStart & minutes<windowEnd;
        if ~any(ix), continue; end
        v=x(ix);
        good=isfinite(v) & v>=0;
        if nnz(good) < ceil(o.MinWindowCoverage*numel(v)), continue; end
        val=sum(v(good));
        d=dayPos(r); s=secPos(r);
        if isnan(Yraw(d,s))
            Yraw(d,s)=val;
        else
            % Duplicate same day/section/type record: average duplicates before
            % any across-lane summation so one raw lane is never double-counted.
            Yraw(d,s)=0.5*(Yraw(d,s)+val);
        end
        usedRecords=usedRecords+1;
    end
end

Yraw=markDeadDetectorsLocal(Yraw,2,o.MaxSensorZeroShare);
Yraw=markImplausibleDetectorsLocal(Yraw,2,o.MaxFlowPerHour*24);   % daily totals
Y=aggregateLaneGroups2DLocal(Yraw,rawSections,sectionGroups,o.LaneMergeRequireAllLanes);
Y=markImplausibleDetectorsLocal(Y,2,o.MaxSectionFlowPerHour*24);   % summed lanes, daily totals
coverage=sum(isfinite(Y),1).'/nd;
dates=allDates;
info=struct('dayStartIndex',dayStart,'dayEndIndex',dayEnd, ...
    'nEligibleDays',nd,'nUsedRecords',usedRecords,'trafficType',o.TrafficType, ...
    'timeWindowHours',o.TimeWindowHours,'weekdaysOnly',o.WeekdaysOnly, ...
    'laneGroupSizes',cellfun(@numel,sectionGroups(:)), ...
    'laneRequireAll',o.LaneMergeRequireAllLanes);
end


% =========================================================================
function [selected,completeRows] = selectObservedSections(matchSub,Y,o)
coverage=sum(isfinite(Y),1).'/size(Y,1);
quality=matchSub.MatchScore(:).*sqrt(max(coverage,0));

if o.UseAllSensors
    % Keep all accepted post-merge measurement rows with adequate coverage.
    % Multiple corsie have already been collapsed upstream; opposite directions
    % remain distinct because their lane-base names/directions differ.
    valid=matchSub.Accepted & coverage>=o.MinSectionCoverage;
    selected=find(valid);
    [~,ord]=sort(quality(selected),'descend');
    selected=selected(ord);

    if numel(selected)<o.MinObservedSections
        % If the preferred coverage cutoff leaves too few sensors, retain all
        % accepted sensors having any usable data rather than reverting to a
        % hard sensor-count cap.
        fallback=matchSub.Accepted & coverage>0;
        selected=find(fallback);
        [~,ord]=sort(quality(selected),'descend');
        selected=selected(ord);
    end

    if isempty(selected)
        completeRows=false(size(Y,1),1);
    else
        completeRows=all(isfinite(Y(:,selected)),2);
    end
    return
end

% Legacy capped mode: greedily retain unique physical segments while
% preserving at least MinCompleteDays common broad-window observations.
valid=coverage>=o.MinSectionCoverage & matchSub.Accepted;
idx=find(valid);
[~,ord]=sort(quality(idx),'descend');
idx=idx(ord);
selected=zeros(0,1);
completeRows=true(size(Y,1),1);
usedSegments=strings(0,1);

for q=1:numel(idx)
    j=idx(q);
    if any(usedSegments==matchSub.SegmentKey(j)), continue; end
    trialRows=completeRows & isfinite(Y(:,j));
    if nnz(trialRows)<o.MinCompleteDays, continue; end
    selected(end+1,1)=j; %#ok<AGROW>
    usedSegments(end+1,1)=matchSub.SegmentKey(j); %#ok<AGROW>
    completeRows=trialRows;
    if numel(selected)>=o.MaxObservedSections, break; end
end

if numel(selected)<o.MinObservedSections
    rest=find(matchSub.Accepted);
    rest=setdiff(rest,selected,'stable');
    while numel(selected)<o.MinObservedSections && ~isempty(rest)
        best=0; bestRows=[]; bestN=-1; bestQuality=-Inf;
        for q=1:numel(rest)
            j=rest(q);
            if any(usedSegments==matchSub.SegmentKey(j)), continue; end
            tr=completeRows & isfinite(Y(:,j));
            n=nnz(tr);
            if n>bestN || (n==bestN && quality(j)>bestQuality)
                best=j; bestRows=tr; bestN=n; bestQuality=quality(j);
            end
        end
        if best==0 || bestN<o.MinCompleteDays, break; end
        selected(end+1,1)=best; %#ok<AGROW>
        usedSegments(end+1,1)=matchSub.SegmentKey(best); %#ok<AGROW>
        completeRows=bestRows;
        rest(rest==best)=[];
    end
end
end


% =========================================================================
function [t,h]=edgeEndsLocal(G,e)
% Tail and head node indices of directed edge e.
en=G.Edges.EndNodes(e,:);
if isnumeric(en), t=double(en(1)); h=double(en(2));
else, t=findnode(G,en(1)); h=findnode(G,en(2)); end
end

function sig = generatePathSignatureClasses(G,measuredEdgeSets,o)
m=numel(measuredEdgeSets);
n=numnodes(G);
deg=indegree(G)+outdegree(G);
lat=G.Nodes.Lat; lon=G.Nodes.Lon;
latSpan=max(lat)-min(lat); lonSpan=max(lon)-min(lon);
nearBoundary = lat<=min(lat)+0.08*latSpan | lat>=max(lat)-0.08*latSpan | ...
               lon<=min(lon)+0.08*lonSpan | lon>=max(lon)-0.08*lonSpan;
if isfield(o,'PathTerminals') && strcmpi(o.PathTerminals,'all')
    terminal=(1:n).';
else
    terminal=find(nearBoundary | deg>=4);
    if numel(terminal)<20, terminal=(1:n).'; end
end
% Stochastic route choice: perturbed edge costs, redrawn every PathNoiseBatch paths.
baseWeight=G.Edges.Weight;
noise=0; batch=100;
if isfield(o,'PathCostNoise'), noise=o.PathCostNoise; end
if isfield(o,'PathNoiseBatch'), batch=o.PathNoiseBatch; end
Gp=G;

% Precompute the measured-edge -> sensor incidence once.  The previous code
% called ISMEMBER once per sensor for every candidate path, which becomes very
% expensive with all sensors and tens of thousands of shortest paths.  With
% this sparse lookup, one indexing operation gives the complete path signature.
nEdges=numedges(G);
nInc=sum(cellfun(@numel,measuredEdgeSets));
edgeIdx=zeros(nInc,1);
sensorIdx=zeros(nInc,1);
pos=0;
for i=1:m
    ids=unique(round(double(measuredEdgeSets{i}(:))));
    ids=ids(isfinite(ids) & ids>=1 & ids<=nEdges);
    ni=numel(ids);
    if ni==0, continue; end
    edgeIdx(pos+(1:ni))=ids;
    sensorIdx(pos+(1:ni))=i;
    pos=pos+ni;
end
edgeIdx=edgeIdx(1:pos);
sensorIdx=sensorIdx(1:pos);
if isempty(edgeIdx)
    sig=struct('Aall',zeros(m,0),'representativeEdgePaths',{{}}, ...
        'representativeNodePaths',{{}},'memberEdgePaths',{{}},'memberPathHits',{{}}, ...
        'hits',zeros(0,1),'acceptedPaths',0,'tries',0);
    return
end
edgeToSensor=sparse(edgeIdx,sensorIdx,1,nEdges,m)>0;
fprintf('Path signatures: sparse edge-to-sensor lookup with %d incidences over %d sensors.\n', ...
    nnz(edgeToSensor),m);

seen=containers.Map('KeyType','char','ValueType','double');        % signature -> class
pathSeen=cell(0,1);   % per class: edge path -> member index (small maps stay fast)
% Storage grows in chunks; member paths are kept in one global list and
% grouped by class at the end (appending to per-class cells copied them).
cap=1024; Aall=false(m,cap); repE=cell(cap,1); repN=cell(cap,1); hits=zeros(cap,1); nQ=0;
mcap=4096; memPath=cell(mcap,1); memClass=zeros(mcap,1); memHits=zeros(mcap,1); nM=0;
perSource=1; if isfield(o,'PathsPerSource'), perSource=max(1,round(o.PathsPerSource)); end
treeSrc=0; treeLeft=0; TRv=[]; TEv=[];                             % current shortest-path tree
accepted=0; tries=0; sinceNew=0; nSearches=0;
while accepted<o.NumPathSamples && tries<o.PathMaxTries && sinceNew<o.SignaturePatience
    tries=tries+1; sinceNew=sinceNew+1;
    if noise>0 && mod(tries-1,batch)==0
        Gp.Edges.Weight=baseWeight.*exp(noise*randn(size(baseWeight)));
        treeLeft=0;                                                 % costs changed
    end
    anchored=isfield(o,'SensorAnchoredShare') && rand<o.SensorAnchoredShare;
    if anchored
        % Sensor-anchored path: through one edge of a random sensor, then to a
        % random node (or from a random node to that edge).
        si=sensorIdx(randi(numel(sensorIdx)));
        cand=edgeIdx(sensorIdx==si);
        e0=cand(randi(numel(cand)));
        [eT,eH]=edgeEndsLocal(G,e0);
        if rand<0.5
            b=terminal(randi(numel(terminal)));
            if b==eH, continue; end
            [pn,dist,pe]=shortestpath(Gp,eH,b,'Method','positive'); nSearches=nSearches+1;
            if isempty(pe) || ~isfinite(dist), continue; end
            pe=[e0;double(pe(:))]; pn=[eT,pn(:).'];
        else
            a=terminal(randi(numel(terminal)));
            if a==eT, continue; end
            [pn,dist,pe]=shortestpath(Gp,a,eT,'Method','positive'); nSearches=nSearches+1;
            if isempty(pe) || ~isfinite(dist), continue; end
            pe=[double(pe(:));e0]; pn=[pn(:).',eH];
        end
    elseif perSource==1
        a=terminal(randi(numel(terminal))); b=terminal(randi(numel(terminal)));
        if a==b, continue; end
        [pn,dist,pe]=shortestpath(Gp,a,b,'Method','positive'); nSearches=nSearches+1;
        if isempty(pe) || ~isfinite(dist), continue; end
    else
        % One shortest-path tree per source serves PathsPerSource destinations.
        if treeLeft<=0
            treeSrc=terminal(randi(numel(terminal)));
            [TRv,~,TEv]=shortestpathtree(Gp,treeSrc,'OutputForm','vector','Method','positive');
            TRv=double(TRv(:)); TEv=double(TEv(:)); nSearches=nSearches+1;
            treeLeft=perSource;
        end
        treeLeft=treeLeft-1;
        b=terminal(randi(numel(terminal)));
        if b==treeSrc || ~(TRv(b)>0), continue; end
        [pn,pe]=treePathLocal(TRv,TEv,treeSrc,b);
        if isempty(pe), continue; end
    end

    % Entire m-dimensional signature in one sparse lookup.
    pe=double(pe(:));
    pe=pe(isfinite(pe) & pe>=1 & pe<=nEdges);
    if isempty(pe), continue; end
    signature=full(any(edgeToSensor(pe,:),1)).';
    signature=logical(signature);
    if ~any(signature), continue; end
    accepted=accepted+1;
    key=char('0'+signature.');
    pathKey=sprintf('%d,',pe);
    if isKey(seen,key)
        q=seen(key); hits(q)=hits(q)+1;
    else
        nQ=nQ+1; q=nQ;
        if q>cap
            cap=2*cap; Aall(:,cap)=false; repE{cap,1}=[]; repN{cap,1}=[]; hits(cap,1)=0;
        end
        seen(key)=q;
        Aall(:,q)=signature; repE{q}=pe; repN{q}=pn; hits(q)=1;
        pathSeen{q,1}=containers.Map('KeyType','char','ValueType','double');
        sinceNew=0;
    end
    mk=pathSeen{q};
    if isKey(mk,pathKey)
        k=mk(pathKey); memHits(k)=memHits(k)+1;
    else
        nM=nM+1;
        if nM>mcap
            mcap=2*mcap; memPath{mcap,1}=[]; memClass(mcap,1)=0; memHits(mcap,1)=0;
        end
        memPath{nM}=pe; memClass(nM)=q; memHits(nM)=1;
        mk(pathKey)=nM;                 % containers.Map is a handle: updates in place
    end
end
Aall=Aall(:,1:nQ); repE=repE(1:nQ); repN=repN(1:nQ); hits=hits(1:nQ);
% Group the member paths by class, in order of first appearance
memPath=memPath(1:nM); memClass=memClass(1:nM); memHits=memHits(1:nM);
memberE=cell(nQ,1); memberHits=cell(nQ,1);
[cls,ordM]=sort(memClass,'ascend');                % stable: keeps first-appearance order
bounds=[0;find(diff(cls));nM];
for t=1:numel(bounds)-1
    r=ordM(bounds(t)+1:bounds(t+1));
    q=cls(bounds(t)+1);
    memberE{q}=memPath(r); memberHits{q}=memHits(r);
end
fprintf('Shortest-path sampling: %d accepted paths in %d attempts (%d shortest-path searches), %d signatures.\n', ...
    accepted,tries,nSearches,nQ);
sig=struct('Aall',double(Aall),'representativeEdgePaths',{repE}, ...
    'representativeNodePaths',{repN},'memberEdgePaths',{memberE},'memberPathHits',{memberHits}, ...
    'hits',hits,'acceptedPaths',accepted,'tries',tries);
end

function [pn,pe]=treePathLocal(TR,TE,src,dst)
% Node and edge path from src to dst in a shortest-path tree given by the
% predecessor vector TR and the edge vector TE (TE(i): edge from TR(i) to i).
pn=dst; pe=zeros(0,1); v=dst; guard=numel(TR)+1;
while v~=src && guard>0
    u=TR(v);
    if ~(u>0), pn=[]; pe=[]; return; end
    pe(end+1,1)=TE(v); %#ok<AGROW>
    pn(end+1)=u; %#ok<AGROW>
    v=u; guard=guard-1;
end
if v~=src, pn=[]; pe=[]; return; end
pn=fliplr(pn); pe=flipud(pe);
end


function [A,sel,muAll,relResidual] = selectPathClasses(Aall,meanY,o)
Q=size(Aall,2);
if Q==0, A=[]; sel=[]; muAll=[]; relResidual=Inf; return; end
muAll=lsqnonneg(double(Aall),double(meanY),optimset('Display','off'));
fit=Aall*muAll;
relResidual=safeRelativeNorm(fit-meanY,meanY);

active=find(muAll>o.MinPathClassPilot);
if isempty(active)
    [~,ord]=sort(muAll,'descend');
    active=ord(1:min(Q,o.MaxPathClasses));
else
    [~,ord]=sort(muAll(active),'descend'); active=active(ord);
end
sel=active(1:min(numel(active),o.MaxPathClasses));

% If NNLS is extremely sparse, retain a few additional high-frequency-compatible
% signatures by their correlation with the positive measured mean.
minClasses=min([4,Q,o.MaxPathClasses]);
if numel(sel)<minClasses
    score=(Aall.'*meanY)./max(sum(Aall,1).',1);
    [~,ord]=sort(score,'descend');
    for q=1:numel(ord)
        if ~ismember(ord(q),sel), sel(end+1,1)=ord(q); end %#ok<AGROW>
        if numel(sel)>=minClasses, break; end
    end
end
% Coverage: a sensor crossed by no kept class has predicted mean zero whatever
% the fit, which inflates the misfit. Add, for each such sensor with positive
% mean, the crossing class with the largest pilot intensity (ties: largest
% mean-correlation score), on top of MaxPathClasses.
if ~isfield(o,'EnsureSensorCoverage') || o.EnsureSensorCoverage
    score=(Aall.'*meanY)./max(sum(Aall,1).',1);
    covered=any(Aall(:,sel)>0,2);
    for i=find(~covered & meanY(:)>0).'
        if any(Aall(i,sel)>0), continue; end          % covered by a class added meanwhile
        cand=find(Aall(i,:)>0);
        if isempty(cand), continue; end
        key=[muAll(cand) score(cand)];
        [~,ord]=sortrows(key,[-1 -2]);
        sel(end+1,1)=cand(ord(1)); %#ok<AGROW>
    end
end
A=double(Aall(:,sel));
% Safety: remove any duplicate or zero columns (should already be impossible).
keep=any(A,1);
A=A(:,keep); sel=sel(keep);
[~,ia]=unique(A.','rows','stable');
A=A(:,ia); sel=sel(ia);
end


% =========================================================================
function plotParmaMapping(G,match,selectedSections,lat,lon,names)
% Plot the downloaded OSM network, traffic sensor GPS positions, matched
% measured edges, and the inferred flow direction.  A V-shaped arrow head is
% drawn along every directed OSM edge assigned to a selected sensor.
try
    if exist('geoaxes','file')~=2 || exist('geoplot','file')~=2 || ...
            exist('geobasemap','file')~=2 || exist('geoscatter','file')~=2
        error('Geographic plotting functions are unavailable.');
    end

    f=figure('Name','Parma OSM map - measured traffic sensors and directions','Color','w');
    gx=geoaxes(f);

    % Prefer genuine OpenStreetMap raster tiles. If custom basemaps are not
    % available, fall back to MATLAB's streets basemap while retaining the
    % downloaded OSM vector graph overlay.
    basemapOK=false;
    if exist('addCustomBasemap','file')==2
        basemapName='parmaNetworkTomographyOSM';
        tileURL='https://a.tile.openstreetmap.org/${z}/${x}/${y}.png';
        attribution=[char(169) ' OpenStreetMap contributors'];
        try
            addCustomBasemap(basemapName,tileURL,'Attribution',attribution);
        catch MEadd
            msg=lower(MEadd.message);
            if ~(contains(msg,'already') || contains(msg,'exist') || ...
                    contains(msg,'duplicate') || contains(msg,'registered'))
                warning('runParmaTrafficTomography:OSMBasemap', ...
                    'Could not register OSM basemap: %s',MEadd.message);
            end
        end
        try
            geobasemap(gx,basemapName);
            basemapOK=true;
        catch
            basemapOK=false;
        end
    end
    if ~basemapOK
        geobasemap(gx,'streets');
    end
    hold(gx,'on');

    % Downloaded road graph.
    allEdges=(1:numedges(G)).';
    [la,lo]=edgeSegmentsLatLon(G,allEdges);
    geoplot(gx,la,lo,'Color',[0.55 0.55 0.55],'LineWidth',0.45, ...
        'HandleVisibility','off');

    % All accepted matches are shown faintly; final measured sections are
    % over-plotted with larger markers and directed road arrows.
    good=find(match.Accepted);
    hAccepted=geoscatter(gx,lat(good),lon(good),18,'o', ...
        'LineWidth',0.8,'DisplayName','Accepted matched sensors');
    %#ok<NASGU>

    hSelected=geoscatter(gx,lat(selectedSections),lon(selectedSections),48, ...
        'filled','DisplayName','Measured flow sensors');
    %#ok<NASGU>

    for k=1:numel(selectedSections)
        i=selectedSections(k);
        edgeIDs=match.DirectedEdgeIDs{i};
        if isempty(edgeIDs), continue; end

        % Plot all directed edge(s) assigned to this sensor. A sensor for
        % which direction is unresolved may legitimately contain both arcs.
        [le,loe]=edgeSegmentsLatLon(G,edgeIDs);
        hEdge=geoplot(gx,le,loe,'LineWidth',3.0, ...
            'DisplayName',sprintf('M%d: %s',k,char(names(i))));
        edgeColor=hEdge.Color;

        for q=1:numel(edgeIDs)
            [ahLat,ahLon]=directedArrowHeadLatLon(G,edgeIDs(q));
            if ~isempty(ahLat)
                geoplot(gx,ahLat,ahLon,'Color',edgeColor,'LineWidth',2.2, ...
                    'HandleVisibility','off');
            end
        end

        directionText=char(match.DirectionRule(i));
        if isempty(directionText), directionText='unknown'; end
        label=sprintf(' M%d: %s [%s]',k,char(names(i)),directionText);
        text(gx,lat(i),lon(i),label,'FontSize',8,'FontWeight','bold', ...
            'Interpreter','none');
    end

    % Tight geographic limits around the downloaded graph/sensors.
    latAll=G.Nodes.Lat; lonAll=G.Nodes.Lon;
    latMin=min([latAll(:);lat(selectedSections)]);
    latMax=max([latAll(:);lat(selectedSections)]);
    lonMin=min([lonAll(:);lon(selectedSections)]);
    lonMax=max([lonAll(:);lon(selectedSections)]);
    latPad=max(0.025*(latMax-latMin),2e-4);
    lonPad=max(0.025*(lonMax-lonMin),2e-4);
    geolimits(gx,[latMin-latPad latMax+latPad],[lonMin-lonPad lonMax+lonPad]);

    title(gx,'Parma: OSM road network, measured flow sensors, and inferred flow direction');
    subtitle(gx,'Arrow heads follow the directed OSM edge assigned from GPS, road name, and section direction cues');
    legend(gx,'Location','best');
    hold(gx,'off');

catch ME
    warning('runParmaTrafficTomography:PlotFallback', ...
        'Geographic map unavailable (%s). Using longitude/latitude fallback.',ME.message);

    figure('Name','Parma measured traffic sensors and directions - fallback','Color','w');
    pG=plot(G,'XData',G.Nodes.Lon,'YData',G.Nodes.Lat, ...
        'NodeLabel',{},'ArrowSize',1.5,'LineWidth',0.35);
    pG.EdgeColor=[0.75 0.75 0.75];
    pG.NodeColor=[0.82 0.82 0.82];
    hold on;
    scatter(lon(selectedSections),lat(selectedSections),45,'filled');

    for k=1:numel(selectedSections)
        i=selectedSections(k);
        edgeIDs=match.DirectedEdgeIDs{i};
        for q=1:numel(edgeIDs)
            [sNode,tNode]=edgeEndpointIndices(G,edgeIDs(q));
            x0=G.Nodes.Lon(sNode); y0=G.Nodes.Lat(sNode);
            dx=G.Nodes.Lon(tNode)-x0; dy=G.Nodes.Lat(tNode)-y0;
            quiver(x0,y0,dx,dy,0,'LineWidth',2.2,'MaxHeadSize',0.9);
        end
        label=sprintf(' M%d: %s [%s]',k,char(names(i)),char(match.DirectionRule(i)));
        text(lon(i),lat(i),label,'FontSize',8,'FontWeight','bold','Interpreter','none');
    end

    xlabel('Longitude'); ylabel('Latitude'); axis equal; grid on;
    title('Parma downloaded OSM graph: measured sensors and inferred traffic direction');
    hold off;
end
end


% =========================================================================
function [latArrow,lonArrow]=directedArrowHeadLatLon(G,edgeID)
% Return a small V-shaped geographic arrow head pointing from the source to
% the target node of one directed graph edge. Geometry is formed in a local
% tangent plane and transformed back to latitude/longitude.
[s,t]=edgeEndpointIndices(G,edgeID);
latS=G.Nodes.Lat(s); lonS=G.Nodes.Lon(s);
latT=G.Nodes.Lat(t); lonT=G.Nodes.Lon(t);

latRef=(latS+latT)/2;
ky=111320;
kx=111320*max(cosd(latRef),0.1);
dx=(lonT-lonS)*kx;
dy=(latT-latS)*ky;
L=hypot(dx,dy);
if ~isfinite(L) || L<1
    latArrow=[]; lonArrow=[]; return
end

ux=dx/L; uy=dy/L;
px=-uy; py=ux;

% Put the arrow head inside the segment rather than at a junction. This
% remains readable when several OSM segments meet at the target node.
tipFrac=0.72;
tipX=tipFrac*dx;
tipY=tipFrac*dy;
headLen=min(35,max(5,0.22*L));
wing=0.55*headLen;
baseX=tipX-headLen*ux;
baseY=tipY-headLen*uy;
leftX=baseX+wing*px;
leftY=baseY+wing*py;
rightX=baseX-wing*px;
rightY=baseY-wing*py;

x=[leftX tipX rightX];
y=[leftY tipY rightY];
latArrow=latS+y/ky;
lonArrow=lonS+x/kx;
end


% =========================================================================
function plotEstimatedClasses(G,pathEdges,mu,measuredEdgeSets,selectedSections,sectionNames)
try
    if exist('geoaxes','file')~=2 || exist('geoplot','file')~=2 || exist('geobasemap','file')~=2
        error('geoaxes unavailable');
    end
    f=figure('Name','Parma estimated path-signature classes','Color','w');
    gx=geoaxes(f); geobasemap(gx,'streets'); hold(gx,'on');
    allEdges=(1:numedges(G)).'; [la,lo]=edgeSegmentsLatLon(G,allEdges);
    geoplot(gx,la,lo,'Color',[0.75 0.75 0.75],'LineWidth',0.35,'HandleVisibility','off');

    [~,ord]=sort(mu,'descend');
    top=ord(1:min(8,numel(ord)));
    lw=1+4*mu(top)/max(max(mu(top)),eps);
    for q=1:numel(top)
        j=top(q); [lp,lop]=edgeSegmentsLatLon(G,pathEdges{j});
        geoplot(gx,lp,lop,'LineWidth',lw(q), ...
            'DisplayName',sprintf('Class %d, mu=%.2f',j,mu(j)));
    end
    obs=unique([measuredEdgeSets{:}]);
    [lm,lom]=edgeSegmentsLatLon(G,obs);
    geoplot(gx,lm,lom,'Color',[0.85 0.05 0.05],'LineWidth',3,'DisplayName','Measured sections');
    title(gx,'Parma: largest estimated path-signature intensities');
    legend(gx,'Location','best'); hold(gx,'off');
catch ME
    warning('runParmaTrafficTomography:PathPlotFallback','Estimated-path map unavailable: %s',ME.message);
end
end


% =========================================================================
function [latLine,lonLine]=edgeSegmentsLatLon(G,edgeIDs)
edgeIDs=edgeIDs(:);
if isempty(edgeIDs), latLine=nan(0,1); lonLine=nan(0,1); return; end
[s,t]=edgeEndpointIndices(G,edgeIDs);
latS=G.Nodes.Lat(s); latT=G.Nodes.Lat(t);
lonS=G.Nodes.Lon(s); lonT=G.Nodes.Lon(t);
latLine=reshape([latS.';latT.';nan(1,numel(edgeIDs))],[],1);
lonLine=reshape([lonS.';lonT.';nan(1,numel(edgeIDs))],[],1);
end


% =========================================================================
function [s,t]=edgeEndpointIndices(G,edgeIDs)
ends=G.Edges.EndNodes(edgeIDs,:);
if isnumeric(ends)
    s=double(ends(:,1)); t=double(ends(:,2));
else
    s=findnode(G,ends(:,1)); t=findnode(G,ends(:,2));
end
s=s(:); t=t(:);
if any(s<=0) || any(t<=0)
    error('runParmaTrafficTomography:EdgeEndpoints','Could not map graph edge endpoints.');
end
end


% =========================================================================
function value = getTag(tags,name,defaultValue)
if isstruct(tags) && isfield(tags,name), value=string(tags.(name));
else, value=string(defaultValue); end
end

function d = haversineKm(lat1,lon1,lat2,lon2)
R=6371.0088;
p1=deg2rad(lat1); p2=deg2rad(lat2);
dp=deg2rad(lat2-lat1); dl=deg2rad(lon2-lon1);
a=sin(dp/2).^2+cos(p1).*cos(p2).*sin(dl/2).^2;
d=2*R*atan2(sqrt(a),sqrt(max(0,1-a)));
end

function n = elementCount(E)
if iscell(E), n=numel(E); else, n=numel(E); end
end

function e = getElement(E,k)
if iscell(E), e=E{k}; else, e=E(k); end
end

function s = idKey(id)
s=sprintf('%.0f',double(id));
end

function h = simpleHash(text)
b=uint8(char(text)); h=uint32(2166136261);
for k=1:numel(b)
    h=bitxor(h,uint32(b(k)));
    h=uint32(mod(uint64(h)*uint64(16777619),2^32));
end
end

function chunks = buildSignalReadChunks(rows,recStart,recLength,chunkMB)
%BUILDSIGNALREADCHUNKS Group record signal slices into bounded contiguous reads.
% signalData is stored as SINGLE in the compact archive. Reading a larger
% contiguous block once is much faster than thousands of tiny matfile/HDF5
% requests. Gaps between selected records are included only while the complete
% span remains below the requested memory target.
rows=rows(:);
emptyChunk=struct('rows',zeros(0,1),'first',0,'last',-1);
chunks=repmat(emptyChunk,0,1);
if isempty(rows), return; end

starts=double(recStart(rows));
[~,ord]=sort(starts,'ascend');
rows=rows(ord);
maxValues=max(1,floor(double(chunkMB)*1024^2/4)); % 4 bytes per SINGLE value

firstRow=rows(1);
gFirst=double(recStart(firstRow));
gLast=gFirst+double(recLength(firstRow))-1;
gRows=firstRow;

for k=2:numel(rows)
    r=rows(k);
    rFirst=double(recStart(r));
    rLast=rFirst+double(recLength(r))-1;
    proposedLast=max(gLast,rLast);
    if proposedLast-gFirst+1 <= maxValues
        gRows(end+1,1)=r; %#ok<AGROW>
        gLast=proposedLast;
    else
        q=numel(chunks)+1;
        chunks(q,1)=struct('rows',gRows,'first',gFirst,'last',gLast); %#ok<AGROW>
        gFirst=rFirst; gLast=rLast; gRows=r;
    end
end
q=numel(chunks)+1;
chunks(q,1)=struct('rows',gRows,'first',gFirst,'last',gLast);
end


function out = mergeNameValueOptions(defaults,user)
if mod(numel(defaults),2)~=0 || mod(numel(user),2)~=0
    error('runParmaTrafficTomography:EstimatorOptions','EstimatorOptions must contain name/value pairs.');
end
out=defaults;
for k=1:2:numel(user)
    name=user{k};
    if ~(ischar(name)||isstring(name))
        error('runParmaTrafficTomography:EstimatorOptions','Estimator option names must be text.');
    end
    hit=[];
    for j=1:2:numel(out)
        if strcmpi(string(out{j}),string(name)), hit=j; break; end
    end
    if isempty(hit), out(end+1:end+2)=user(k:k+1); %#ok<AGROW>
    else, out{hit+1}=user{k+1}; end
end
end


function [estMask,valMask,info]=chooseSectionHoldoutLocal(Aall,m,o,hits)
if nargin<4 || isempty(hits), hits=ones(size(Aall,2),1); end
hits=double(hits(:));
minCoupling=0; if isfield(o,'MinValidationCoupling'), minCoupling=o.MinValidationCoupling; end
%CHOOSESECTIONHOLDOUTLOCAL Reproducible physical-section hold-out split.
% The automatic split uses topology only (Aall), never held-out flow values.
% It greedily removes candidate rows while preserving enough distinct nonzero
% signatures to support the configured path-class budget when feasible.

estMask=true(m,1); valMask=false(m,1);
fullDistinct=countDistinctNonzeroSignaturesLocal(Aall);
maxValidation=max(0,m-o.MinObservedSections);
if ~o.SectionHoldoutValidation || m<2 || maxValidation<1
    info=struct('Enabled',logical(o.SectionHoldoutValidation),'TargetValidationCount',0, ...
        'ActualValidationCount',0,'FullDistinctSignatures',fullDistinct, ...
        'EstimationDistinctSignatures',fullDistinct,'RequiredDistinctSignatures',min(fullDistinct,max(2,min(o.MaxPathClasses,o.MinObservedSections))), ...
        'RequestedRows',zeros(0,1),'Seed',o.SectionSplitSeed,'Reason',"hold-out disabled or too few measurements");
    return
end

if ~isempty(o.SectionValidationRows)
    requested=unique(double(o.SectionValidationRows(:)),'stable');
    if any(requested<1 | requested>m)
        error('runParmaTrafficTomography:SectionValidationRows', ...
            'SectionValidationRows must refer to measured rows 1..%d.',m);
    end
    if numel(requested)>maxValidation
        error('runParmaTrafficTomography:SectionValidationRows', ...
            ['The requested %d validation rows leave fewer than MinObservedSections=%d ', ...
             'measurements for estimation.'],numel(requested),o.MinObservedSections);
    end
    valMask(requested)=true; estMask=~valMask;
    estDistinct=countDistinctNonzeroSignaturesLocal(Aall(estMask,:));
    info=struct('Enabled',true,'TargetValidationCount',numel(requested), ...
        'ActualValidationCount',nnz(valMask),'FullDistinctSignatures',fullDistinct, ...
        'EstimationDistinctSignatures',estDistinct,'RequiredDistinctSignatures',min(fullDistinct,max(2,min(o.MaxPathClasses,o.MinObservedSections))), ...
        'RequestedRows',requested,'Seed',o.SectionSplitSeed,'Reason',"explicit SectionValidationRows");
    info.Coupling=validationCouplingLocal(Aall,hits,estMask,(1:m).');
    weak=requested(info.Coupling(requested)<minCoupling);
    if ~isempty(weak)
        warning('runParmaTrafficTomography:WeakValidationCoupling', ...
            ['Requested validation rows %s share less than %.0f%% of their sampled paths with ' ...
             'estimation sensors; their errors mostly test the path model.'],mat2str(weak(:).'),100*minCoupling);
    end
    return
end

target=max(o.SectionValidationMinCount,round(o.SectionValidationFraction*m));
target=min(maxValidation,max(0,target));
if target<1
    info=struct('Enabled',true,'TargetValidationCount',0,'ActualValidationCount',0, ...
        'FullDistinctSignatures',fullDistinct,'EstimationDistinctSignatures',fullDistinct, ...
        'RequiredDistinctSignatures',min(fullDistinct,max(2,min(o.MaxPathClasses,o.MinObservedSections))), ...
        'RequestedRows',zeros(0,1),'Seed',o.SectionSplitSeed, ...
        'Reason',"requested fraction yields no feasible validation rows");
    return
end

requiredDistinct=min(fullDistinct,max(2,min(o.MaxPathClasses,o.MinObservedSections)));
state=rng;
cleanup=onCleanup(@()rng(state)); %#ok<NASGU>
rng(o.SectionSplitSeed,'twister');
order=randperm(m);
for k=1:numel(order)
    if nnz(valMask)>=target, break; end
    r=order(k);
    trialEst=estMask; trialEst(r)=false;
    if nnz(trialEst)<o.MinObservedSections, continue; end
    % The candidate, and every sensor already held out, must keep its traffic
    % linked to the estimation sensors that remain.
    trialVal=valMask; trialVal(r)=true;
    cpl=validationCouplingLocal(Aall,hits,trialEst,find(trialVal));
    if any(cpl<minCoupling), continue; end
    nDistinct=countDistinctNonzeroSignaturesLocal(Aall(trialEst,:));
    if nDistinct>=requiredDistinct
        estMask=trialEst; valMask(r)=true;
    end
end

actual=nnz(valMask);
estDistinct=countDistinctNonzeroSignaturesLocal(Aall(estMask,:));
reason="automatic topology-preserving split";
if actual<target
    reason=sprintf('topology-preserving split reached %d/%d requested validation rows',actual,target);
    warning('runParmaTrafficTomography:SectionHoldoutReduced', ...
        ['Requested %d validation measurements, but only %d could be held out while ', ...
         'preserving at least %d distinct candidate signatures and %d estimation sensors.'], ...
        target,actual,requiredDistinct,o.MinObservedSections);
end
info=struct('Enabled',true,'TargetValidationCount',target, ...
    'ActualValidationCount',actual,'FullDistinctSignatures',fullDistinct, ...
    'EstimationDistinctSignatures',estDistinct,'RequiredDistinctSignatures',requiredDistinct, ...
    'RequestedRows',zeros(0,1), ...
    'Seed',o.SectionSplitSeed,'Reason',string(reason));
info.Coupling=validationCouplingLocal(Aall,hits,estMask,(1:m).');
end

function c=validationCouplingLocal(Aall,hits,estMask,rows)
% Share of the sampled paths through each sensor in ROWS (weighted by the
% number of times each path class was sampled) whose class also crosses at
% least one estimation sensor other than the sensor itself.
c=nan(numel(rows),1);
for k=1:numel(rows)
    j=rows(k);
    through=Aall(j,:)>0;
    if ~any(through), continue; end
    others=estMask; others(j)=false;
    linked=through & any(Aall(others,:)>0,1);
    c(k)=sum(hits(linked))/max(sum(hits(through)),eps);
end
end

function n=countDistinctNonzeroSignaturesLocal(A)
if isempty(A), n=0; return; end
sig=double(A.');
nonzero=any(sig~=0,2);
sig=sig(nonzero,:);
if isempty(sig), n=0; else, n=size(unique(sig,'rows','stable'),1); end
end

function [tf,reason] = sliceMeasurementIdentifiableLocal(Aq)
% Check whether the supplied dictionary has distinct nonzero signatures.
% Used both to audit original classes and to check the reduced fitting model.
tf=false; reason="";
if isempty(Aq) || size(Aq,1)==0 || size(Aq,2)==0
    reason="no retained sensors or no observable classes"; return
end
if any(sum(Aq~=0,1)==0)
    reason="one or more path classes are unobserved by the complete sensors"; return
end
sig=double(Aq.');
if size(unique(sig,'rows','stable'),1) < size(sig,1)
    reason="two or more path classes become indistinguishable after sensor exclusion"; return
end
tf=true;
end

function [mae,rmse,nrmse,bias,biasPercent] = flowErrorMetrics(predicted,measured)
%FLOWERRORMETRICS Error between predicted and measured mean sensor flows.
predicted=double(predicted(:));
measured=double(measured(:));
ok=isfinite(predicted) & isfinite(measured);
if ~any(ok)
    mae=NaN; rmse=NaN; nrmse=NaN; bias=NaN; biasPercent=NaN;
    return
end
predicted=predicted(ok); measured=measured(ok);
err=predicted-measured;
mae=mean(abs(err));
rmse=sqrt(mean(err.^2));
rmsMeasured=sqrt(mean(measured.^2));
nrmse=rmse/max(rmsMeasured,eps);
bias=mean(err);
biasPercent=100*bias/max(abs(mean(measured)),eps);
end


% =========================================================================
function e = safeRelativeNorm(delta,reference)
den=norm(reference);
if ~isfinite(den)||den<=eps, e=NaN; else, e=norm(delta)/den; end
end


% =========================================================================
function reduction=reduceSliceClassesLocal(Aq)
% Group original classes ONLY by their signatures on retained estimation rows.
% Validation observations and signatures never choose the allocation weights.
p=size(Aq,2); observed=any(Aq~=0,1).';
groupIndex=zeros(p,1); groupSize=zeros(p,1);
if any(observed)
    [signatures,~,g]=unique(double(Aq(:,observed).'),'rows');
    B=signatures.'; k=size(B,2);
    sizes=accumarray(g,1,[k 1]);
    originals=find(observed);
    groupIndex(originals)=g; groupSize(originals)=sizes(g);
    expansion=sparse(originals,g,1./sizes(g),p,k);
else
    B=zeros(size(Aq,1),0); expansion=sparse(p,0);
end
assert(isequal(size(Aq*expansion),size(B)) && ...
    norm(Aq*expansion-B,'fro')<=1e-10*max(1,norm(B,'fro')), ...
    'Invalid reduced-class expansion.');
reduction=struct('A',B,'Expansion',expansion, ...
    'GroupIndex',groupIndex,'GroupSize',groupSize,'Unobserved',~observed, ...
    'AllocationPolicy',"Equal original-class shares; unobserved classes set to zero by convention");
end


function out=fitStratumLocal(q,stratumQ,ctx)
% Fit one stratum (generated from the body of the stratum loop of
% runParmaTrafficTomography so that strata can be fitted with parfor).
% Per-stratum results are returned in OUT, one column or element each.
A=ctx.A;
Aq=ctx.Aq;
classReductionByStratum=ctx.classReductionByStratum;
estOptions=ctx.estOptions;
fq=ctx.fq;
i=ctx.i;
identifiableByStratum=ctx.identifiableByStratum;
m=ctx.m;
mEst=ctx.mEst;
mask=ctx.mask;
nativeInfo=ctx.nativeInfo;
o=ctx.o;
p=ctx.p;
projectionMatrixByStratum=ctx.projectionMatrixByStratum;
reduction=ctx.reduction;
sectionEstimationMask=ctx.sectionEstimationMask;
sectionValidationMask=ctx.sectionValidationMask;
Snum=1;
muHatByStratum=nan(p,Snum);
muRatePerHourByStratum=nan(p,Snum);
infoByStratum=cell(Snum,1);
trainMeanRelErr=nan(Snum,1);
testMeanRelErr=nan(Snum,1);
allMeanRelErr=nan(Snum,1);
firstOrderOptimality=nan(Snum,1);
estimationTime=nan(Snum,1);
nTrainVec=zeros(Snum,1);
nTestVec=zeros(Snum,1);
testDaysByStratum=zeros(Snum,1);     % distinct test days in the stratum
medianFanoByStratum=nan(Snum,1);

% Measured/predicted mean-flow diagnostics.  Each column is one
% day-type/time/date stratum and each row is one measured road section.
predMeanByStratum=nan(m,Snum);
trainMeasuredMeanByStratum=nan(m,Snum);
testMeasuredMeanByStratum=nan(m,Snum);
allMeasuredMeanByStratum=nan(m,Snum);
sectionValidationMeasuredMeanByStratum=nan(m,Snum);
sectionValidationPredictedMeanByStratum=nan(m,Snum);
sectionValidationObservationCountByStratum=zeros(m,Snum);
% Diagnostics added for the paper (baseline, Poisson checks, identifiability).
sectionValidationPredictedMeanNNLSByStratum=nan(m,Snum); % NNLS pilot (means only)
lag1AutocorrByStratum=nan(Snum,1);       % median lag-one autocorrelation over sensors
projectionRankByStratum=nan(Snum,1);     % jump-class rank of the projections used
rankConditionByStratum=false(Snum,1);    % rank == number of reduced classes
reducedClassCountFitted=nan(Snum,1);     % number of reduced classes in the fit
sensitivityConditionByStratum=nan(Snum,1); % condition number of the criterion Hessian
medianRelativeSEByStratum=nan(Snum,1);   % median sandwich SE / estimate (interior classes)
trainMisfitByStratum=nan(Snum,1);        % ||A mu - ybar||/||ybar|| of the estimate
nnlsMisfitByStratum=nan(Snum,1);         % same for the NNLS fit (best the classes allow)
exitflagByStratum=nan(Snum,1);           % trust-region exit flag (0: iteration limit)
downweightedFractionByStratum=nan(Snum,1); % share of projected counts downweighted
multistartAtBestByStratum=nan(Snum,1);   % share of starts reaching the best maximum
multistartDistinctByStratum=nan(Snum,1); % number of distinct local maxima found
multistartMaxRelDistByStratum=nan(Snum,1); % largest ||mu_k-muHat||/||muHat||
multistartGapByStratum=nan(Snum,1);      % criterion gap per observation, best - worst start
reliableByStratum=false(Snum,1);         % fitted and passing all reliability checks
reliabilityReason=strings(Snum,1);       % why a stratum is not reliable
muHatUnreliableByStratum=nan(p,Snum);    % estimates of rejected fits (audit only)

trainMAE=nan(Snum,1); trainRMSE=nan(Snum,1); trainNRMSE=nan(Snum,1);
trainBias=nan(Snum,1); trainBiasPercent=nan(Snum,1);
testMAE=nan(Snum,1); testRMSE=nan(Snum,1); testNRMSE=nan(Snum,1);
testBias=nan(Snum,1); testBiasPercent=nan(Snum,1);
allMAE=nan(Snum,1); allRMSE=nan(Snum,1); allNRMSE=nan(Snum,1);
allBias=nan(Snum,1); allBiasPercent=nan(Snum,1);
sectionValidationMAE=nan(Snum,1); sectionValidationRMSE=nan(Snum,1); sectionValidationNRMSE=nan(Snum,1);
sectionValidationBias=nan(Snum,1); sectionValidationBiasPercent=nan(Snum,1);

sensorCompleteByStratum=false(m,Snum);
sensorUsedByStratum=false(m,Snum);
predMeanNNLSByStratum=nan(m,Snum);       % NNLS prediction of all sensors (baseline)
nnlsRelResidualByStratum=nan(m,Snum);    % (predicted-observed)/observed of the stratum NNLS fit
for onceStratum_=1:1   % "continue" below skips the rest of this stratum
    completeMask=logical(stratumQ.SensorMask(:));
    sensorCompleteByStratum(:,1)=completeMask;
    mask=completeMask & sectionEstimationMask(:);
    YsFull=stratumQ.YFull;
    Ys=YsFull(:,mask);
    Ns=size(Ys,1);
    if nnz(mask)<o.MinObservedSections
        reliabilityReason(1)=sprintf('only %d complete estimation sensors',nnz(mask));
        warning('runParmaTrafficTomography:SliceSensors', ...
            'Skipping stratum %d: only %d/%d estimation sensors are complete in this slice.',q,nnz(mask),mEst);
        continue
    end
    if ~identifiableByStratum(q)
        reliabilityReason(1)="not identifiable after missing-sensor exclusion";
        warning('runParmaTrafficTomography:SliceIdentifiability', ...
            'Skipping stratum %d after missing-sensor exclusion: %s.',q,identifiabilityReason(q));
        continue
    end
    reduction=classReductionByStratum{q};
    Aq=reduction.A;
    if Ns<o.MinSamplesPerStratum
        reliabilityReason(1)=sprintf('N=%d < MinSamplesPerStratum=%d',Ns,o.MinSamplesPerStratum);
        warning('runParmaTrafficTomography:StratumSampleSize', ...
            'Skipping stratum %d: N=%d < MinSamplesPerStratum=%d.',q,Ns,o.MinSamplesPerStratum);
        continue
    end
    % Temporal split: hold out whole random days (default) or the last intervals.
    isTestRow=false(Ns,1);
    dq=stratumQ.Dates(:);
    if strcmpi(o.TemporalTestBlocks,'days') && isfield(nativeInfo,'outerTestDates') && numel(dq)==Ns
        % Common test days, fixed before stratification (same for both criteria)
        isTestRow=ismember(dq,nativeInfo.outerTestDates);
    elseif strcmpi(o.TemporalTestBlocks,'days') && numel(dq)==Ns && numel(unique(dq))>=5
        [~,~,gd]=unique(dq);
        nd=max(gd); nTestDays=max(1,round((1-o.TrainFraction)*nd));
        rsd=RandStream('mt19937ar','Seed',o.Seed+5000+q);
        isTestRow=ismember(gd,randperm(rsd,nd,nTestDays));
    else
        nLast=Ns-min(max(ceil(o.TrainFraction*Ns),max(5,ceil(0.6*Ns))),Ns-1);
        isTestRow(end-nLast+1:end)=true;
    end
    trainRows=find(~isTestRow); testRows=find(isTestRow);
    nTrain=numel(trainRows);
    if nTrain<2 || isempty(testRows)
        warning('runParmaTrafficTomography:StratumSplit', ...
            'Skipping stratum %d because it cannot be split into train/test data.',q);
        continue
    end
    YTrain=Ys(trainRows,:);
    YTest=Ys(testRows,:);

    % NNLS and likelihood use only complete sensors for this slice.
    muNNLSfull=lsqnonneg(Aq,mean(YTrain,1).',optimset('Display','off'));
    % Fit only as many classes as the stratum's data support: rank by the
    % stratum NNLS pilot, then restore coverage of every sensor with traffic.
    nAll=size(Aq,2);
    if isfinite(o.SamplesPerFittedClass), nKeep=max(1,floor(nTrain/o.SamplesPerFittedClass)); else, nKeep=nAll; end   % Inf: all classes
    if nKeep<nAll
        [~,ordC]=sort(muNNLSfull,'descend');
        keepC=false(nAll,1); keepC(ordC(1:nKeep))=true;
        mtr=mean(YTrain,1).';
        for i=find(~any(Aq(:,keepC)>0,2) & mtr>0).'
            if any(Aq(i,keepC)>0), continue; end
            cand=find(Aq(i,:)>0);
            if isempty(cand), continue; end
            [~,jb]=max(muNNLSfull(cand)); keepC(cand(jb))=true;
        end
        % Also keep, for every validation sensor, the crossing class with the
        % largest pilot intensity: otherwise its traffic may lie entirely on
        % dropped classes and be predicted as zero. Only the network structure
        % and the estimation-sensor pilot are used, not the validation counts.
        if any(sectionValidationMask)
            Aval=double(A(sectionValidationMask,:))*reduction.Expansion>0;   % validation rows x reduced classes
            for i=1:size(Aval,1)
                if any(Aval(i,keepC)), continue; end
                cand=find(Aval(i,:) & muNNLSfull(:).'>0);
                if isempty(cand), continue; end
                [~,jb]=max(muNNLSfull(cand)); keepC(cand(jb))=true;
            end
        end
    else
        keepC=true(nAll,1);
    end
    fittedClassIdx=find(keepC);
    AqFit=Aq(:,keepC);
    if nTrain<o.MinSamplesPerClass*size(AqFit,2)
        reliabilityReason(1)=sprintf('%d training intervals < %g x %d fitted classes',nTrain,o.MinSamplesPerClass,size(AqFit,2));
        warning('runParmaTrafficTomography:StratumSamplesPerClass', ...
            'Skipping stratum %d: %s.',q,reliabilityReason(1));
        continue
    end
    sensorUsedByStratum(:,1)=mask;
    nTrainVec(1)=size(YTrain,1); nTestVec(1)=size(YTest,1);
    if numel(dq)==Ns, testDaysByStratum(1)=numel(unique(dq(testRows))); end
    % Means-only baseline on the fitted classes (zero elsewhere).
    muNNLSraw=zeros(nAll,1);
    muNNLSraw(keepC)=lsqnonneg(AqFit,mean(YTrain,1).',optimset('Display','off'));
    mtrq=mean(YTrain,1).';
    nnlsRelResidualByStratum(mask,1)=(Aq*muNNLSraw-mtrq)./max(mtrq,eps);
    muInit=max(muNNLSraw(keepC),0.05);
    localOpts=estOptions;
    localOpts=mergeNameValueOptions(localOpts,{'PilotMu',muInit,'Mu0',muInit, ...
                                               'CandidateSeed',o.Seed+100+q,'DesignOnly',false});
    if ~isempty(projectionMatrixByStratum{q})
        localOpts=mergeNameValueOptions(localOpts,{'ProjectionMatrix',projectionMatrixByStratum{q}});
    end
    fprintf('  [%2d/%2d] %-14s %-8s %s-%s | %s-%s, N=%d, sensors=%d/%d ... ', ...
        q,Snum,stratumQ.DayType,stratumQ.YearLabel, ...
        minuteLabel(stratumQ.StartMinute),minuteLabel(stratumQ.EndMinute), ...
        stratumQ.SeasonStart,stratumQ.SeasonEnd,Ns,nnz(mask),mEst);
    try
        tEst=tic;
        [muFit,infoq]=estimateMuOptimizedProjections(AqFit,YTrain,localOpts{:});
        muFit=double(muFit(:));
        infoq.stratumEstimator=string(o.StratumEstimator); infoq.pclEstimate=muFit;
        if ~strcmpi(o.StratumEstimator,'pcl')
            tG=tic;
            try
                gOpts={};
                if strcmpi(o.GMMProjections,'pclDesign') && isfield(infoq,'C') && ~isempty(infoq.C)
                    Cg=unique([eye(size(AqFit,1)); double(infoq.C)],'rows','stable');
                    gOpts={'ProjectionMatrix',Cg,'MomentOrder',2,'CumulantStart',false};
                end
                gOpts=[gOpts {'WeightIterations',o.GMMWeightIterations}];
                if strcmpi(o.StratumEstimator,'robustGMM')
                    [muG,ig]=estimateMuEfficientGMM(AqFit,YTrain,'Robust',true,'Epsilon',o.GMMEpsilon,gOpts{:});
                else
                    [muG,ig]=estimateMuEfficientGMM(AqFit,YTrain,'Mu0',muFit,gOpts{:});   % PCL estimate as pilot
                end
                muG=double(muG(:));
                if numel(muG)==numel(muFit) && all(isfinite(muG)) && all(muG>=0)
                    muFit=muG; infoq.gmm=rmfield(ig,intersect(fieldnames(ig),{'C'}));
                    infoq.gmmExitflag=ig.exitflag;
                else
                    infoq.gmmError="invalid GMM estimate; PCL estimate kept";
                end
            catch MEg
                infoq.gmmError=string(MEg.message);
            end
            infoq.gmmSeconds=toc(tG);
        end
        assert(numel(muFit)==size(AqFit,2)&&all(isfinite(muFit)&muFit>=0), ...
            'Reduced estimator returned invalid class intensities.');
        muReduced=zeros(nAll,1); muReduced(keepC)=muFit;
        infoq.fittedClassIndices=fittedClassIdx;
        muq=reduction.Expansion*muReduced;
        % Equal splitting preserves predictions on the retained sensors.
        assert(norm(A(mask,:)*muq-Aq*muReduced,Inf)<=1e-9*max(1,norm(Aq*muReduced,Inf)), ...
            'Reduced-class expansion did not preserve fitted sensor predictions.');
        infoq.reducedMu=muReduced;
        infoq.expandedMu=muq;
        % Held-out covariance check of the latent flows: under the model, the
        % covariance of the counts at sensors i and j is the total rate of the
        % classes crossing both, (Aq*diag(mu)*Aq')_ij. NNLS fits the same means but
        % splits the flow differently, so the off-diagonal covariances of the test
        % counts test the split itself. Error: ||pred - obs||/||obs|| over pairs i<j.
        % The level of the observed covariances is inflated by overdispersion, which
        % no Poisson prediction reproduces; the scale-free error, with one scale
        % factor per stratum fitted by least squares, compares the pattern of which
        % sensor pairs covary, the part that depends on the split of the flow.
        cv=struct('pcl',NaN,'nnls',NaN,'pclScaled',NaN,'nnlsScaled',NaN,'pclScale',NaN, ...
            'nnlsScale',NaN,'pairs',0,'testRows',size(YTest,1));
        if size(YTest,1)>=3 && size(Aq,1)>=2
            So=cov(YTest); off=triu(true(size(So)),1); no=norm(So(off));
            if no>0
                Sp=Aq*diag(muReduced)*Aq.'; Sn=Aq*diag(muNNLSraw)*Aq.';
                p_=Sp(off); n_=Sn(off); o_=So(off);
                cv.pcl=norm(p_-o_)/no; cv.nnls=norm(n_-o_)/no; cv.pairs=nnz(off);
                if any(p_), cv.pclScale=(p_.'*o_)/(p_.'*p_); cv.pclScaled=norm(cv.pclScale*p_-o_)/no; end
                if any(n_), cv.nnlsScale=(n_.'*o_)/(n_.'*n_); cv.nnlsScaled=norm(cv.nnlsScale*n_-o_)/no; end
            end
        end
        infoq.covValidation=cv;
        infoq.classReduction=reduction;
        infoq.parameterIndexing="Reduced classes; expandedMu uses original classes";
        estimationTime(1)=toc(tEst);
        muHatByStratum(:,1)=muq;
        muRatePerHourByStratum(:,1)=muq/(stratumQ.ObservationMinutes/60);
        infoq.sensorMask=mask;
        infoq.sensorRows=find(mask);
        infoByStratum{1}=infoq;
        firstOrderOptimality(1)=infoq.firstOrderOptimality;
        projectionRankByStratum(1)=infoq.projectionDesignRank;
        reducedClassCountFitted(1)=size(AqFit,2);
        rankConditionByStratum(1)=infoq.projectionDesignRank>=size(AqFit,2);
        if isfield(infoq,'sensitivityCondition')
            sensitivityConditionByStratum(1)=infoq.sensitivityCondition;
        end
        if isfield(infoq,'downweightedFraction')
            downweightedFractionByStratum(1)=infoq.downweightedFraction;
        end
        if isfield(infoq,'multistart')
            ms=infoq.multistart;
            multistartAtBestByStratum(1)=ms.fractionAtBest;
            multistartDistinctByStratum(1)=ms.nDistinctMaxima;
            multistartMaxRelDistByStratum(1)=ms.maxRelativeDistance;
            multistartGapByStratum(1)=ms.criterionGapPerObservation;
        end
        why=strings(0,1);
        if o.RequireRankCondition && ~rankConditionByStratum(1)
            why(end+1)=sprintf('jump-class rank %d < %d classes',infoq.projectionDesignRank,size(AqFit,2));
        end
        if ~(sensitivityConditionByStratum(1)<=o.MaxSensitivityCondition)
            why(end+1)=sprintf('Hessian condition %.3g > %.3g',sensitivityConditionByStratum(1),o.MaxSensitivityCondition);
        end
        mtTrain=mean(YTrain,1).';
        fitRel=norm(Aq*muReduced-mtTrain)/max(norm(mtTrain),eps);
        nnlsRel=norm(Aq*muNNLSraw-mtTrain)/max(norm(mtTrain),eps);
        trainMisfitByStratum(1)=fitRel; nnlsMisfitByStratum(1)=nnlsRel;
        if isfield(infoq,'exitflag'), exitflagByStratum(1)=infoq.exitflag; end
        if ~(fitRel-nnlsRel<=o.MaxExcessMisfit)
            why(end+1)=sprintf('training-mean misfit %.3g vs %.3g for NNLS (excess > %.3g)', ...
                fitRel,nnlsRel,o.MaxExcessMisfit);
        end
        if ~(fitRel<=o.MaxTrainMeanRelError)
            why(end+1)=sprintf('training-mean misfit %.3g > %.3g',fitRel,o.MaxTrainMeanRelError);
        end
        if ~isempty(why) && exitflagByStratum(1)==0
            why(end+1)="trust region stopped at the iteration limit";
        end
        if isfield(infoq,'standardErrors') && ~isempty(infoq.standardErrors)
            se=infoq.standardErrors(:); in=isfinite(se) & muFit>0;   % fitted classes only
            if any(in), medianRelativeSEByStratum(1)=median(se(in)./muFit(in)); end
        end
        if ~isempty(why)
            % Keep the estimate for auditing, but treat the stratum as not
            % estimated in every summary, figure and validation statistic.
            reliabilityReason(1)=strjoin(why,'; ');
            muHatUnreliableByStratum(:,1)=muq;
            muHatByStratum(:,1)=NaN;
            muRatePerHourByStratum(:,1)=NaN;
            warning('runParmaTrafficTomography:UnreliableStratum', ...
                'Stratum %d is not reliable (%s); it is reported as not estimated.',q,reliabilityReason(1));
        else
            reliableByStratum(1)=true;
        end

        mt=mean(YTrain,1).'; mh=mean(YTest,1).'; ma=mean(Ys,1).'; mpEst=Aq*muReduced;
        mpAll=A*muq;
        trainMeasuredMeanByStratum(mask,1)=mt;
        testMeasuredMeanByStratum(mask,1)=mh;
        allMeasuredMeanByStratum(:,1)=mean(YsFull,1,'omitnan').';
        predMeanByStratum(:,1)=mpAll;
        predMeanNNLSByStratum(:,1)=A*(reduction.Expansion*muNNLSraw);

        trainMeanRelErr(1)=safeRelativeNorm(mpEst-mt,mt);
        testMeanRelErr(1)=safeRelativeNorm(mpEst-mh,mh);
        allMeanRelErr(1)=safeRelativeNorm(mpEst-ma,ma);
        [trainMAE(1),trainRMSE(1),trainNRMSE(1),trainBias(1),trainBiasPercent(1)] = ...
            flowErrorMetrics(mpEst,mt);
        [testMAE(1),testRMSE(1),testNRMSE(1),testBias(1),testBiasPercent(1)] = ...
            flowErrorMetrics(mpEst,mh);
        [allMAE(1),allRMSE(1),allNRMSE(1),allBias(1),allBiasPercent(1)] = ...
            flowErrorMetrics(mpEst,ma);

        % True section hold-out validation: these measured values were never used
        % in class screening, projection design, initialization, or likelihood fit.
        if any(sectionValidationMask)
            vrows=find(sectionValidationMask);
            V=YsFull(:,vrows);
            vn=sum(isfinite(V) & V>=0,1).';
            vm=mean(V,1,'omitnan').';
            vm(vn==0)=NaN;
            vp=mpAll(vrows);
            sectionValidationObservationCountByStratum(vrows,1)=vn;
            sectionValidationMeasuredMeanByStratum(vrows,1)=vm;
            sectionValidationPredictedMeanByStratum(vrows,1)=vp;
            % Baseline: the same held-out prediction from the NNLS fit of the
            % training means, which uses no distributional information.
            sectionValidationPredictedMeanNNLSByStratum(vrows,1)= ...
                A(vrows,:)*(reduction.Expansion*muNNLSraw);
            [sectionValidationMAE(1),sectionValidationRMSE(1),sectionValidationNRMSE(1), ...
                sectionValidationBias(1),sectionValidationBiasPercent(1)] = flowErrorMetrics(vp,vm);
        end

        fq=var(Ys,0,1).'./max(ma,eps);
        fq(ma<o.FanoMinMeanCount)=NaN;
        medianFanoByStratum(1)=median(fq,'omitnan');
        lag1AutocorrByStratum(1)=stratumLag1AutocorrLocal(Ys,stratumQ,ma,o.FanoMinMeanCount);
        if any(sectionValidationMask) && isfinite(sectionValidationNRMSE(1))
            fprintf('ok, temporal-test NRMSE %.3g; section-validation NRMSE %.3g\n', ...
                testNRMSE(1),sectionValidationNRMSE(1));
        else
            fprintf('ok, temporal-test RMSE %.3g, NRMSE %.3g, bias %.2f%%\n', ...
                testRMSE(1),testNRMSE(1),testBiasPercent(1));
        end
    catch ME
        reliabilityReason(1)="estimator error: "+string(ME.message);
        warning('runParmaTrafficTomography:StratumEstimationFailed', ...
            'Stratum %d (%s %s %s-%s | %s-%s) failed: %s',q,stratumQ.DayType,stratumQ.YearLabel, ...
            minuteLabel(stratumQ.StartMinute),minuteLabel(stratumQ.EndMinute), ...
            stratumQ.SeasonStart,stratumQ.SeasonEnd,ME.message);
        fprintf('FAILED\n');
    end
end
out=struct();
out.muHatByStratum=muHatByStratum;
out.muRatePerHourByStratum=muRatePerHourByStratum;
out.infoByStratum=infoByStratum;
out.trainMeanRelErr=trainMeanRelErr;
out.testMeanRelErr=testMeanRelErr;
out.allMeanRelErr=allMeanRelErr;
out.firstOrderOptimality=firstOrderOptimality;
out.estimationTime=estimationTime;
out.nTrainVec=nTrainVec;
out.nTestVec=nTestVec;
out.testDaysByStratum=testDaysByStratum;
out.medianFanoByStratum=medianFanoByStratum;
out.predMeanByStratum=predMeanByStratum;
out.trainMeasuredMeanByStratum=trainMeasuredMeanByStratum;
out.testMeasuredMeanByStratum=testMeasuredMeanByStratum;
out.allMeasuredMeanByStratum=allMeasuredMeanByStratum;
out.sectionValidationMeasuredMeanByStratum=sectionValidationMeasuredMeanByStratum;
out.sectionValidationPredictedMeanByStratum=sectionValidationPredictedMeanByStratum;
out.sectionValidationObservationCountByStratum=sectionValidationObservationCountByStratum;
out.sectionValidationPredictedMeanNNLSByStratum=sectionValidationPredictedMeanNNLSByStratum;
out.lag1AutocorrByStratum=lag1AutocorrByStratum;
out.projectionRankByStratum=projectionRankByStratum;
out.rankConditionByStratum=rankConditionByStratum;
out.reducedClassCountFitted=reducedClassCountFitted;
out.sensitivityConditionByStratum=sensitivityConditionByStratum;
out.medianRelativeSEByStratum=medianRelativeSEByStratum;
out.trainMisfitByStratum=trainMisfitByStratum;
out.nnlsMisfitByStratum=nnlsMisfitByStratum;
out.exitflagByStratum=exitflagByStratum;
out.downweightedFractionByStratum=downweightedFractionByStratum;
out.multistartAtBestByStratum=multistartAtBestByStratum;
out.multistartDistinctByStratum=multistartDistinctByStratum;
out.multistartMaxRelDistByStratum=multistartMaxRelDistByStratum;
out.multistartGapByStratum=multistartGapByStratum;
out.reliableByStratum=reliableByStratum;
out.reliabilityReason=reliabilityReason;
out.muHatUnreliableByStratum=muHatUnreliableByStratum;
out.trainMAE=trainMAE;
out.trainRMSE=trainRMSE;
out.trainNRMSE=trainNRMSE;
out.trainBias=trainBias;
out.trainBiasPercent=trainBiasPercent;
out.testMAE=testMAE;
out.testRMSE=testRMSE;
out.testNRMSE=testNRMSE;
out.testBias=testBias;
out.testBiasPercent=testBiasPercent;
out.allMAE=allMAE;
out.allRMSE=allRMSE;
out.allNRMSE=allNRMSE;
out.allBias=allBias;
out.allBiasPercent=allBiasPercent;
out.sectionValidationMAE=sectionValidationMAE;
out.sectionValidationRMSE=sectionValidationRMSE;
out.sectionValidationNRMSE=sectionValidationNRMSE;
out.sectionValidationBias=sectionValidationBias;
out.sectionValidationBiasPercent=sectionValidationBiasPercent;
out.sensorCompleteByStratum=sensorCompleteByStratum;
out.sensorUsedByStratum=sensorUsedByStratum;
out.predMeanNNLSByStratum=predMeanNNLSByStratum;
out.nnlsRelResidualByStratum=nnlsRelResidualByStratum;
end

function [sig,nAdded] = addLocalSensorClasses(sig,measuredEdgeSets,G)
% Append, for every sensor without one, a class crossing only that sensor. Its
% representative path is the sensor's own directed edge(s); it has one hit.
m = size(sig.Aall,1); nAdded = 0;
single = sum(sig.Aall~=0,1)==1;
have = false(m,1);
if any(single), [~,r] = max(sig.Aall(:,single)~=0,[],1); have(r) = true; end
for i = find(~have(:)).'
    e = unique(round(double(measuredEdgeSets{i}(:)))).';
    if isempty(e), continue; end
    col = zeros(m,1); col(i) = 1;
    sig.Aall(:,end+1) = col;
    nodes = [];
    try
        nodes = unique(reshape(double(G.Edges.EndNodes(e,:)).',1,[]),'stable');
    catch
    end
    sig.representativeEdgePaths{end+1} = e;
    sig.representativeNodePaths{end+1} = nodes;
    sig.memberEdgePaths{end+1} = {e};
    sig.memberPathHits{end+1} = 1;
    sig.hits(end+1) = 1;
    nAdded = nAdded + 1;
end
if size(sig.hits,1) > 1 && size(sig.hits,2) == 1, sig.hits = sig.hits(:); end
end

function match = applyEdgeOverridesLocal(match,cleanNames,file,G)
% Manual corrections of the sensor-to-edge matching (see EdgeOverrideCSV).
if ~isfile(file)
    error('runParmaTrafficTomography:EdgeOverrideFile','Edge override table not found: %s',file);
end
T = readtable(file,'TextType','string','Delimiter',',');
need = {'SectionName','Action'};
assert(all(ismember(need,T.Properties.VariableNames)), ...
    'runParmaTrafficTomography:EdgeOverrideColumns','The table needs the columns SectionName and Action (and EdgeIDs for "edges").');
if ~ismember('EdgeIDs',T.Properties.VariableNames), T.EdgeIDs = strings(height(T),1); end
if ~ismember('Street',T.Properties.VariableNames), T.Street = strings(height(T),1); end
if ~ismember('Direction',T.Properties.VariableNames), T.Direction = strings(height(T),1); end
T.EdgeIDs = string(T.EdgeIDs); T.EdgeIDs(ismissing(T.EdgeIDs)) = "";
T.Street(ismissing(T.Street)) = ""; T.Direction(ismissing(T.Direction)) = "";
key = lower(strtrim(string(cleanNames)));
en = G.Edges.EndNodes; nE = size(en,1);
fprintf('Edge overrides from %s:\n',file);
T.Action(ismissing(T.Action)) = "";
for k = 1:height(T)
    if strlength(strtrim(T.Action(k))) == 0, continue; end      % row not decided yet
    nm = lower(strtrim(string(cleanSectionNamesLocal(T.SectionName(k)))));
    rows = find(key == nm);
    if isempty(rows)
        warning('runParmaTrafficTomography:EdgeOverrideName','Override row %d: section "%s" not found.',k,T.SectionName(k));
        continue
    end
    act = lower(strtrim(T.Action(k)));
    for r = rows(:).'
        old = double(match.DirectedEdgeIDs{r}(:)).';
        switch act
            case "exclude"
                match.Accepted(r) = false; new = [];
            case "reverse"
                new = arrayfun(@(e) oppositeEdgeLocal(G,e), old);
                if any(new == 0)
                    warning('runParmaTrafficTomography:EdgeOverrideReverse', ...
                        'No opposite edge found for "%s"; sensor excluded.',T.SectionName(k));
                    match.Accepted(r) = false; new = [];
                else
                    match.DirectedEdgeIDs{r} = new; match.Accepted(r) = true;
                end
            case "edges"
                new = str2double(split(strtrim(string(T.EdgeIDs(k))))).';
                assert(~isempty(new) && all(isfinite(new)) && all(new>=1 & new<=nE & new==round(new)), ...
                    'runParmaTrafficTomography:EdgeOverrideIDs','Override row %d: invalid EdgeIDs "%s".',k,T.EdgeIDs(k));
                match.DirectedEdgeIDs{r} = new; match.Accepted(r) = true;
            case "street"
                new = streetEdgeLocal(G,double(match.Lat(r)),double(match.Lon(r)),T.Street(k),T.Direction(k));
                if isempty(new)
                    warning('runParmaTrafficTomography:EdgeOverrideStreet', ...
                        'No edge of "%s" heading %s within 200 m of "%s"; sensor excluded.',T.Street(k),T.Direction(k),T.SectionName(k));
                    match.Accepted(r) = false;
                else
                    match.DirectedEdgeIDs{r} = new; match.Accepted(r) = true;
                end
            otherwise
                error('runParmaTrafficTomography:EdgeOverrideAction', ...
                    'Override row %d: unknown action "%s" (use reverse, exclude, edges or street).',k,T.Action(k));
        end
        try, match.AssignmentSource(r) = "manual override"; catch, end
        fprintf('  %-45s %-8s %s -> %s\n',char(strtrim(string(match.SectionName(r)))),char(act), ...
            mat2str(old),mat2str(new));
    end
end
end

function r = oppositeEdgeLocal(G,e)
% Opposite direction of edge e: the reverse edge if it exists, otherwise the
% nearest edge (midpoint within 60 m) whose bearing differs by 180 deg +- 35 deg.
en = G.Edges.EndNodes; a = en(e,1); b = en(e,2);
r = find(en(:,1)==b & en(:,2)==a,1);
if ~isempty(r), return; end
lat = double(G.Nodes.Lat); lon = double(G.Nodes.Lon);
kx = cosd(mean(lat(en(e,:))));
mx = (lon(en(:,1))+lon(en(:,2)))/2; my = (lat(en(:,1))+lat(en(:,2)))/2;
d = 111320*hypot((mx-mx(e))*kx, my-my(e));                        % metres
br = atan2d(lat(en(:,2))-lat(en(:,1)), (lon(en(:,2))-lon(en(:,1)))*kx);
db = abs(mod(br - br(e) + 180, 360) - 180);                        % 0..180
cand = find(d <= 60 & db >= 145 & (1:numel(d)).' ~= e);
if isempty(cand), r = 0; return; end
[~,i] = min(d(cand)); r = cand(i);
end

function e = streetEdgeLocal(G,lat0,lon0,street,direction)
% Nearest edge (midpoint within 200 m, at least 10 m long, so that junction
% stubs are skipped) whose road name contains STREET (case-insensitive) and
% whose bearing is within 50 deg of DIRECTION.
comp = struct('N',0,'NE',45,'E',90,'SE',135,'S',180,'SW',225,'W',270,'NW',315);
d = upper(strtrim(char(direction)));
assert(isfield(comp,d),'runParmaTrafficTomography:EdgeOverrideDirection', ...
    'Direction "%s" is not one of N, NE, E, SE, S, SW, W, NW.',direction);
target = comp.(d);
en = G.Edges.EndNodes; if ~isnumeric(en), en = [findnode(G,en(:,1)) findnode(G,en(:,2))]; end
en = double(en); lat = double(G.Nodes.Lat); lon = double(G.Nodes.Lon);
nm = lower(string(G.Edges.RoadName)); nm(ismissing(nm)) = "";
kx = cosd(lat0);
mx = (lon(en(:,1))+lon(en(:,2)))/2; my = (lat(en(:,1))+lat(en(:,2)))/2;
dist = 111320*hypot((mx-lon0)*kx, my-lat0);
br = mod(atan2d((lon(en(:,2))-lon(en(:,1)))*kx, lat(en(:,2))-lat(en(:,1))),360);
db = abs(mod(br - target + 180,360) - 180);
len = 111320*hypot((lon(en(:,2))-lon(en(:,1)))*kx, lat(en(:,2))-lat(en(:,1)));
ok = contains(nm,lower(strtrim(string(street)))) & db <= 50 & dist <= 200 & len >= 10;
e = [];
if any(ok), c = find(ok); [~,i] = min(dist(c)); e = c(i); end
end

function testDates = systematicTestDates(allDates,o,rs)
% Systematic sample within each day type: every k-th date (in date order) from a
% random start, so that test days are spread evenly over the year in every day
% type. k = round(1/(1-TrainFraction)) for working-day types and
% k = round(1/WeekendTestFraction) for types made mostly of Saturdays, Sundays
% or holidays.
allDates = allDates(:);
kWork = max(2,round(1/(1-o.TrainFraction)));
kWeekend = max(2,round(1/o.WeekendTestFraction));
[code,isHoliday] = classifyTrafficDayType(allDates,o);
weekendDay = ismember(weekday(allDates),[1 7]) | isHoliday(:);
pick = false(numel(allDates),1);
for c = unique(code(:)).'
    idx = find(code == c);                       % dates of this day type, in date order
    if isempty(idx), continue; end
    k = kWork; if mean(weekendDay(idx)) > 0.5, k = kWeekend; end
    start = randi(rs,min(k,numel(idx)));
    pick(idx(start:k:end)) = true;
end
testDates = allDates(pick);
fprintf('Test days: %d of %d (systematic within day types: every %d-th working day, every %d-th weekend day).\n', ...
    nnz(pick),numel(pick),kWork,kWeekend);
end

function [A,sel] = selectPathClassesGiven(Aall,meanY,sel,muAll,o)
% Dictionary from a given class list SEL, with the same coverage completion and
% clean-up as selectPathClasses (used by ScreeningMethod 'frequency').
sel = sel(:);
if ~isfield(o,'EnsureSensorCoverage') || o.EnsureSensorCoverage
    score = (Aall.'*meanY)./max(sum(Aall,1).',1);
    for i = find(~any(Aall(:,sel)>0,2) & meanY(:)>0).'
        if any(Aall(i,sel)>0), continue; end
        cand = find(Aall(i,:)>0);
        if isempty(cand), continue; end
        [~,ord] = sortrows([muAll(cand) score(cand)],[-1 -2]);
        sel(end+1,1) = cand(ord(1)); %#ok<AGROW>
    end
end
A = double(Aall(:,sel));
keep = any(A,1); A = A(:,keep); sel = sel(keep);
[~,ia] = unique(A.','rows','stable'); A = A(:,ia); sel = sel(ia);
end
