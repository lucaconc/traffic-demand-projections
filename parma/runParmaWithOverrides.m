function out = runParmaWithOverrides(varargin)
%RUNPARMAWITHOVERRIDES Parma study with corrected sensor matching, plus diagnostics.
%
%   out = runParmaWithOverrides()
%   out = runParmaWithOverrides('OverrideCSV','/path/parma_edge_overrides.csv', ...
%                               'Archive','/path/traffic_Parma_compact.mat')
%
%   Steps
%     1. Checks: the library's own pipeline is on the path (no old copies), the
%        archive and the correction table exist.
%     2. runParmaTrafficTomography with EdgeOverrideCSV, saved to a new file.
%     3. Verification: number of manual overrides applied against the number of
%        rows with an action in the table.
%     4. Diagnostics: reliability tally, error decomposition (day-to-day floor,
%        projected CL, NNLS, class-model misfit), edge candidates of the corrected
%        sensors (listEdgeCandidates) and edge maps of the remaining groups
%        (generateSensorEdgeMaps).
%   Everything is written to OutputDir (outside the library), with a log.
%
%   Name/value (defaults from paperConfig):
%     Archive      cfg.parma.compactFile, else
%                  <paper>/matlab/parma_realdata_tomography/traffic_data_Parma/traffic_Parma_compact.mat
%     OverrideCSV  <paper>/matlab/parma_edge_overrides.csv if it exists, else the
%                  library's reference copy data/parma_edge_overrides.csv
%     OutputDir    <paper>/matlab/parma_override_run
%     UseParallel  true
%     ParmaOptions {} extra options of runParmaTrafficTomography (added after
%                  cfg.parma.options)
%
%   out: struct with the results file, the counts and the decomposition.
cfg = paperConfig();
libRoot = fileparts(fileparts(mfilename('fullpath')));
defArchive = cfg.parma.compactFile;
if isempty(defArchive)
    defArchive = fullfile(cfg.paperRoot,'matlab','parma_realdata_tomography','traffic_data_Parma', ...
        'traffic_Parma_compact.mat');
end
defOvr = fullfile(cfg.paperRoot,'matlab','parma_edge_overrides.csv');
if ~isfile(defOvr), defOvr = fullfile(libRoot,'data','parma_edge_overrides.csv'); end
p = inputParser;
addParameter(p,'Archive',defArchive); addParameter(p,'OverrideCSV',defOvr);
addParameter(p,'OutputDir',fullfile(cfg.outputRoot,'parma_override_run'));
addParameter(p,'UseParallel',true); addParameter(p,'ParmaOptions',{});
parse(p,varargin{:}); o = p.Results;
archive = char(o.Archive); ovr = char(o.OverrideCSV);
outDir = char(safeOutputPath(o.OutputDir)); outMat = fullfile(outDir,'parma_with_overrides.mat');

%% 1. Checks
setupPaperLib();
w = which('runParmaTrafficTomography');
assert(startsWith(w,libRoot), ['An old copy of runParmaTrafficTomography is on the path (%s). ' ...
    'Remove it, or remove its folder from the path.'],w);
assert(isfile(archive),'Traffic archive not found: %s',archive);
assert(isfile(ovr),'Correction table not found: %s',ovr);
if ~isfolder(outDir), mkdir(outDir); end
here = pwd; back = onCleanup(@() cd(here)); cd(outDir);       % CSV diagnostics of the run land here
diary(fullfile(outDir,'parma_override_run_log.txt')); stopLog = onCleanup(@() diary('off'));
fprintf('Library:     %s\nArchive:     %s\nCorrections: %s\nOutputs:     %s\n\n',libRoot,archive,ovr,outDir);
T = readtable(ovr,'TextType','string','Delimiter',',');
T.Action(ismissing(T.Action)) = "";
nPlanned = nnz(strlength(strtrim(T.Action)) > 0);
fprintf('Correction table: %d rows, %d with an action.\n\n',height(T),nPlanned);

%% 2. Parma with the corrections
t0 = tic;
R = runParmaTrafficTomography(archive, cfg.parma.options{:}, o.ParmaOptions{:}, ...
    'EdgeOverrideCSV',ovr, 'SaveFile',outMat, 'UseParallel',logical(o.UseParallel), 'Plot',false);
fprintf('\nParma run finished in %.1f min; results saved to %s\n\n',toc(t0)/60,outMat);

%% 3. Were the corrections applied?
nApplied = nnz(string(R.sectionMatch.AssignmentSource) == "manual override");
fprintf('Manual overrides applied: %d of %d planned\n',nApplied,nPlanned);
if nApplied < nPlanned
    warning('runParmaWithOverrides:NotApplied', ...
        'Not all corrections were applied: look for "not found" warnings in the log.');
end

%% 4a. Reliability tally
rel = R.reliableByStratum(:); why = R.reliabilityReason(~rel);
fprintf('\nReliable strata: %d of %d\n',nnz(rel),numel(rel));
pats  = ["jump-class rank","Hessian condition","excess >","iteration limit","estimator error"];
names = ["rank deficient","Hessian condition","excess misfit vs NNLS","iteration limit","estimator error"];
tally = zeros(1,numel(pats)+1);
for k = 1:numel(pats)
    tally(k) = sum(~cellfun(@isempty,regexp(cellstr(why),pats(k),'once')));
    fprintf('  %-24s %4d strata\n',names(k),tally(k));
end
tally(end) = sum(strlength(why) == 0);
fprintf('  %-24s %4d strata\n','not fitted (no data)',tally(end));

%% 4b. Error decomposition (reliable strata, >= 2 test days)
te = R.testMeasuredMeanByStratum;  tr = R.trainMeasuredMeanByStratum;
pr = R.predictedMeanByStratum;     pn = R.predictedMeanNNLSByStratum;
ok = logical(R.sensorUsedByStratum) & repmat(rel.',size(te,1),1) ...
   & repmat(R.testDaysByStratum(:).' >= 2,size(te,1),1) ...
   & isfinite(te) & isfinite(tr) & isfinite(pr) & isfinite(pn);
nr = @(a,b) sqrt(mean((a(ok)-b(ok)).^2)) / sqrt(mean(b(ok).^2));
dec = struct('floor',nr(tr,te),'pcl',nr(pr,te),'nnls',nr(pn,te),'misfitPCL',nr(pr,tr),'misfitNNLS',nr(pn,tr));
fprintf('\nError decomposition\n');
fprintf('  Day-to-day floor   (training mean -> test day mean): %.3f\n',dec.floor);
fprintf('  Projected CL       (prediction    -> test day mean): %.3f\n',dec.pcl);
fprintf('  NNLS               (prediction    -> test day mean): %.3f\n',dec.nnls);
fprintf('  Class-model misfit (prediction -> training mean):    PCL %.3f, NNLS %.3f\n',dec.misfitPCL,dec.misfitNNLS);

%% 4c. Edge candidates of the corrected sensors
candFile = fullfile(outDir,'parma_edge_candidates.csv');
listEdgeCandidates(outMat, ovr, 'OutputFile', candFile);

%% 4d. Edge maps of the remaining groups
mapDir = fullfile(outDir,'edge_maps');
try
    generateSensorEdgeMaps(outMat,'OutputDir',mapDir);
catch ME
    warning('runParmaWithOverrides:Maps','Edge maps skipped: %s',ME.message);
end
fprintf('\nDone. Outputs in %s: parma_edge_candidates.csv, edge_maps/, parma_override_run_log.txt\n',outDir);

out = struct('resultsFile',outMat,'outputDir',outDir,'overridesPlanned',nPlanned, ...
    'overridesApplied',nApplied,'reliableStrata',nnz(rel),'strata',numel(rel), ...
    'reliabilityTally',tally,'decomposition',dec,'candidatesFile',candFile,'mapsFolder',mapDir);
end
