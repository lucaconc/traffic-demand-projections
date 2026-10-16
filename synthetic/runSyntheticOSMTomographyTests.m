function results=runSyntheticOSMTomographyTests(varargin)
%RUNSYNTHETICOSMTOMOGRAPHYTESTS Single-stratum Monte Carlo OSM benchmark.
% results=runSyntheticOSMTomographyTests;
% results=runSyntheticOSMTomographyTests('Graph',parmaResults.G);
% results=runSyntheticOSMTomographyTests('Replicates',30,'OutputDir','synthetic_tests');
%
% Dependencies: MATLAB R2020b+, Statistics and Machine Learning Toolbox
% (poissrnd/gamrnd), Optimization Toolbox (lsqnonneg), and your existing
% estimateMuOptimizedProjections.m with its dependencies on the MATLAB path.
% Uses the SAME estimator entry point/options as the Parma pipeline.
% CompareMethods ({'NNLS (means)','Cumulant matching','Bootstrap GLS','Projected GMM',
% 'Pairwise CL'}): also fitted on the same data and classes, with errors and
% times in columns Err_<tag>, Sec_<tag> (tags NNLS, Cum, GLS, GMM, Pair); cumulant
% matching and GLS use the identifying order up to MaxCumulantOrder (4), GLS uses
% BootstrapResamples (200), the pairwise likelihood stops after PairwiseMaxSeconds
% (30), the cost-aware selection uses CostExponent (0.5) and integer groups of
% at most CostMaxGroup (2) sensors (tag CostCL). Use {} to fit only the
% projected composite likelihood. NNLS is
% initialization only: the full reduced candidate dictionary is retained so
% false-positive flow on inactive classes contributes to the error metric.
% No temporal partitioning, outlier cleaning, or validation-based tuning.
%
% Default city: central Parma OSM bbox [south west north east]. Supply BBox
% for another city, or Graph to reuse the saved OSM digraph from Parma.
% Downloads/cache an OSM major-road graph; respects one-way/access tags but
% does not model turn restrictions, temporary restrictions, or traffic lights.
%
% Experiment design:
% * Sample distinct shortest paths; first NumTruePaths have positive random
%   Poisson rates. The remaining candidates have exactly zero true flow.
% * Place sensors to cover every candidate path, then group identical sensor
%   signatures. Ground truth is summed over member paths in each class.
% * Rates are random Uniform(RateRange) per TRUE PATH per native interval.
%   Rates/dictionary are held fixed across conditions within a replicate.
% * Sample-size sweep: prefixes of one clean sample, no contamination.
% * Outlier sweep: at fixed BaseObservations, add positive count spikes to
%   exactly K distinct observation rows, one randomly selected sensor each.
%   Spike=max(1,Poisson(OutlierMultiplier * that sensor's mean clean count)).
%   K means rows contaminated, not extra dictionary paths; sample N is fixed.
% * Impurity sweep: independently replace each latent path count with a
%   Gamma-Poisson count with probability rho. The mean is unchanged;
%   Var(X)=mu+rho*ImpurityCV^2*mu^2, hence Fano=1+rho*ImpurityCV^2*mu.
%   At rho=0 observations are pure Poisson; only overdispersion is tested.
% * Main error: norm(muHat-muTrue)/norm(muTrue), over ALL candidate classes.
%   Unidentifiable individual path rates are never used as scoring targets.
% * Sensor-count sweep: SensorCounts=[5 10 15 20 25 30], fixed BaseObservations,
%   clean counts, nested prefixes of the sensor selection order. Regroup
%   candidate classes for each prefix; omit zero signatures from fitting.
%   Main error compares observable merged-class rates. Record observable true
%   flow fraction and expanded full-class error (equal-share convention).
% * Same topology/sensors for all trials; uncertainty is conditional on this
%   graph/dictionary. Random rates/counts vary between independent replicates.
%   Default: 50 independent trials.
% * Curves: median and 10th--90th replicate percentiles (NOT confidence bands).
%   Failed fits are NaN and reported in CSVs; effective sample sizes are recorded in summary.csv.
%
% API references:
% https://www.mathworks.com/help/stats/gamrnd.html
% https://wiki.openstreetmap.org/wiki/Overpass_API/Overpass_QL
%
% Outputs: four violin error PDFs/PNGs/FIGs through the companion plotter,
% an OSM sensor/path map, raw_trials.csv, summary.csv, and a MAT
% checkpoint after every replicate. Figures have no titles, with black text.
% Re-run generates new results from the same seed; checkpoint is not a resume.
%
% Paper results (added): every fit also records the error on the active
% classes, the fraction of estimated flow on inactive classes, the error of
% the NNLS fit of the means (baseline), whether the jump-class rank condition
% holds, and the design/fit times.  At the end, writeSyntheticPaperResults
% writes the \ResSyn... macros of the paper to PaperResultsFile, including the
% composite-versus-exact-likelihood comparison on the running example
% (FullLikelihoodReplicates>0).  Requires syntheticPaperMetrics.m,
% writeSyntheticPaperResults.m, runningExampleFullLikelihood.m and
% writePaperMacros.m on the path.

p=inputParser;
addParameter(p,'Graph',[],@(x)isempty(x)||isa(x,'digraph'));
addParameter(p,'BBox',[44.77 10.29 44.83 10.37],@(x)isnumeric(x)&&numel(x)==4&&all(isfinite(x)));
addParameter(p,'OutputDir','synthetic_osm_tests',@(x)ischar(x)||isstring(x));
addParameter(p,'Replicates',50,@positiveInteger);
addParameter(p,'NumTruePaths',12,@positiveInteger);
addParameter(p,'NumCandidatePaths',60,@positiveInteger);
addParameter(p,'NumSensors',30,@positiveInteger);
addParameter(p,'SensorCounts',[5 10 15 20 25 30],@integerVector);
addParameter(p,'MinPathEdges',5,@positiveInteger);
addParameter(p,'PathMaxTries',30000,@positiveInteger);
addParameter(p,'RateRange',[5 30],@(x)isnumeric(x)&&numel(x)==2&&all(isfinite(x))&&x(1)>0&&x(2)>=x(1));
addParameter(p,'ObservationMinutes',15,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0);
addParameter(p,'ObservationCounts',[25 50 100 250 500 1000],@integerVector);
addParameter(p,'BaseObservations',500,@positiveInteger);
addParameter(p,'OutlierCounts',[0 2 5 10 20 40],@nonnegativeIntegerVector);
addParameter(p,'OutlierMultiplier',5,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0);
addParameter(p,'ImpurityFractions',[0 .05 .1 .2 .4 .6],@(x)isnumeric(x)&&isvector(x)&&~isempty(x)&&all(isfinite(x))&&all(x>=0&x<=1));
addParameter(p,'ImpurityCV',0.5,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0);
addParameter(p,'Seed',20260930,@positiveInteger);
addParameter(p,'EstimatorOptions',{},@iscell);
% Comparison with the methods of the benchmarks, on the same data and classes
addParameter(p,'CompareMethods',{'NNLS (means)','Cumulant matching','Bootstrap GLS','Projected GMM','Pairwise CL','Cost-aware CL','Prior GLS'},@(x)iscell(x)||isstring(x));
addParameter(p,'PriorCV',0.5,@(x)isnumeric(x)&&isscalar(x)&&x>0);   % relative error of the prior demand of 'Prior GLS'
addParameter(p,'CostExponent',0.5,@(x)isnumeric(x)&&isscalar(x)&&x>=0);   % cost-aware selection
addParameter(p,'CostMaxGroup',2,@positiveInteger);
addParameter(p,'ParallelReplicates',true,@(x)islogical(x)&&isscalar(x));   % parfor over replicates when a pool is open
addParameter(p,'MaxCumulantOrder',4,@positiveInteger);
addParameter(p,'BootstrapResamples',200,@positiveInteger);
addParameter(p,'PairwiseMaxSeconds',30,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'Visible','on',@(x)any(strcmpi(x,{'on','off'})));
addParameter(p,'PaperResultsFile','',@(x)ischar(x)||isstring(x));
addParameter(p,'PaperResults',true,@(x)islogical(x)&&isscalar(x));   % false: no paper macros
addParameter(p,'FullLikelihoodReplicates',200,@(x)isnumeric(x)&&isscalar(x)&&x>=0&&x==floor(x));
addParameter(p,'FullLikelihoodN',500,@positiveInteger);
addParameter(p,'NumStarts',3,@positiveInteger);
addParameter(p,'UseParallel',false,@(x)islogical(x)&&isscalar(x));
% Travel-time sweep: vehicles reach successive sensors after the travel time along
% the path at this speed (km/h); Inf = all sensors of a path in the same interval.
addParameter(p,'TravelSpeeds',[Inf 50 30 20 10],@(x)isnumeric(x)&&isvector(x)&&all(x>0));
parse(p,varargin{:}); o=p.Results;
% All output paths lie outside the library (see safeOutputPath).
o=guardOutputFields(o,{'OutputDir'});
if strlength(string(o.PaperResultsFile))==0
    o.PaperResultsFile=fullfile(o.OutputDir,'paper_results_synthetic.tex');
end
for field={'ObservationCounts','OutlierCounts','ImpurityFractions','SensorCounts'}
    o.(field{1})=unique(o.(field{1})(:)).';
end
assert(o.NumCandidatePaths>o.NumTruePaths,'Candidate set must be larger than true set.');
assert(all(o.SensorCounts<=o.NumSensors),'SensorCounts must not exceed NumSensors.');
assert(all(o.OutlierCounts<=o.BaseObservations),'Outlier count cannot exceed BaseObservations.');
assert(o.BBox(1)<o.BBox(3)&&o.BBox(2)<o.BBox(4),'BBox must be [south west north east].');
deps={'estimateMuOptimizedProjections','generateSyntheticTomographyFigures','poissrnd','gamrnd','binornd','lsqnonneg', ...
        'syntheticPaperMetrics','runningExampleFullLikelihood','collectSyntheticPaperResults'};
if o.PaperResults, deps=[deps,{'writeSyntheticPaperResults','writePaperMacros'}]; end   % paper text only
for name=deps
    assert(exist(name{1},'file')~=0,'Missing dependency: %s. Add the original estimator and required toolboxes to your MATLAB path.',name{1});
end
if ~isfolder(o.OutputDir), mkdir(o.OutputDir); end
oldRng=rng; cleanup=onCleanup(@()rng(oldRng)); %#ok<NASGU>
rng(o.Seed,'twister');
[G,source]=loadGraph(o);
[paths,Apath,sensorEdges]=makeDictionary(G,o);
[signature,~,pathClass]=unique(Apath.','rows'); A=double(signature.');
pClass=size(A,2); trueA=Apath(:,1:o.NumTruePaths);
% Distance (km) from the start of each true path to each sensor it crosses.
distKm=nan(size(trueA));
for j=1:o.NumTruePaths
    w=G.Edges.Weight(paths{j}); cum=[0;cumsum(w(1:end-1))];
    for i=find(trueA(:,j)).'
        pos=find(paths{j}==sensorEdges(i),1);
        if ~isempty(pos), distKm(i,j)=cum(pos); end
    end
end
gaps=[];
for j=1:o.NumTruePaths
    d=sort(distKm(isfinite(distKm(:,j)),j)); gaps=[gaps;diff(d)]; %#ok<AGROW>
end
sensorSpacingKm=median(gaps(gaps>0));
assert(all(any(A,1))&&size(unique(A.','rows'),1)==pClass,'Invalid signature dictionary.');
assert(pClass>=2,'Dictionary has fewer than two distinguishable classes; increase sensors.');
assert(all(o.SensorCounts<=size(A,1)),'Not enough selected edges for SensorCounts.');
sensorModels=cell(numel(o.SensorCounts),1);
for k=1:numel(o.SensorCounts)
    sensorModels{k}=sensorReduction(A(1:o.SensorCounts(k),:));
end
K=numel(o.ObservationCounts)+numel(o.OutlierCounts)+numel(o.ImpurityFractions)+numel(o.SensorCounts)+numel(o.TravelSpeeds);
rows=K;                                  % one replicate; replicate tables are concatenated
raw=table('Size',[rows 12],'VariableTypes', ...
    {'double','string','double','double','double','double','double','double','double','double','double','string'}, ...
    'VariableNames',{'Replicate','Experiment','Level','N','OutlierRows','ActualImpurity', ...
    'RelativeRateError','SensorMeanError','FalsePositiveMass','Seconds','FirstOrderOptimality','Failure'});
raw.RelativeRateError(:)=NaN; raw.SensorMeanError(:)=NaN; raw.FalsePositiveMass(:)=NaN;
raw.FirstOrderOptimality(:)=NaN;
raw.SensorCount=zeros(rows,1); raw.FittedClassCount=zeros(rows,1);
raw.ObservableFlowFraction=nan(rows,1); raw.ExpandedClassRateError=nan(rows,1);
% Added for the paper: detection vs estimation, NNLS baseline, rank, times.
raw.ActiveRateError=nan(rows,1);        % error on classes with positive true rate
raw.InactiveFlowFraction=nan(rows,1);   % estimated flow on inactive classes / total estimated flow
raw.NNLSRateError=nan(rows,1);          % relative error of the NNLS fit of the means
raw.RankConditionHolds=nan(rows,1);     % 1 if jump-class rank equals the number of classes
raw.NumProjections=nan(rows,1);
raw.DesignSeconds=nan(rows,1); raw.FitSeconds=nan(rows,1);
raw.MultiStartAtBest=nan(rows,1);       % share of starts reaching the best maximum
raw.MultiStartDistinct=nan(rows,1);     % distinct local maxima found
raw.MultiStartMaxRelDist=nan(rows,1);   % largest ||mu_k-muHat||/||muHat||
% Per-method errors and times (projected CL: RelativeRateError, Sec_PCL)
cmpNames=cellstr(o.CompareMethods); cmpTags=cellfun(@methodTag,cmpNames,'UniformOutput',false);
raw.Sec_PCL=nan(rows,1);
for t=1:numel(cmpTags)
    raw.(['Err_' cmpTags{t}])=nan(rows,1); raw.(['Sec_' cmpTags{t}])=nan(rows,1);
    raw.(['Def_' cmpTags{t}])=nan(rows,1);   % 1 if the method's system is rank deficient
end
% Identifying order of the cumulant system (full design and each sensor prefix)
ordFull=identifyingOrderLocal(A,o.MaxCumulantOrder);
ordSensor=cellfun(@(M)identifyingOrderLocal(M.A,o.MaxCumulantOrder),sensorModels,'UniformOutput',false);
results=struct('settings',o,'G',G,'mapSource',source,'candidatePaths',{paths}, ...
    'truePathIndices',(1:o.NumTruePaths).','sensorEdges',sensorEdges, ...
    'Apath',Apath,'A',A,'pathClass',pathClass,'completedReplicates',0);
results.requestedCandidatePaths=o.NumCandidatePaths;
results.inactiveClassCount=pClass-numel(unique(pathClass(1:o.NumTruePaths)));
results.rateUnits='vehicles per native observation interval';
results.hourlyRateFactor=60/o.ObservationMinutes;
results.metric='Relative L2 error across observable candidate classes at the tested sensor set';
results.sensorSweepModels=sensorModels;
results.sensorSweepCounts=o.SensorCounts(:);
results.sensorSweepPolicy='Nested prefixes of topological sensor order; merge duplicates and omit zero signatures';
results.travelDistanceKm=distKm;
results.sensorSpacingKm=sensorSpacingKm;
fprintf('OSM graph: %d nodes, %d edges. %d true / %d candidate paths, %d sensors, %d classes.\n', ...
    numnodes(G),numedges(G),o.NumTruePaths,numel(paths),numel(sensorEdges),pClass);
ctx=struct('o',o,'A',A,'trueA',trueA,'pClass',pClass,'pathClass',pathClass, ...
    'sensorModels',{sensorModels},'distKm',distKm,'cmpNames',{cmpNames},'cmpTags',{cmpTags}, ...
    'ordFull',ordFull,'ordSensor',{ordSensor},'rawTemplate',raw,'K',K);
usePar=o.ParallelReplicates && ~isempty(ver('parallel')) && ~isempty(gcp('nocreate')) && o.Replicates>1;
out=cell(o.Replicates,1);
if usePar
    fprintf('Running %d replicates in parallel.\n',o.Replicates);
    parfor rep=1:o.Replicates
        out{rep}=syntheticReplicate(rep,ctx);
    end
    results=assembleReplicates(results,out,o.Replicates);
    save(fullfile(o.OutputDir,'synthetic_results.mat'),'results','-v7.3');
else
    for rep=1:o.Replicates
        out{rep}=syntheticReplicate(rep,ctx);
        results=assembleReplicates(results,out,rep);              % checkpoint after every replicate
        save(fullfile(o.OutputDir,'synthetic_results.mat'),'results','-v7.3');
    end
end
raw=results.raw; row=height(raw);
results.summary=summarizeTrials(results.raw);
results.paperResults=collectSyntheticPaperResults(results.raw,o,A,sensorModels);
results.paperResults.sensorSpacingKm=sensorSpacingKm;
if exist('machineSummary','file')==2, results.machine=machineSummary(); end
if o.FullLikelihoodReplicates>0
    fprintf('Composite versus exact likelihood on the running example (%d replicates)...\n', ...
        o.FullLikelihoodReplicates);
    rng(o.Seed+999,'twister');
    results.paperResults.fullLik=runningExampleFullLikelihood('N',o.FullLikelihoodN, ...
        'R',o.FullLikelihoodReplicates,'Seed',o.Seed+999);
end
% A failure here must not lose the simulation: save and plot first if needed.
results.paperMacros=struct();
if o.PaperResults
    try
        results.paperMacros=writeSyntheticPaperResults(results.paperResults,char(o.PaperResultsFile));
        fprintf('Paper result macros written to %s\n',char(o.PaperResultsFile));
    catch ME
        warning('SyntheticOSM:PaperMacros', ...
            ['Could not write the paper macros (%s). The simulation is saved; ' ...
             'rerun writeSyntheticPaperResults(results.paperResults,file) later.'],ME.message);
    end
end
writetable(results.raw,fullfile(o.OutputDir,'raw_trials.csv'));
writetable(results.summary,fullfile(o.OutputDir,'summary.csv'));
% Save before plotting, so that a figure problem cannot cost the simulation.
save(fullfile(o.OutputDir,'synthetic_results.mat'),'results','-v7.3');
try
    results.figures=plotResults(results,o);
catch MEplot
    warning('SyntheticOSM:Figures', ['Figures not produced (%s). The results are saved; ' ...
        'rerun generateSyntheticTomographyFigures on synthetic_results.mat.'],MEplot.message);
    results.figures=gobjects(0);
end
% Do not save graphics handles in the MAT checkpoint.
figures=results.figures; results=rmfield(results,'figures');
save(fullfile(o.OutputDir,'synthetic_results.mat'),'results','-v7.3');
results.figures=figures;

fprintf('Finished. Outputs: %s\n',char(o.OutputDir));
end

function Y=travelTimeCountsLocal(X,trueA,distKm,speedKmh,dt,N)
% Counts when vehicles reach successive sensors after the travel time along
% their path. Vehicles of interval r enter the path at uniform times within it;
% at a sensor with offset delta=k+f intervals, a Binomial(X,f) share is counted
% in interval r+k+1 and the rest in r+k. Earlier intervals feed the first rows,
% so warm-up rows are drawn from X and discarded.
m=size(trueA,1);
delta=60*distKm/speedKmh/dt;                 % offsets in intervals (NaN if not crossed)
K=ceil(max(delta(:),[],'omitnan'))+1;
R=N+K;
assert(size(X,1)>=R,'Not enough clean rows for the travel-time warm-up (need %d).',R);
Yf=zeros(R+K+1,m);
for j=1:size(trueA,2)
    x=X(1:R,j);
    for i=find(trueA(:,j)).'
        d=delta(i,j); if ~isfinite(d), d=0; end
        k=floor(d); f=d-k;
        b=binornd(x,f);
        idx=(1:R).'+k;
        Yf(idx,i)=Yf(idx,i)+(x-b);
        Yf(idx+1,i)=Yf(idx+1,i)+b;
    end
end
Y=Yf(K+1:K+N,:);
end

function [mu,info]=fitRates(A,Y,o,rep)
pilot=max(lsqnonneg(A,mean(Y,1).',optimset('Display','off')),0.05);
base={'MaxCombinationSize',min(3,size(A,1)),'MaxCandidates',max(600,size(A,1)), ...
    'NumProjections',[],'RequireGCDOne',true,'TrimToIdentifiable',true, ...
    'ComputeCandidateFisherRanks',false,'MaxIter',60,'MaxSubproblemIter',40, ...
    'TolGradient',1e-5,'CandidateSeed',o.Seed+1000*rep+777,'Display','off', ...
    'NumStarts',o.NumStarts,'MultiStartSeed',o.Seed+1000*rep+555,'UseParallel',o.UseParallel};
opts=mergeNameValueOptions(base,o.EstimatorOptions);
% Never let caller-supplied truth-dependent pilots or design-only settings
% bypass this experiment's data-only initialization and fitting contract.
opts=mergeNameValueOptions(opts,{'PilotMu',pilot,'Mu0',pilot,'DesignOnly',false});
[mu,info]=estimateMuOptimizedProjections(A,Y,opts{:}); mu=double(mu(:));
assert(numel(mu)==size(A,2)&&all(isfinite(mu)&mu>=0), ...
    'Estimator returned invalid or differently indexed rates.');
end

function [G,source]=loadGraph(o)
if ~isempty(o.Graph)
    G=o.Graph; source='User-supplied OSM digraph';
else
    b=o.BBox; cache=fullfile(o.OutputDir,sprintf('osm_%.4f_%.4f_%.4f_%.4f.mat',b));
    if isfile(cache)
        d=load(cache,'G','source'); G=d.G; source=d.source;
    else
        query=sprintf(['[out:json][timeout:120];(way["highway"~"^(motorway|trunk|primary|secondary|tertiary|unclassified|residential|motorway_link|trunk_link|primary_link|secondary_link|tertiary_link)$"]' ...
            '(%.6f,%.6f,%.6f,%.6f););(._;>;);out body;'],b);
        endpoints={'https://overpass-api.de/api/interpreter','https://overpass.kumi.systems/api/interpreter'};
        osm=[]; messages=strings(0,1);
        for j=1:numel(endpoints)
            try
                candidate=webread(endpoints{j},'data',query,weboptions('Timeout',180,'ContentType','json'));
                assert(isstruct(candidate)&&isfield(candidate,'elements'),'Invalid OSM response.');
                if isfield(candidate,'remark')&&contains(lower(string(candidate.remark)),["error","timeout"])
                    error('SyntheticOSM:PartialResponse','%s',char(string(candidate.remark)));
                end
                osm=candidate; source=['OpenStreetMap via ' endpoints{j}]; break
            catch ME, messages(end+1)=string(ME.message); end %#ok<AGROW>
        end
        assert(~isempty(osm),'OSM download failed. Retry or supply Graph=parmaResults.G. Details: %s',strjoin(messages,' | '));
        G=buildDirectedRoadGraph(osm); save(cache,'G','source','osm','-v7.3');
    end
end
assert(all(ismember({'Lat','Lon'},G.Nodes.Properties.VariableNames)),'Graph nodes require Lat and Lon.');
assert(ismember('Weight',G.Edges.Properties.VariableNames)&&all(G.Edges.Weight>0),'Graph requires positive edge weights.');
% Work inside the largest strongly connected component for reachable OD pairs.
component=conncomp(G,'Type','strong'); counts=accumarray(component(:),1);
[~,big]=max(counts); G=subgraph(G,find(component==big));
assert(numnodes(G)>10,'Largest strongly connected component is too small.');
end

function [paths,Apath,sensors]=makeDictionary(G,o)
paths=cell(o.NumCandidatePaths,1); found=0;
seen=containers.Map('KeyType','char','ValueType','logical');
terminals=find(indegree(G)+outdegree(G)>=3);
if numel(terminals)<10, terminals=(1:numnodes(G)).'; end
for attempt=1:o.PathMaxTries
    od=terminals(randperm(numel(terminals),2));
    [~,distance,edges]=shortestpath(G,od(1),od(2),'Method','positive');
    if ~isfinite(distance)||numel(edges)<o.MinPathEdges, continue; end
    key=sprintf('%d,',edges);
    if isKey(seen,key), continue; end
    seen(key)=true; found=found+1; paths{found}=double(edges(:));
    if found==o.NumCandidatePaths, break; end
end
assert(found==o.NumCandidatePaths,'Not enough paths; enlarge BBox or reduce NumCandidatePaths/MinPathEdges.');
% Random order defines the true subset independently of rate/count realizations.
paths=paths(randperm(numel(paths)));
incidence=sparse(numedges(G),numel(paths));
for j=1:numel(paths), incidence(paths{j},j)=1; end
% Greedy topological coverage only, with no knowledge of which paths are active.
uncovered=true(1,numel(paths)); sensors=zeros(0,1);
while any(uncovered)
    score=full(sum(incidence(:,uncovered),2)); score(sensors)=-Inf;
    [best,e]=max(score); assert(best>0,'Cannot cover candidate paths.');
    sensors(end+1,1)=e; uncovered(logical(full(incidence(e,:))))=false; %#ok<AGROW>
end
assert(numel(sensors)<=o.NumSensors,'Increase NumSensors to at least %d to observe every candidate path.',numel(sensors));
% Prefer extra measurements that divide the current signature groups.
available=find(any(incidence,2));
while numel(sensors)<min(o.NumSensors,numel(available))
    available=setdiff(available,sensors,'stable');
    current=full(incidence(sensors,:)); [~,~,group]=unique(current.','rows');
    gain=zeros(numel(available),1);
    for g=1:max(group)
        ix=find(group==g); onesCount=full(sum(incidence(available,ix),2));
        gain=gain+onesCount.*(numel(ix)-onesCount);
    end
    [~,j]=max(gain); sensors(end+1,1)=available(j); %#ok<AGROW>
end
Apath=full(incidence(sensors,:));
assert(all(any(Apath,1)),'Unobserved candidate path.');
end

function S=summarizeTrials(T)
S=table('Size',[0 8],'VariableTypes',{'string','double','double','double','double','double','double','double'}, ...
    'VariableNames',{'Experiment','Level','MedianError','P10','P90','Successful','Failed','MedianOptimality'});
for name=["observations","outliers","impurity","sensors"]
    levels=unique(T.Level(T.Experiment==name));
    for x=levels.'
        rows=T.Experiment==name&T.Level==x; v=T.RelativeRateError(rows); good=isfinite(v);
        if any(good), q=prctile(v(good),[10 50 90]); else, q=[NaN NaN NaN]; end
        opt=T.FirstOrderOptimality(rows); opt=opt(isfinite(opt));
        if isempty(opt), optimality=NaN; else, optimality=median(opt); end
        S(end+1,:)={name,x,q(2),q(1),q(3),nnz(good),nnz(~good),optimality}; %#ok<AGROW>
    end
end
end

function figs=plotResults(R,o)
plots=generateSyntheticTomographyFigures(R,'OutputDir', ...
    fullfile(fileparts(char(o.OutputDir)),'synthetic_paper_figures'),'Visible','off');
figs=[plots.figures;gobjects(1)];
figs(5)=figure('Color','w','Visible','off','WindowStyle','normal','Position',[50 50 1000 800]);
ax=axes(figs(5)); hold(ax,'on'); G=R.G;
e=G.Edges.EndNodes; if isnumeric(e), s=e(:,1); t=e(:,2); else, s=findnode(G,e(:,1));t=findnode(G,e(:,2));end
lon=G.Nodes.Lon; lat=G.Nodes.Lat;
x=[lon(s) lon(t) nan(numel(s),1)].'; y=[lat(s) lat(t) nan(numel(s),1)].';
plot(ax,x(:),y(:),'Color',[.75 .75 .75],'HandleVisibility','off');
for j=[o.NumTruePaths+1:numel(R.candidatePaths),1:o.NumTruePaths]
    ids=R.candidatePaths{j}; color=[.65 .78 .88]; width=.8;
    if j<=o.NumTruePaths, color=[.9 .3 .08]; width=1.8; end
    xx=[lon(s(ids)) lon(t(ids)) nan(numel(ids),1)].'; yy=[lat(s(ids)) lat(t(ids)) nan(numel(ids),1)].';
    plot(ax,xx(:),yy(:),'Color',color,'LineWidth',width,'HandleVisibility','off');
end
ids=R.sensorEdges;
scatter(ax,(lon(s(ids))+lon(t(ids)))/2,(lat(s(ids))+lat(t(ids)))/2,35,'k','filled');
xlabel(ax,'Longitude [degrees]'); ylabel(ax,'Latitude [degrees]');
daspect(ax,[1/cosd(mean(lat)) 1 1]); axis(ax,'tight'); styleAxes(ax);
exportPlot(figs(5),ax,fullfile(o.OutputDir,'05_osm_paths_sensors'));
for j=1:numel(figs), set(figs(j),'Visible',o.Visible); end
end

function styleAxes(ax)
set(ax,'FontSize',15,'XColor','k','YColor','k','Color','w');
set([ax.XLabel ax.YLabel],'FontSize',17,'Color','k');
end

function exportPlot(f,target,base)
drawnow nocallbacks;
try
    exportgraphics(target,[base '.png'],'Resolution',300,'BackgroundColor','white');
    exportgraphics(target,[base '.pdf'],'ContentType','vector','BackgroundColor','white');
catch ME
    warning('SyntheticOSM:ExportFallback','Using print: %s',ME.message);
    set(f,'Units','inches'); pos=f.Position;
    set(f,'PaperUnits','inches','PaperPosition',[0 0 pos(3:4)],'PaperSize',pos(3:4),'InvertHardcopy','off');
    print(f,[base '.png'],'-dpng','-r300'); print(f,[base '.pdf'],'-dpdf','-painters');
end
savefig(f,[base '.fig']);
end

function tf=positiveInteger(x)
tf=isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=1&&x==fix(x);
end
function tf=integerVector(x)
tf=isnumeric(x)&&isvector(x)&&~isempty(x)&&all(isfinite(x))&&all(x>=1&x==fix(x));
end
function tf=nonnegativeIntegerVector(x)
tf=isnumeric(x)&&isvector(x)&&~isempty(x)&&all(isfinite(x))&&all(x>=0&x==fix(x));
end


% OSM graph utilities and option merger adapted from the supplied Parma code.
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


function model=sensorReduction(A)
p=size(A,2); observed=any(A,1).'; group=zeros(p,1);
if any(observed)
    [sig,~,g]=unique(A(:,observed).','rows'); B=sig.';
    sizes=accumarray(g,1); group(observed)=g;
    expansion=sparse(find(observed),g,1./sizes(g),p,size(B,2));
else
    B=zeros(size(A,1),0); expansion=sparse(p,0);
end
assert(norm(A*expansion-B,'fro')<1e-9,'Invalid sensor reduction.');
model=struct('A',B,'Observed',observed,'GroupIndex',group,'Expansion',expansion);
end

function tag=methodTag(name)
% Short field-safe tag of a method name.
switch char(name)
    case 'NNLS (means)', tag='NNLS';
    case 'Cumulant matching', tag='Cum';
    case 'Bootstrap GLS', tag='GLS';
    case 'Projected GMM', tag='GMM';
    case 'Pairwise CL', tag='Pair';
    case 'Cost-aware CL', tag='CostCL';
    case 'Prior GLS', tag='PriorGLS';
    otherwise, tag=matlab.lang.makeValidName(char(name));
end
end

function ord=identifyingOrderLocal(A,maxOrder)
% Smallest r <= maxOrder whose cumulant system identifies all classes ([] if none).
ord=[]; n=size(A,2);
for r=1:min(maxOrder,size(A,1))
    if rank(cumulantSystem(A,r))==n, ord=r; return; end
end
end

function C=compareMethodsLocal(A,Y,truth,ord,names,o,rep)
% Fit the comparison methods on the same data; errors are relative L2 errors
% on the fitted classes, as for the projected composite likelihood.
k=numel(names); C=struct('err',nan(1,k),'sec',nan(1,k),'def',nan(1,k));
n=size(A,2); nt=max(norm(truth),eps);
if isempty(ord), ordUse=min(o.MaxCumulantOrder,size(A,1)); else, ordUse=ord; end
rsB=RandStream('mt19937ar','Seed',o.Seed+1000*rep+999);
try, fitCumulantNNLS(A,Y,1); catch, end          % untimed warm-up
for t=1:k
    t0=tic;
    try
        switch names{t}
            case 'NNLS (means)', est=fitCumulantNNLS(A,Y,1); C.def(t)=rank(A)<n;
            case 'Cumulant matching', est=fitCumulantNNLS(A,Y,ordUse); C.def(t)=isempty(ord);
            case 'Bootstrap GLS', est=fitBootstrapGLS(A,Y,ordUse,o.BootstrapResamples,rsB); C.def(t)=isempty(ord);
            case 'Projected GMM', est=estimateMuEfficientGMM(A,Y); C.def(t)=isempty(ord);
            case 'Pairwise CL'
                est=fitPairwiseCompositeLikelihood(A,Y,fitCumulantNNLS(A,Y,min(2,size(A,1))), ...
                    'MaxSeconds',o.PairwiseMaxSeconds);
                C.def(t)=rank(cumulantSystem(A,min(2,size(A,1))))<n;
            case 'Cost-aware CL'
                pilot=max(fitCumulantNNLS(A,Y,min(size(A,1),4)),1e-3);
                Csel=selectProjectionsByCost(A,Y,pilot,'CostExponent',o.CostExponent,'MaxGroup',o.CostMaxGroup);
                est=estimateMuOptimizedProjections(A,Y,'Display','off','ProjectionMatrix',Csel,'Mu0',pilot);
                C.def(t)=0;
            case 'Prior GLS'
                % Prior demand: the true rates with a lognormal relative error of size PriorCV
                % (an outdated or survey-based matrix); inactive classes get a small positive
                % prior, a tenth of the mean active rate, so the prior does not reveal them.
                rsP=RandStream('mt19937ar','Seed',o.Seed+1000*rep+777);
                act=truth>0; base=truth; if any(act), base(~act)=0.1*mean(truth(act)); else, base(:)=1; end
                prior=base.*exp(o.PriorCV*randn(rsP,n,1)-o.PriorCV^2/2);
                est=fitPriorGLS(A,Y,prior,o.PriorCV); C.def(t)=0;
            otherwise, error('SyntheticOSM:UnknownMethod','Unknown method %s.',names{t});
        end
        C.sec(t)=toc(t0);
        est=double(est(:));
        if numel(est)==n && all(isfinite(est)), C.err(t)=norm(est-truth)/nt; end
    catch ME
        C.sec(t)=toc(t0);
        warning('SyntheticOSM:CompareFailed','%s: %s',names{t},ME.message);
    end
end
end

function R=syntheticReplicate(rep,ctx)
% One replicate of the five experiments (all levels); returns its rows of the
% trial table, its estimates and its true rates. Data and fits use seeds that
% depend only on rep, so results do not depend on the order of execution.
o=ctx.o; A=ctx.A; trueA=ctx.trueA; pClass=ctx.pClass; pathClass=ctx.pathClass;
sensorModels=ctx.sensorModels; distKm=ctx.distKm; cmpNames=ctx.cmpNames; cmpTags=ctx.cmpTags;
ordFull=ctx.ordFull; ordSensor=ctx.ordSensor; raw=ctx.rawTemplate; K=ctx.K;
muHat=nan(pClass,K); sensorEstimates=cell(numel(o.SensorCounts),1); row=0;
    rng(o.Seed+1000*rep,'twister');
    muPath=o.RateRange(1)+diff(o.RateRange)*rand(o.NumTruePaths,1);
    truth=accumarray(pathClass(1:o.NumTruePaths),muPath,[pClass 1]);
    Nmax=max([o.ObservationCounts(:);o.BaseObservations]);
    cleanX=poissrnd(repmat(muPath.',Nmax,1)); cleanY=cleanX*trueA.';
    N=o.BaseObservations; meanY=trueA*muPath;
    order=randperm(N); sensors=randi(size(A,1),N,1);
    spike=max(1,poissrnd(o.OutlierMultiplier*meanY(sensors)));
    u=rand(N,o.NumTruePaths);
    shape=1/o.ImpurityCV^2;
    multiplier=gamrnd(shape,1/shape,N,o.NumTruePaths);
    impureX=poissrnd(multiplier.*muPath.');
    cachedBaseline=[]; cachedBaselineInfo=[]; cachedCompare=[];
    for experiment=1:5
        if experiment==1, levels=o.ObservationCounts; label="observations";
        elseif experiment==2, levels=o.OutlierCounts; label="outliers";
        elseif experiment==3, levels=o.ImpurityFractions; label="impurity";
        elseif experiment==4, levels=o.SensorCounts; label="sensors";
        else, levels=o.TravelSpeeds; label="traveltime"; end
        for level=levels(:).'
            row=row+1; kOut=0; actual=0;
            fitA=A; fitTruth=truth; expansion=speye(pClass); observable=true(pClass,1);
            if experiment==1
                Y=cleanY(1:level,:);
            elseif experiment==2
                Y=cleanY(1:N,:); kOut=level; rr=order(1:kOut);
                for j=rr, Y(j,sensors(j))=Y(j,sensors(j))+spike(j); end
            elseif experiment==3
                X=cleanX(1:N,:); replace=u<level; X(replace)=impureX(replace);
                Y=X*trueA.'; actual=mean(replace(:));
            elseif experiment==5
                if isinf(level)
                    Y=cleanY(1:N,:);                  % same as the clean baseline
                else
                    rng(o.Seed+1000*rep+333+find(o.TravelSpeeds==level,1),'twister');
                    Y=travelTimeCountsLocal(cleanX,trueA,distKm,level,o.ObservationMinutes,N);
                end
            else
                modelIndex=find(o.SensorCounts==level,1);
                model=sensorModels{modelIndex};
                fitA=model.A; expansion=model.Expansion; observable=model.Observed;
                if any(observable)
                    fitTruth=accumarray(model.GroupIndex(observable),truth(observable),[size(fitA,2) 1]);
                else, fitTruth=zeros(0,1); end
                Y=cleanY(1:N,1:level);
            end
            raw.Replicate(row)=rep; raw.Experiment(row)=label; raw.Level(row)=level;
            raw.SensorCount(row)=size(fitA,1); raw.FittedClassCount(row)=size(fitA,2);
            raw.ObservableFlowFraction(row)=sum(truth(observable))/sum(truth);
            raw.N(row)=size(Y,1); raw.OutlierRows(row)=kOut; raw.ActualImpurity(row)=actual;
            % Reset fitting RNG separately: estimator randomness must not alter
            % subsequent synthetic observations or differ across paired levels.
            rng(o.Seed+1000*rep+777,'twister');
            t0=tic;
            try
                assert(~isempty(fitA)&&norm(fitTruth)>0,'No observable positive true flow for this sensor set.');
                baseline=isequal(fitA,A)&&isequal(Y,cleanY(1:N,:));
                if baseline && ~isempty(cachedBaseline)
                    estimate=cachedBaseline; info=cachedBaselineInfo;
                else
                    tP=tic; [estimate,info]=fitRates(fitA,Y,o,rep); info.pclSeconds=toc(tP);
                    if baseline, cachedBaseline=estimate; cachedBaselineInfo=info; end
                end
                if isfield(info,'pclSeconds'), raw.Sec_PCL(row)=info.pclSeconds; end
                expanded=expansion*estimate; muHat(:,row)=expanded;
                raw.RelativeRateError(row)=norm(estimate-fitTruth)/norm(fitTruth);
                raw.SensorMeanError(row)=norm(fitA*(estimate-fitTruth))/max(norm(fitA*fitTruth),eps);
                raw.FalsePositiveMass(row)=sum(estimate(fitTruth==0))/sum(fitTruth);
                raw.ExpandedClassRateError(row)=norm(expanded-truth)/norm(truth);
                if experiment==4, sensorEstimates{modelIndex}=estimate; end
                if isfield(info,'firstOrderOptimality'), raw.FirstOrderOptimality(row)=info.firstOrderOptimality; end
                M=syntheticPaperMetrics(estimate,fitTruth,fitA,Y);
                raw.ActiveRateError(row)=M.Eact;
                raw.InactiveFlowFraction(row)=M.F0;
                raw.NNLSRateError(row)=M.ENNLS;
                if isfield(info,'projectionDesignRank')
                    raw.RankConditionHolds(row)=double(info.projectionDesignRank>=size(fitA,2));
                end
                if isfield(info,'nSelectedProjections'), raw.NumProjections(row)=info.nSelectedProjections; end
                if isfield(info,'designRuntimeSeconds'), raw.DesignSeconds(row)=info.designRuntimeSeconds; end
                if isfield(info,'multistart')
                    raw.MultiStartAtBest(row)=info.multistart.fractionAtBest;
                    raw.MultiStartDistinct(row)=info.multistart.nDistinctMaxima;
                    raw.MultiStartMaxRelDist(row)=info.multistart.maxRelativeDistance;
                end
                if isfield(info,'runtimeSeconds'), raw.FitSeconds(row)=info.runtimeSeconds; end
                % Other methods on the same data, sensor set and classes
                if ~isempty(cmpNames)
                    if experiment==4, ordHere=ordSensor{modelIndex}; else, ordHere=ordFull; end
                    if baseline && ~isempty(cachedCompare)
                        C=cachedCompare;
                    else
                        rng(o.Seed+1000*rep+888,'twister');
                        C=compareMethodsLocal(fitA,Y,fitTruth,ordHere,cmpNames,o,rep);
                        if baseline, cachedCompare=C; end
                    end
                    for t=1:numel(cmpTags)
                        raw.(['Err_' cmpTags{t}])(row)=C.err(t);
                        raw.(['Sec_' cmpTags{t}])(row)=C.sec(t);
                        raw.(['Def_' cmpTags{t}])(row)=C.def(t);
                    end
                end
            catch ME
                raw.Failure(row)=string(ME.identifier)+": "+string(ME.message);
                warning('SyntheticOSM:FitFailed','Trial %d (%s=%g): %s',rep,label,level,ME.message);
            end
            raw.Seconds(row)=toc(t0);
            fprintf('Replicate %d/%d | %s=%g | relative error %.4g\n', ...
                rep,o.Replicates,label,level,raw.RelativeRateError(row));
        end
    end
R=struct('raw',raw(1:row,:),'muHat',muHat(:,1:row),'sensorEstimates',{sensorEstimates}, ...
    'truth',truth,'muPath',muPath);
end

function results=assembleReplicates(results,out,nDone)
% Concatenate the first nDone replicates into the results structure.
parts=out(1:nDone);
results.raw=vertcat(cellfun(@(r) r.raw,parts,'UniformOutput',false));
results.raw=vertcat(results.raw{:});
results.trueClassRates=cell2mat(cellfun(@(r) r.truth(:),parts(:).','UniformOutput',false));
results.truePathRates=cell2mat(cellfun(@(r) r.muPath(:),parts(:).','UniformOutput',false));
results.sensorSweepEstimates=[parts{1}.sensorEstimates];
for k=2:nDone, results.sensorSweepEstimates=[results.sensorSweepEstimates, parts{k}.sensorEstimates]; end
results.estimatedClassRates=cell2mat(cellfun(@(r) r.muHat,parts(:).','UniformOutput',false));
results.completedReplicates=nDone;
end
