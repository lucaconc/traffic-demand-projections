function S=collectSyntheticPaperResults(T,o,A,sensorModels)
%COLLECTSYNTHETICPAPERRESULTS Arrange synthetic trials for writeSyntheticPaperResults.
%
%   S = collectSyntheticPaperResults(results.raw,results.settings, ...
%           results.A,results.sensorSweepModels)
%
%   Works on the output of runSyntheticOSMTomographyTests or on the saved
%   synthetic_results.mat checkpoint (variable results).
% Failed fits are NaN.
ok=ismissing(T.Failure)|strlength(T.Failure)==0;
S.sampleGrid=o.ObservationCounts;
S.E_N=sweepMatrix(T,ok,"observations",S.sampleGrid,'RelativeRateError');
S.ENNLS_N=sweepMatrix(T,ok,"observations",S.sampleGrid,'NNLSRateError');
% Baseline condition N0: clean data at BaseObservations (outlier sweep, K=0).
base=T.Experiment=="outliers" & T.Level==0 & ok;
if ~any(base)
    base=T.Experiment=="observations" & T.N==o.BaseObservations & ok;
end
S.Eact_N0=T.ActiveRateError(base);
S.F0_N0=T.InactiveFlowFraction(base);
S.outlierGrid=o.OutlierCounts;
S.E_K=sweepMatrix(T,ok,"outliers",S.outlierGrid,'RelativeRateError');
S.impurityGrid=o.ImpurityFractions;
S.E_rho=sweepMatrix(T,ok,"impurity",S.impurityGrid,'RelativeRateError');
S.sensorGrid=o.SensorCounts;
S.E_m=sweepMatrix(T,ok,"sensors",S.sensorGrid,'RelativeRateError');
S.nClasses_m=cellfun(@(M)size(M.A,2),sensorModels(:)).';
if isfield(o,'TravelSpeeds') && any(T.Experiment=="traveltime")
    S.travelGrid=o.TravelSpeeds;
    S.E_tt=sweepMatrix(T,ok,"traveltime",S.travelGrid,'RelativeRateError');
    S.ENNLS_tt=sweepMatrix(T,ok,"traveltime",S.travelGrid,'NNLSRateError');
end
rk=T.RankConditionHolds(ok);
S.rankHolds=rk(isfinite(rk))>0;
S.designSeconds=T.DesignSeconds(ok);
S.fitSeconds=T.FitSeconds(ok);
if ismember('MultiStartAtBest',T.Properties.VariableNames)
    S.msAtBest=T.MultiStartAtBest(ok); S.msDistinct=T.MultiStartDistinct(ok);
    S.msMaxRelDist=T.MultiStartMaxRelDist(ok);
end
% Comparison with the other methods (columns Err_<tag>, Sec_<tag>; projected CL in RelativeRateError)
vn=T.Properties.VariableNames; tags=regexprep(vn(startsWith(vn,'Err_')),'^Err_','');
S.cmpTags=[{'PCL'},tags];
exps={'observations','outliers','impurity','sensors','traveltime'};
grids={o.ObservationCounts,o.OutlierCounts,o.ImpurityFractions,o.SensorCounts,[]};
if isfield(o,'TravelSpeeds'), grids{5}=o.TravelSpeeds; end
for t=1:numel(S.cmpTags)
    col='RelativeRateError'; sc='Sec_PCL';
    if t>1, col=['Err_' S.cmpTags{t}]; sc=['Sec_' S.cmpTags{t}]; end
    for e=1:numel(exps)
        if isempty(grids{e}) || ~any(T.Experiment==exps{e}), continue; end
        S.cmp.(exps{e}).(S.cmpTags{t})=sweepMatrix(T,ok,exps{e},grids{e},col);
    end
    v=T.(col)(base); S.cmpBase.(S.cmpTags{t})=median(v(isfinite(v)));
    if ismember(sc,vn), w=T.(sc)(base); S.cmpSec.(S.cmpTags{t})=median(w(isfinite(w))); end
end
S.m=size(A,1);
S.n=size(A,2);
h=T.NumProjections(base);
S.H=median(h(isfinite(h)));
end

function E=sweepMatrix(T,ok,experiment,levels,column)
reps=unique(T.Replicate);
E=nan(numel(reps),numel(levels));
for k=1:numel(levels)
    for r=1:numel(reps)
        use=T.Experiment==experiment & T.Level==levels(k) & T.Replicate==reps(r) & ok;
        if any(use), E(r,k)=T.(column)(find(use,1)); end
    end
end
end

