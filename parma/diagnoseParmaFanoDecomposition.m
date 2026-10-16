function D = diagnoseParmaFanoDecomposition(R,varargin)
%DIAGNOSEPARMAFANODECOMPOSITION Explain residual Poisson overdispersion.
%
%   D = diagnoseParmaFanoDecomposition(results)
%   D = diagnoseParmaFanoDecomposition('parma_tomography_results.mat')
%
% The diagnostic uses the already-cleaned native-interval observations saved
% in results.strata. It does NOT rerun OSM processing or tomography.
%
% For each stratum and measured sensor it reports the pooled Fano factor and
% exact sample-variance decompositions obtained by grouping observations by:
%   1) calendar year,
%   2) actual weekday,
%   3) calendar year x actual weekday.
%
% For a grouping G,
%   total sample variance = within-group component + between-group component,
% so, after division by the pooled mean count,
%   pooled Fano = within Fano contribution + between Fano contribution.
%
% OPTIONS
%   MinGroupSamples  5      minimum observations for subgroup Fano summaries
%   MinMeanCount     []     uses results.settings.FanoMinMeanCount or 1
%   Plot             true
%   Visible          'on'
%   SummaryCSV       'parma_fano_decomposition.csv' ('' disables)
%   SensorCSV        'parma_fano_decomposition_sensor.csv' ('' disables)
%   FigureFile       ''     optional PDF/PNG/etc. filename via exportgraphics
%   Resolution       300    raster export resolution
%
% OUTPUT D
%   summaryTable     one row per Poisson stratum
%   sensorTable      one row per sensor x stratum
%   global           aggregate medians and counts
%   figure           diagnostic figure handle (empty when Plot=false)

ip=inputParser;
ip.FunctionName=mfilename;
addParameter(ip,'MinGroupSamples',5,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=2&&x==floor(x));
addParameter(ip,'MinMeanCount',[],@(x)isempty(x)||(isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=0));
addParameter(ip,'Plot',true,@isScalarLogicalLikeLocal);
addParameter(ip,'Visible','on',@(x)ischar(x)||isstring(x));
addParameter(ip,'SummaryCSV','parma_fano_decomposition.csv',@(x)ischar(x)||isstring(x));
addParameter(ip,'SensorCSV','parma_fano_decomposition_sensor.csv',@(x)ischar(x)||isstring(x));
addParameter(ip,'FigureFile','',@(x)ischar(x)||isstring(x));
addParameter(ip,'Resolution',300,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=72);
parse(ip,varargin{:});
o=ip.Results;
% All output paths lie outside the library (see safeOutputPath).
o=guardOutputFields(o,{'SummaryCSV','SensorCSV','FigureFile'});
o.Plot=logical(o.Plot);
o.Visible=char(string(o.Visible));
o.SummaryCSV=char(string(o.SummaryCSV));
o.SensorCSV=char(string(o.SensorCSV));
o.FigureFile=char(string(o.FigureFile));

[R,resultFile]=loadResultsLocal(R);
if ~isfield(R,'strata') || isempty(R.strata)
    error('diagnoseParmaFanoDecomposition:MissingStrata', ...
        'The results structure must contain nonempty field strata.');
end
strata=R.strata(:);
S=numel(strata);
if isempty(o.MinMeanCount)
    minMean=1;
    if isfield(R,'settings') && isfield(R.settings,'FanoMinMeanCount') && ...
            isfinite(R.settings.FanoMinMeanCount)
        minMean=double(R.settings.FanoMinMeanCount);
    end
else
    minMean=double(o.MinMeanCount);
end

% v22 strata may use a different complete-sensor subset in every slice.
% Recover the global measured-sensor count from SensorMask when available.
if isfield(strata,'SensorMask') && ~isempty(strata(1).SensorMask)
    m=numel(strata(1).SensorMask);
elseif isfield(R,'selectedSections') && ~isempty(R.selectedSections)
    m=numel(R.selectedSections);
else
    m=size(strata(1).Y,2);
end
for q=1:S
    if numel(strata(q).Dates)~=size(strata(q).Y,1)
        error('diagnoseParmaFanoDecomposition:Dimensions', ...
            'Stratum %d has inconsistent Y/Dates dimensions.',q);
    end
    if isfield(strata,'SensorMask') && numel(strata(q).SensorMask)~=m
        error('diagnoseParmaFanoDecomposition:SensorMask', ...
            'Stratum %d has an inconsistent SensorMask length.',q);
    end
end

nRows=S*m;
Stratum=zeros(nRows,1); Sensor=zeros(nRows,1); UsedInEstimation=false(nRows,1); DayType=strings(nRows,1); Year=zeros(nRows,1); YearLabel=strings(nRows,1);
StartTime=strings(nRows,1); EndTime=strings(nRows,1); SeasonStart=strings(nRows,1); SeasonEnd=strings(nRows,1);
N=zeros(nRows,1); MeanCount=nan(nRows,1); PooledFano=nan(nRows,1);
YearGroups=zeros(nRows,1); MedianWithinYearFano=nan(nRows,1); YearWithinFanoContribution=nan(nRows,1); YearBetweenFanoContribution=nan(nRows,1); YearBetweenVarianceFraction=nan(nRows,1);
WeekdayGroups=zeros(nRows,1); MedianWithinWeekdayFano=nan(nRows,1); WeekdayWithinFanoContribution=nan(nRows,1); WeekdayBetweenFanoContribution=nan(nRows,1); WeekdayBetweenVarianceFraction=nan(nRows,1);
YearWeekdayGroups=zeros(nRows,1); MedianWithinYearWeekdayFano=nan(nRows,1); YearWeekdayWithinFanoContribution=nan(nRows,1); YearWeekdayBetweenFanoContribution=nan(nRows,1); YearWeekdayBetweenVarianceFraction=nan(nRows,1);

r=0;
for q=1:S
    st=strata(q);
    dates=st.Dates(:);
    yy=year(dates);
    wd=weekday(dates); % Sunday=1,...,Saturday=7 (same convention used in pipeline)
    yw=yy*10+wd;
    if isfield(st,'SensorMask') && ~isempty(st.SensorMask)
        sensorMask=logical(st.SensorMask(:));
    else
        sensorMask=true(m,1);
    end
    activeRows=find(sensorMask);
    for j=1:m
        r=r+1;
        Stratum(r)=q; Sensor(r)=j; UsedInEstimation(r)=sensorMask(j);
        DayType(r)=getStringFieldLocal(st,'DayType','');
        Year(r)=getNumericFieldLocal(st,'Year',0);
        YearLabel(r)=getStringFieldLocal(st,'YearLabel','AllYears');
        StartTime(r)=minuteLabelLocal(getNumericFieldLocal(st,'StartMinute',NaN));
        EndTime(r)=minuteLabelLocal(getNumericFieldLocal(st,'EndMinute',NaN));
        SeasonStart(r)=getStringFieldLocal(st,'SeasonStart','');
        SeasonEnd(r)=getStringFieldLocal(st,'SeasonEnd','');
        if ~sensorMask(j)
            N(r)=0;
            continue
        end
        localCol=find(activeRows==j,1);
        if isempty(localCol) || localCol>size(st.Y,2)
            N(r)=0;
            continue
        end
        y=double(st.Y(:,localCol));
        good=isfinite(y) & y>=0 & ~isnat(dates);
        y=y(good); yg=yy(good); wg=wd(good); ywg=yw(good);
        N(r)=numel(y);
        if numel(y)<2, continue; end
        mu=mean(y); MeanCount(r)=mu;
        if ~(isfinite(mu)&&mu>=minMean), continue; end
        tv=var(y,0); PooledFano(r)=tv/max(mu,eps);

        a=groupDecompLocal(y,yg,mu,o.MinGroupSamples,minMean);
        YearGroups(r)=a.nGroups; MedianWithinYearFano(r)=a.medianGroupFano;
        YearWithinFanoContribution(r)=a.withinVar/max(mu,eps);
        YearBetweenFanoContribution(r)=a.betweenVar/max(mu,eps);
        YearBetweenVarianceFraction(r)=a.betweenFraction;

        a=groupDecompLocal(y,wg,mu,o.MinGroupSamples,minMean);
        WeekdayGroups(r)=a.nGroups; MedianWithinWeekdayFano(r)=a.medianGroupFano;
        WeekdayWithinFanoContribution(r)=a.withinVar/max(mu,eps);
        WeekdayBetweenFanoContribution(r)=a.betweenVar/max(mu,eps);
        WeekdayBetweenVarianceFraction(r)=a.betweenFraction;

        a=groupDecompLocal(y,ywg,mu,o.MinGroupSamples,minMean);
        YearWeekdayGroups(r)=a.nGroups; MedianWithinYearWeekdayFano(r)=a.medianGroupFano;
        YearWeekdayWithinFanoContribution(r)=a.withinVar/max(mu,eps);
        YearWeekdayBetweenFanoContribution(r)=a.betweenVar/max(mu,eps);
        YearWeekdayBetweenVarianceFraction(r)=a.betweenFraction;
    end
end

sensorTable=table(Stratum,Sensor,UsedInEstimation,DayType,Year,YearLabel,StartTime,EndTime,SeasonStart,SeasonEnd,N,MeanCount,PooledFano, ...
    YearGroups,MedianWithinYearFano,YearWithinFanoContribution,YearBetweenFanoContribution,YearBetweenVarianceFraction, ...
    WeekdayGroups,MedianWithinWeekdayFano,WeekdayWithinFanoContribution,WeekdayBetweenFanoContribution,WeekdayBetweenVarianceFraction, ...
    YearWeekdayGroups,MedianWithinYearWeekdayFano,YearWeekdayWithinFanoContribution,YearWeekdayBetweenFanoContribution,YearWeekdayBetweenVarianceFraction);

% One-row-per-stratum summary, taking robust medians across measured sensors.
StratumS=(1:S).'; DayTypeS=strings(S,1); YearS=zeros(S,1); YearLabelS=strings(S,1); StartTimeS=strings(S,1); EndTimeS=strings(S,1); SeasonStartS=strings(S,1); SeasonEndS=strings(S,1);
Observations=zeros(S,1); SensorsUsed=zeros(S,1);
MedianPooledFano=nan(S,1); MedianWithinYearFanoS=nan(S,1); MedianWithinWeekdayFanoS=nan(S,1); MedianWithinYearWeekdayFanoS=nan(S,1);
MedianYearBetweenFraction=nan(S,1); MedianWeekdayBetweenFraction=nan(S,1); MedianYearWeekdayBetweenFraction=nan(S,1);
MedianYearBetweenFanoContribution=nan(S,1); MedianWeekdayBetweenFanoContribution=nan(S,1); MedianYearWeekdayBetweenFanoContribution=nan(S,1);
for q=1:S
    idx=find(sensorTable.Stratum==q);
    st=strata(q);
    DayTypeS(q)=getStringFieldLocal(st,'DayType','');
    YearS(q)=getNumericFieldLocal(st,'Year',0);
    YearLabelS(q)=getStringFieldLocal(st,'YearLabel','AllYears');
    StartTimeS(q)=minuteLabelLocal(getNumericFieldLocal(st,'StartMinute',NaN));
    EndTimeS(q)=minuteLabelLocal(getNumericFieldLocal(st,'EndMinute',NaN));
    SeasonStartS(q)=getStringFieldLocal(st,'SeasonStart','');
    SeasonEndS(q)=getStringFieldLocal(st,'SeasonEnd','');
    Observations(q)=size(st.Y,1);
    SensorsUsed(q)=nnz(sensorTable.UsedInEstimation(idx));
    MedianPooledFano(q)=medianFiniteLocal(sensorTable.PooledFano(idx));
    MedianWithinYearFanoS(q)=medianFiniteLocal(sensorTable.MedianWithinYearFano(idx));
    MedianWithinWeekdayFanoS(q)=medianFiniteLocal(sensorTable.MedianWithinWeekdayFano(idx));
    MedianWithinYearWeekdayFanoS(q)=medianFiniteLocal(sensorTable.MedianWithinYearWeekdayFano(idx));
    MedianYearBetweenFraction(q)=medianFiniteLocal(sensorTable.YearBetweenVarianceFraction(idx));
    MedianWeekdayBetweenFraction(q)=medianFiniteLocal(sensorTable.WeekdayBetweenVarianceFraction(idx));
    MedianYearWeekdayBetweenFraction(q)=medianFiniteLocal(sensorTable.YearWeekdayBetweenVarianceFraction(idx));
    MedianYearBetweenFanoContribution(q)=medianFiniteLocal(sensorTable.YearBetweenFanoContribution(idx));
    MedianWeekdayBetweenFanoContribution(q)=medianFiniteLocal(sensorTable.WeekdayBetweenFanoContribution(idx));
    MedianYearWeekdayBetweenFanoContribution(q)=medianFiniteLocal(sensorTable.YearWeekdayBetweenFanoContribution(idx));
end
summaryTable=table(StratumS,DayTypeS,YearS,YearLabelS,StartTimeS,EndTimeS,SeasonStartS,SeasonEndS,Observations,SensorsUsed, ...
    MedianPooledFano,MedianWithinYearFanoS,MedianWithinWeekdayFanoS,MedianWithinYearWeekdayFanoS, ...
    MedianYearBetweenFraction,MedianWeekdayBetweenFraction,MedianYearWeekdayBetweenFraction, ...
    MedianYearBetweenFanoContribution,MedianWeekdayBetweenFanoContribution,MedianYearWeekdayBetweenFanoContribution, ...
    'VariableNames',{'Stratum','DayType','Year','YearLabel','StartTime','EndTime','SeasonStart','SeasonEnd','Observations','SensorsUsed', ...
    'MedianPooledFano','MedianWithinYearFano','MedianWithinWeekdayFano','MedianWithinYearWeekdayFano', ...
    'MedianYearBetweenVarianceFraction','MedianWeekdayBetweenVarianceFraction','MedianYearWeekdayBetweenVarianceFraction', ...
    'MedianYearBetweenFanoContribution','MedianWeekdayBetweenFanoContribution','MedianYearWeekdayBetweenFanoContribution'});

G=struct();
G.medianPooledFano=medianFiniteLocal(sensorTable.PooledFano);
G.medianWithinYearFano=medianFiniteLocal(sensorTable.MedianWithinYearFano);
G.medianWithinWeekdayFano=medianFiniteLocal(sensorTable.MedianWithinWeekdayFano);
G.medianWithinYearWeekdayFano=medianFiniteLocal(sensorTable.MedianWithinYearWeekdayFano);
G.medianYearBetweenVarianceFraction=medianFiniteLocal(sensorTable.YearBetweenVarianceFraction);
G.medianWeekdayBetweenVarianceFraction=medianFiniteLocal(sensorTable.WeekdayBetweenVarianceFraction);
G.medianYearWeekdayBetweenVarianceFraction=medianFiniteLocal(sensorTable.YearWeekdayBetweenVarianceFraction);
G.nStrata=S; G.nSensors=m; G.minGroupSamples=o.MinGroupSamples; G.minMeanCount=minMean;

fprintf('\n================ FANO DECOMPOSITION ================\n');
fprintf('Median pooled Fano                         : %.4g\n',G.medianPooledFano);
fprintf('Median within-year subgroup Fano           : %.4g\n',G.medianWithinYearFano);
fprintf('Median within-weekday subgroup Fano        : %.4g\n',G.medianWithinWeekdayFano);
fprintf('Median within year x weekday subgroup Fano : %.4g\n',G.medianWithinYearWeekdayFano);
fprintf('Median variance fraction between years     : %.1f%%\n',100*G.medianYearBetweenVarianceFraction);
fprintf('Median variance fraction between weekdays  : %.1f%%\n',100*G.medianWeekdayBetweenVarianceFraction);
fprintf('Median fraction between year x weekday     : %.1f%%\n',100*G.medianYearWeekdayBetweenVarianceFraction);
fprintf('====================================================\n');

if ~isempty(o.SummaryCSV)
    writetable(summaryTable,o.SummaryCSV);
    fprintf('Fano decomposition summary written to: %s\n',o.SummaryCSV);
end
if ~isempty(o.SensorCSV)
    writetable(sensorTable,o.SensorCSV);
    fprintf('Fano decomposition sensor table written to: %s\n',o.SensorCSV);
end

f=[];
if o.Plot
    f=figure('Name','Parma Fano decomposition','Color','w','Visible',o.Visible, ...
        'Units','inches','Position',[1 1 10.5 7.0]);
    tl=tiledlayout(f,2,2,'TileSpacing','compact','Padding','compact');

    ax=nexttile(tl); hold(ax,'on');
    [~,ord]=sort(summaryTable.MedianPooledFano,'descend','MissingPlacement','last');
    x=(1:S).';
    plot(ax,x,summaryTable.MedianPooledFano(ord),'-o','MarkerSize',3,'DisplayName','Pooled');
    plot(ax,x,summaryTable.MedianWithinYearFano(ord),'-','DisplayName','Within year');
    plot(ax,x,summaryTable.MedianWithinWeekdayFano(ord),'-','DisplayName','Within weekday');
    plot(ax,x,summaryTable.MedianWithinYearWeekdayFano(ord),'-','DisplayName','Within year x weekday');
    yline(ax,1,'k--','Poisson'); grid(ax,'on'); box(ax,'on');
    xlabel(ax,'Strata sorted by pooled Fano'); ylabel(ax,'Median sensor Fano');
    title(ax,'How subdivision changes dispersion'); legend(ax,'Location','best'); paperAxesLocal(ax);

    ax=nexttile(tl);
    vals=100*[summaryTable.MedianYearBetweenVarianceFraction summaryTable.MedianWeekdayBetweenVarianceFraction];
    good=all(isfinite(vals),2);
    if any(good)
        scatter(ax,vals(good,1),vals(good,2),24,'filled'); hold(ax,'on');
        plot(ax,[0 100],[0 100],'k--'); xlim(ax,[0 100]); ylim(ax,[0 100]);
    end
    grid(ax,'on'); box(ax,'on'); xlabel(ax,'Variance between years [%]'); ylabel(ax,'Variance between weekdays [%]');
    title(ax,'Alternative sources of between-group variance'); paperAxesLocal(ax);

    ax=nexttile(tl);
    yy=100*summaryTable.MedianYearBetweenVarianceFraction; yy=yy(isfinite(yy));
    histogram(ax,yy,'BinLimits',[0 100]); xlabel(ax,'Between-year share of variance [%]'); ylabel(ax,'Number of strata');
    title(ax,'Year-to-year drift'); grid(ax,'on'); box(ax,'on'); paperAxesLocal(ax);

    ax=nexttile(tl);
    globalVals=[G.medianPooledFano G.medianWithinYearFano G.medianWithinWeekdayFano G.medianWithinYearWeekdayFano];
    bar(ax,globalVals); hold(ax,'on'); yline(ax,1,'k--','Poisson');
    xticks(ax,1:4); xticklabels(ax,{'Pooled','Within year','Within weekday','Within year x weekday'}); xtickangle(ax,18);
    ylabel(ax,'Median Fano across sensor-strata'); title(ax,'Global residual dispersion'); grid(ax,'on'); box(ax,'on'); paperAxesLocal(ax);
    title(tl,'Fano decomposition: temporal mixing versus residual overdispersion');

    if ~isempty(o.FigureFile)
        [~,~,ext]=fileparts(o.FigureFile);
        ext=lower(ext);
        if ismember(ext,{'.png','.jpg','.jpeg','.tif','.tiff'})
            exportgraphics(f,o.FigureFile,'Resolution',o.Resolution,'BackgroundColor','white');
        else
            exportgraphics(f,o.FigureFile,'ContentType','vector','BackgroundColor','white');
        end
    end
end

D=struct();
D.version='2026-09-15-fano-decomposition-v17';
D.resultFile=string(resultFile);
D.summaryTable=summaryTable;
D.sensorTable=sensorTable;
D.global=G;
D.figure=f;
end

function a=groupDecompLocal(y,g,mu,minGroupSamples,minMean)
y=double(y(:)); g=double(g(:));
good=isfinite(y) & isfinite(g); y=y(good); g=g(good);
a=struct('nGroups',0,'withinVar',NaN,'betweenVar',NaN,'betweenFraction',NaN,'medianGroupFano',NaN);
N=numel(y); if N<2, return; end
ug=unique(g,'stable'); a.nGroups=numel(ug);
withinSS=0; betweenSS=0; gf=nan(numel(ug),1);
for k=1:numel(ug)
    v=y(g==ug(k)); ng=numel(v); mg=mean(v);
    if ng>=2, withinSS=withinSS+(ng-1)*var(v,0); end
    betweenSS=betweenSS+ng*(mg-mu)^2;
    if ng>=minGroupSamples && mg>=minMean && ng>=2
        gf(k)=var(v,0)/max(mg,eps);
    end
end
a.withinVar=withinSS/(N-1);
a.betweenVar=betweenSS/(N-1);
tv=a.withinVar+a.betweenVar;
a.betweenFraction=a.betweenVar/max(tv,eps);
a.medianGroupFano=medianFiniteLocal(gf);
end

function v=medianFiniteLocal(x)
x=double(x(:)); x=x(isfinite(x));
if isempty(x), v=NaN; else, v=median(x); end
end

function s=getStringFieldLocal(st,name,default)
if isfield(st,name), s=string(st.(name)); else, s=string(default); end
end

function x=getNumericFieldLocal(st,name,default)
if isfield(st,name), x=double(st.(name)); else, x=default; end
end

function s=minuteLabelLocal(m)
if ~isfinite(m), s=""; return; end
m=round(m); hh=floor(m/60); mm=mod(m,60);
if hh>=24 && mm==0, s="24:00"; else, s=string(sprintf('%02d:%02d',mod(hh,24),mm)); end
end

function [R,resultFile]=loadResultsLocal(input)
resultFile='';
if isstruct(input), R=input; return; end
if ~(ischar(input)||isstring(input))
    error('diagnoseParmaFanoDecomposition:Input','Input must be results struct or MAT filename.');
end
resultFile=char(string(input));
if exist(resultFile,'file')~=2
    error('diagnoseParmaFanoDecomposition:MissingFile','Result MAT file not found: %s',resultFile);
end
S=load(resultFile);
if isfield(S,'results') && isstruct(S.results), R=S.results; return; end
fn=fieldnames(S); isS=false(numel(fn),1);
for k=1:numel(fn), isS(k)=isstruct(S.(fn{k})); end
idx=find(isS);
if numel(idx)==1, R=S.(fn{idx}); return; end
error('diagnoseParmaFanoDecomposition:Variable', ...
    'MAT file must contain variable "results" or exactly one compatible result struct.');
end

function paperAxesLocal(ax)
set(ax,'FontName','Helvetica','FontSize',9,'LineWidth',0.8,'TickDir','out');
end

function tf=isScalarLogicalLikeLocal(x)
tf=(islogical(x)||isnumeric(x)) && isscalar(x) && any(x==[0 1]);
end
