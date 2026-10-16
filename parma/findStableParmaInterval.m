function D = findStableParmaInterval(input,varargin)
%FINDSTABLEPARMAINTERVAL Find a high-quality stable interval of >= 90 days.
%
%   D = findStableParmaInterval(results)
%   D = findStableParmaInterval('parma_tomography_results.mat')
%   D = findStableParmaInterval(struct('dates',dates,'Y',Y))
%   D = findStableParmaInterval(dates,Y,'Name',value,...)
%
% The diagnostic is intended for choosing a calibration/estimation period.
% Candidate intervals are constrained to a single calendar year by default.
% Quality combines:
%   * observed sensor-day coverage;
%   * lower-quartile per-sensor coverage;
%   * robust variability after removing sensor x weekday baselines;
%   * robust long-term drift after the same weekday normalization;
%   * a small preference for intervals longer than MinDays.
%
% If a RESULTS structure contains sectionEstimationMask, only estimation
% sections are used by default so held-out section-validation data cannot
% influence interval selection.
%
% OPTIONS
%   MinDays                 90
%   CandidateLengthStepDays 7
%   MinCoverage             0.90
%   MinSensorCoverage       0.80  (25th percentile of sensor coverage)
%   MaxRelativeDeviation    0.30  robust weekday-adjusted relative spread
%   MaxTrendFraction        0.20  maximum robust end-to-end drift
%   MinPositiveFraction     0.70
%   UseEstimationSensorsOnly true
%   SingleYear              true
%   TopCandidates           250
%   CSV                     'parma_stable_interval_candidates.csv'
%   Plot                    true
%   Visible                 'on'
%
% OUTPUT D fields include startDate, endDate, days, score, coverage,
% sensorCoverage25, relativeDeviation, trendFraction, passesThresholds,
% rankingTable and figure.

% Resolve input before parsing name/value options.
if isdatetime(input)
    if isempty(varargin) || ~isnumeric(varargin{1})
        error('findStableParmaInterval:Input','When the first input is datetime, the second input must be the data matrix Y.');
    end
    dates=input(:);
    Y=double(varargin{1});
    varargin=varargin(2:end);
    sourceLabel='dates/Y';
else
    [dates,Y,sourceLabel,defaultMask]=resolveStableInputLocal(input);
end

ip=inputParser;
ip.FunctionName=mfilename;
addParameter(ip,'MinDays',90,@isPositiveIntegerLocal);
addParameter(ip,'CandidateLengthStepDays',7,@isPositiveIntegerLocal);
addParameter(ip,'MinCoverage',0.90,@isUnitIntervalLocal);
addParameter(ip,'MinSensorCoverage',0.80,@isUnitIntervalLocal);
addParameter(ip,'MaxRelativeDeviation',0.30,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=0);
addParameter(ip,'MaxTrendFraction',0.20,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=0);
addParameter(ip,'MinPositiveFraction',0.70,@isUnitIntervalLocal);
addParameter(ip,'UseEstimationSensorsOnly',true,@isScalarLogicalLikeLocal);
addParameter(ip,'SingleYear',true,@isScalarLogicalLikeLocal);
addParameter(ip,'TopCandidates',250,@isPositiveIntegerLocal);
addParameter(ip,'CSV','parma_stable_interval_candidates.csv',@(x)ischar(x)||isstring(x));
addParameter(ip,'Plot',true,@isScalarLogicalLikeLocal);
addParameter(ip,'Visible','on',@(x)ischar(x)||isstring(x));
parse(ip,varargin{:});
o=ip.Results;
% All output paths lie outside the library (see safeOutputPath).
o=guardOutputFields(o,{'CSV'});
o.UseEstimationSensorsOnly=logical(o.UseEstimationSensorsOnly);
o.SingleYear=logical(o.SingleYear);
o.Plot=logical(o.Plot);
o.CSV=char(string(o.CSV));
o.Visible=char(string(o.Visible));

if exist('defaultMask','var') && o.UseEstimationSensorsOnly && ~isempty(defaultMask)
    defaultMask=logical(defaultMask(:));
    if numel(defaultMask)==size(Y,2)
        Y=Y(:,defaultMask);
    end
end

if size(Y,1)~=numel(dates)
    error('findStableParmaInterval:Size','Y must have one row for every date.');
end
if isempty(Y) || size(Y,2)<1
    error('findStableParmaInterval:Data','No sensor columns are available for interval selection.');
end

% Sort chronologically and collapse duplicate dates conservatively by taking the
% mean of finite values for each sensor on that date.
[dates,ord]=sort(dates(:));
Y=Y(ord,:);
[uniqueDates,~,g]=unique(dates,'stable');
if numel(uniqueDates)<numel(dates)
    Yu=nan(numel(uniqueDates),size(Y,2));
    for k=1:numel(uniqueDates)
        Yu(k,:)=mean(Y(g==k,:),1,'omitnan');
    end
    dates=uniqueDates; Y=Yu;
end
Y(~isfinite(Y) | Y<0)=NaN;

if numel(dates)<o.MinDays
    error('findStableParmaInterval:TooShort','Only %d dates are available; at least %d are required.',numel(dates),o.MinDays);
end

% Candidate windows. Evaluate starts daily and lengths on a weekly grid, always
% including MinDays and the full available span of each year/group.
rows=struct('StartDate',{},'EndDate',{},'Days',{},'NRows',{},'Coverage',{}, ...
    'SensorCoverage25',{},'PositiveFraction',{},'RelativeDeviation',{}, ...
    'TrendFraction',{},'CompleteSensors90',{},'Score',{},'PassesThresholds',{});

if o.SingleYear
    groupValues=unique(year(dates),'stable').';
else
    groupValues=0;
end

for gv=groupValues
    if gv==0
        idx=(1:numel(dates)).';
    else
        idx=find(year(dates)==gv);
    end
    if isempty(idx), continue; end
    dg=dates(idx); Yg=Y(idx,:);
    % Require actual calendar continuity. Missing archive days remain represented
    % by rows when using the Parma compact archive; this guard also supports
    % externally supplied irregular date vectors.
    nCalendar=days(dg(end)-dg(1))+1;
    if nCalendar<o.MinDays, continue; end
    maxLen=min(numel(dg),nCalendar);
    lengths=unique([o.MinDays, (o.MinDays+o.CandidateLengthStepDays):o.CandidateLengthStepDays:maxLen, maxLen]);
    lengths=lengths(lengths>=o.MinDays & lengths<=maxLen);
    for L=lengths
        for s=1:(numel(dg)-L+1)
            e=s+L-1;
            if days(dg(e)-dg(s))+1 ~= L
                continue
            end
            [met,pass]=windowMetricsLocal(dg(s:e),Yg(s:e,:),o);
            q=numel(rows)+1;
            rows(q).StartDate=dg(s); %#ok<AGROW>
            rows(q).EndDate=dg(e);
            rows(q).Days=L;
            rows(q).NRows=L;
            rows(q).Coverage=met.Coverage;
            rows(q).SensorCoverage25=met.SensorCoverage25;
            rows(q).PositiveFraction=met.PositiveFraction;
            rows(q).RelativeDeviation=met.RelativeDeviation;
            rows(q).TrendFraction=met.TrendFraction;
            rows(q).CompleteSensors90=met.CompleteSensors90;
            rows(q).Score=met.Score;
            rows(q).PassesThresholds=pass;
        end
    end
end

if isempty(rows)
    error('findStableParmaInterval:NoCandidate','No contiguous candidate interval of at least %d days was found.',o.MinDays);
end

T=struct2table(rows);
% Any passing interval outranks any failing interval. Score is primary within
% each group; a tiny duration term breaks near-ties in favour of longer spans.
priority=10*double(T.PassesThresholds)+T.Score+1e-7*T.Days;
[~,ord]=sort(priority,'descend');
T=T(ord,:);
best=T(1,:);

D=struct();
D.source=string(sourceLabel);
D.startDate=best.StartDate(1);
D.endDate=best.EndDate(1);
D.days=best.Days(1);
D.score=best.Score(1);
D.coverage=best.Coverage(1);
D.sensorCoverage25=best.SensorCoverage25(1);
D.positiveFraction=best.PositiveFraction(1);
D.relativeDeviation=best.RelativeDeviation(1);
D.trendFraction=best.TrendFraction(1);
D.completeSensors90=best.CompleteSensors90(1);
D.passesThresholds=logical(best.PassesThresholds(1));
D.nSensors=size(Y,2);
D.settings=o;
D.rankingTable=T(1:min(o.TopCandidates,height(T)),:);
D.figure=[];

fprintf('\n================ STABLE DATA INTERVAL ================\n');
fprintf('Best interval       : %s to %s (%d days)\n',datestr(D.startDate,'yyyy-mm-dd'),datestr(D.endDate,'yyyy-mm-dd'),D.days);
fprintf('Sensors scored      : %d\n',D.nSensors);
fprintf('Cell coverage       : %.2f%%\n',100*D.coverage);
fprintf('25th pct sensor cov.: %.2f%%\n',100*D.sensorCoverage25);
fprintf('Positive fraction   : %.2f%%\n',100*D.positiveFraction);
fprintf('Weekday-adj. spread : %.2f%%\n',100*D.relativeDeviation);
fprintf('Robust drift        : %.2f%%\n',100*D.trendFraction);
fprintf('Composite score     : %.4f\n',D.score);
fprintf('Passes thresholds   : %d\n',D.passesThresholds);
if ~D.passesThresholds
    fprintf('NOTE: no candidate met every requested threshold; the best available interval is shown.\n');
end
fprintf('======================================================\n');

if ~isempty(o.CSV)
    writetable(D.rankingTable,o.CSV);
    fprintf('Stable-interval candidate table written to: %s\n',o.CSV);
end

if o.Plot
    f=figure('Name','Parma stable-data interval diagnostic','Color','w','Visible',o.Visible);
    tl=tiledlayout(f,2,2,'TileSpacing','compact','Padding','compact'); %#ok<NASGU>
    inBest=dates>=D.startDate & dates<=D.endDate;
    nexttile;
    plot(dates,mean(isfinite(Y),2),'-','LineWidth',0.8); hold on;
    xline(D.startDate,'--'); xline(D.endDate,'--');
    ylim([0 1.02]); grid on; ylabel('Fraction sensors observed'); xlabel('Date');
    title('Daily measurement coverage');

    nexttile;
    net=weekdayAdjustedNetworkResidualLocal(dates,Y);
    plot(dates,100*(exp(net)-1),'-','LineWidth',0.8); hold on;
    xline(D.startDate,'--'); xline(D.endDate,'--');
    grid on; ylabel('Weekday-adjusted network deviation [%]'); xlabel('Date');
    title('Network-level stability');

    nexttile;
    scorePlot=T.Score; startPlot=T.StartDate;
    scatter(startPlot,scorePlot,12,T.Days,'filled'); hold on;
    scatter(D.startDate,D.score,70,'filled');
    grid on; xlabel('Candidate start'); ylabel('Quality score');
    title('Candidate windows (marker color = duration)'); colorbar;

    nexttile;
    Yb=Y(inBest,:);
    sensorCov=mean(isfinite(Yb),1);
    bar(sensorCov); ylim([0 1.02]); grid on;
    xlabel('Sensor'); ylabel('Coverage'); title(sprintf('Best interval sensor coverage (%d days)',D.days));
    D.figure=f;
end
end

function [dates,Y,label,mask]=resolveStableInputLocal(input)
mask=[];
if ischar(input)||isstring(input)
    file=char(string(input));
    if exist(file,'file')~=2, error('findStableParmaInterval:MissingFile','File not found: %s',file); end
    S=load(file);
    if isfield(S,'results') && isstruct(S.results)
        R=S.results;
    else
        fn=fieldnames(S); hit=false(numel(fn),1);
        for k=1:numel(fn), hit(k)=isstruct(S.(fn{k})); end
        ix=find(hit,1);
        if isempty(ix), error('findStableParmaInterval:MAT','MAT file does not contain a results structure.'); end
        R=S.(fn{ix});
    end
    label=file;
elseif isstruct(input)
    R=input; label='structure';
else
    error('findStableParmaInterval:Input','Input must be results, a MAT filename, a dates/Y structure, or a datetime vector followed by Y.');
end
if isfield(R,'selectionDates') && isfield(R,'Yselection')
    dates=R.selectionDates(:); Y=double(R.Yselection);
    if isfield(R,'sectionEstimationMask'), mask=R.sectionEstimationMask; end
elseif isfield(R,'dates') && isfield(R,'Y')
    dates=R.dates(:); Y=double(R.Y);
    if isfield(R,'sensorMask'), mask=R.sensorMask; end
else
    error('findStableParmaInterval:Fields','Structure must contain selectionDates/Yselection or dates/Y.');
end
end

function [M,pass]=windowMetricsLocal(dates,Y,o)
valid=isfinite(Y) & Y>=0;
M.Coverage=nnz(valid)/max(numel(valid),1);
sensorCoverage=mean(valid,1);
M.SensorCoverage25=percentileLocal(sensorCoverage,25);
vals=Y(valid);
if isempty(vals), M.PositiveFraction=0; else, M.PositiveFraction=nnz(vals>0)/numel(vals); end
M.CompleteSensors90=nnz(sensorCoverage>=0.90);
R=weekdayResidualMatrixLocal(dates,Y);
r=R(isfinite(R));
if isempty(r)
    M.RelativeDeviation=Inf;
else
    med=median(r,'omitnan');
    sigma=1.4826*median(abs(r-med),'omitnan');
    if ~isfinite(sigma), sigma=Inf; end
    M.RelativeDeviation=max(0,expm1(sigma));
end
net=median(R,2,'omitnan');
good=isfinite(net);
if nnz(good)>=max(14,ceil(0.5*numel(dates)))
    t=(0:numel(dates)-1).';
    pp=polyfit(t(good),net(good),1);
    driftSlope=abs(expm1(pp(1)*max(numel(dates)-1,1)));
    k=max(7,floor(numel(dates)/3));
    a=median(net(1:k),'omitnan'); b=median(net(end-k+1:end),'omitnan');
    if isfinite(a)&&isfinite(b), driftStep=abs(expm1(b-a)); else, driftStep=Inf; end
    M.TrendFraction=max(driftSlope,driftStep);
else
    M.TrendFraction=Inf;
end
stableScore=exp(-2*min(M.RelativeDeviation,3));
trendScore=exp(-3*min(M.TrendFraction,3));
lengthScore=min(1,max(0,(numel(dates)-o.MinDays)/max(365-o.MinDays,1)));
M.Score=0.45*M.Coverage + 0.20*M.SensorCoverage25 + 0.20*stableScore + 0.10*trendScore + 0.05*lengthScore;
pass=M.Coverage>=o.MinCoverage && M.SensorCoverage25>=o.MinSensorCoverage && ...
    M.PositiveFraction>=o.MinPositiveFraction && M.RelativeDeviation<=o.MaxRelativeDeviation && ...
    M.TrendFraction<=o.MaxTrendFraction;
end

function R=weekdayResidualMatrixLocal(dates,Y)
R=nan(size(Y)); wd=weekday(dates(:));
for j=1:size(Y,2)
    good=isfinite(Y(:,j)) & Y(:,j)>=0;
    if ~any(good), continue; end
    fallback=median(Y(good,j),'omitnan');
    for w=1:7
        ix=good & wd==w;
        if nnz(ix)>=2, base=median(Y(ix,j),'omitnan'); else, base=fallback; end
        if ~isfinite(base), continue; end
        R(ix,j)=log1p(Y(ix,j))-log1p(max(base,0));
    end
end
end

function net=weekdayAdjustedNetworkResidualLocal(dates,Y)
R=weekdayResidualMatrixLocal(dates,Y);
net=median(R,2,'omitnan');
end

function q=percentileLocal(x,p)
x=sort(double(x(isfinite(x))));
if isempty(x), q=NaN; return; end
if numel(x)==1, q=x; return; end
u=1+(numel(x)-1)*p/100;
a=floor(u); b=ceil(u);
if a==b, q=x(a); else, q=x(a)+(u-a)*(x(b)-x(a)); end
end

function tf=isPositiveIntegerLocal(x)
tf=isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=1&&x==floor(x);
end
function tf=isUnitIntervalLocal(x)
tf=isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=0&&x<=1;
end
function tf=isScalarLogicalLikeLocal(x)
tf=(islogical(x)||isnumeric(x))&&isscalar(x)&&isfinite(double(x));
end
