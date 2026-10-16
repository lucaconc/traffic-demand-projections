function out = generateParmaRequestedFigures(R,varargin)
%GENERATEPARMAREQUESTEDFIGURES Five figures from runParmaTrafficTomography.
% out = generateParmaRequestedFigures(results)
% out = generateParmaRequestedFigures('parma_tomography_results.mat')
% Requires MATLAB R2020b or newer. No Mapping/Statistics toolbox, internet,
% original traffic archive, or estimator is required. Maps use the saved OSM
% road graph (longitude/latitude with local aspect correction).
%
% Name/value options:
%   OutputDir       'parma_requested_figures' (PNG, vector PDF, FIG, CSV)
%   Visible         'on' or 'off'
%   PathWeighting   'equal' (default) or 'hits' (sampling multiplicities)
%   Validation      'section' (default) or 'temporal' (explicit alternative)
%   LabelSensors    false; label points with saved SectionID
%   LabelStrata     false; optionally label rectangles with stratum IDs
%   SensorArrowScale 0.04 (fraction of map span; formerly 0.018)
%   ProfileSensorRow [] selects highest-flow usable validation sensor, or
%                    highest-flow estimation sensor if none is usable.
%                    Set a measured-row index to choose a specific sensor.
%   FontSize        16 points for tick labels, legends, and annotations
%   LabelFontSize   18 points for axis and colorbar labels
%   Resolution      300 DPI
%   DisplayFence    3: the error panels of Figure 4 show the values within the far-outlier
%                   fences [Q1-3*IQR, Q3+3*IQR] (quartiles of all values); values beyond
%                   are genuine outliers, drawn as triangles on the limits with their
%                   number. 0 disables the fences (then DisplayPercentiles applies).
%   DisplayPercentiles [0 100]: alternative percentile limits for the error panels of Figure 4: the values between these
%                   percentiles (pooled over all day types); values beyond are drawn as
%                   markers on the axis limit, with their number printed. Statistics
%                   (medians, means, CSV, paper numbers) always use all values.
%                   [0 100] shows everything.
%   MinTestDays     2: in the temporal validation, only strata with at least two test days
%                   are shown; a single day's mean mostly reflects that day's events
%   MinObservedPerHour 30: the local relative error is computed only for sensor-stratum
%                   pairs observing at least 30 vehicles/hour; below that, dividing by a
%                   tiny count produces huge percentages that say nothing about the fit
%   FlowBands       8 discrete color bands in the peak-flow maps (Figure 2)
%   FlowBandScale   'quantile' (default) or 'linear'. Band boundaries are
%                   common to the three maps and are printed on the color
%                   bar, so every edge color matches exactly one band.
%   FlowDirection   'combined' (default): one line per road link, colored by
%                   the sum of its two directed flows, so every drawn line has
%                   exactly one color of the color bar. 'directed': one line
%                   per directed edge, the two directions drawn side by side.
%   DirectionOffset true; with FlowDirection='directed', offset the two
%                   directions to the right of travel
%   ShowNNLS        true; in Figure 4, mark the medians of the NNLS pilot
%                   prediction (red diamonds) when the results contain
%                   sectionValidationPredictedMeanNNLSByStratum
%
% Definitions:
% * Figure 2 shows separate weekday, Saturday, and Sunday/holiday peaks,
%   each maximizing total latent hourly flow within that day category.
%   The legacy overall-peak output fields/CSVs are also retained.
% * A class rate is SPLIT across ALL SAVED distinct sampled member paths;
%   this is a visualization assumption, not an identified within-class fit.
%   The saved paths are not an exhaustive enumeration of graph paths.
% * Shared directed edges accumulate the contributions of all member paths.
% * Validation samples in the violins are sensor-by-stratum MEAN errors,
%   predicted minus observed, equally weighted; not individual vehicle counts.
% * Figure 4 shows LOCAL relative errors, 100*(predicted-observed)/observed
%   at each sensor (paper: xi^loc).  The CSV also keeps the errors divided by
%   total LATENT predicted flow in the stratum (paper: xi^rel), which are
%   much smaller and not comparable across sensors.
% * Sunday includes holidays if the estimation run grouped them together.
% * Sensor positions/directions are the saved, possibly lane-merged metadata.
%   Missing coordinates cannot be mapped; rejected matches get no arrows.
%
% out contains figure handles, peak stratum, edge rates, path allocations,
% validation values, and exported filenames. Inputs are never modified.

if nargin==0, R='parma_tomography_results.mat'; end
p=inputParser;
addParameter(p,'OutputDir','parma_requested_figures',@(x)ischar(x)||isstring(x));
addParameter(p,'Visible','on',@(x)any(strcmpi(x,{'on','off'})));
addParameter(p,'PathWeighting','equal',@(x)any(strcmpi(x,{'equal','hits'})));
addParameter(p,'Validation','auto',@(x)any(strcmpi(x,{'auto','section','temporal'})));
addParameter(p,'LabelSensors',false,@(x)islogical(x)&&isscalar(x));
addParameter(p,'LabelStrata',false,@(x)islogical(x)&&isscalar(x));
addParameter(p,'SensorArrowScale',0.04,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0);
addParameter(p,'ProfileSensorRow',[],@(x)isempty(x)||(isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=1&&x==fix(x)));
addParameter(p,'FontSize',16,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0);
addParameter(p,'LabelFontSize',18,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0);
addParameter(p,'Resolution',300,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'ShowNNLS',true,@(x)islogical(x)&&isscalar(x));
addParameter(p,'DisplayFence',3,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
addParameter(p,'DisplayPercentiles',[0 100],@(x)isnumeric(x)&&numel(x)==2&&x(1)>=0&&x(2)<=100&&x(1)<x(2));
addParameter(p,'MinTestDays',2,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
addParameter(p,'MinObservedPerHour',30,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
addParameter(p,'FlowBands',8,@(x)isnumeric(x)&&isscalar(x)&&x>=2&&x==fix(x));
addParameter(p,'FlowBandScale','quantile',@(x)any(strcmpi(x,{'quantile','linear'})));
addParameter(p,'DirectionOffset',true,@(x)islogical(x)&&isscalar(x));
addParameter(p,'FlowDirection','combined',@(x)any(strcmpi(x,{'combined','directed'})));
parse(p,varargin{:}); o=p.Results;
% All output paths lie outside the library (see safeOutputPath).
o=guardOutputFields(o,{'OutputDir'});
if ischar(R)||isstring(R)
    data=load(R);
    if isfield(data,'results'), R=data.results;
    else, error('ParmaFigures:Input','MAT file must contain variable results.'); end
end
need={'G','sectionMatch','strataTable','muRatePerHourByStratum','classEdgePaths'};
assert(isstruct(R)&&all(isfield(R,need)), ...
    'Input must be the results structure from the supplied tomography file.');
G=R.G; T=R.strataTable; rates=double(R.muRatePerHourByStratum);
assert(size(rates,2)==height(T),'Strata and rate columns do not agree.');
good=all(isfinite(rates)&rates>=0,1).';
if ismember('Estimated',T.Properties.VariableNames), good=good & logical(T.Estimated); end
assert(any(good),'No successfully fitted stratum is available.');
total=sum(rates,1).'; total(~good)=NaN;
score=total; score(~good)=-Inf; [~,peak]=max(score);
[s,t]=endpoints(G); lon=double(G.Nodes.Lon); lat=double(G.Nodes.Lat);
category=dayCategory(T.DayType);
[edgeRate,pathTable]=allocatePaths(R,rates(:,peak),o.PathWeighting);
if strcmpi(o.Validation,'auto')
    % Section hold-out when sensors were held out, otherwise the temporal hold-out.
    if isfield(R,'sectionValidationMask') && any(R.sectionValidationMask), o.Validation='section';
    else, o.Validation='temporal'; end
end
[V,validationLabel]=validationValues(R,T,total,good,o.Validation,category,o.MinObservedPerHour,o.MinTestDays);
if ~isfolder(o.OutputDir), mkdir(o.OutputDir); end
out=struct(); out.peakStratum=T(peak,:); out.peakTotalPerHour=total(peak);
out.pathAllocation=pathTable; out.validation=V; out.totalPerHour=total;
out.edgeFlows=table((1:numedges(G)).',s,t,edgeRate, ...
    'VariableNames',{'EdgeID','SourceNode','TargetNode','VehiclesPerHour'});
figs=gobjects(5,1);

%% 1. Selected estimation/validation sensors with accepted inferred directions.
figs(1)=newFigure('Parma sensors and flow directions',o,[1200 850]);
ax=axes(figs(1)); hold(ax,'on'); roadMap(ax,lon,lat,s,t,1.2,[.48 .51 .54]);
M=R.sectionMatch; x=double(M.Lon); y=double(M.Lat);
gps=isfinite(x)&isfinite(y)&abs(x)<=180&abs(y)<=90;
% Keep geographically invalid/outlier metadata out of map extent.
inMap=gps & x>=min(lon)-0.02 & x<=max(lon)+0.02 & ...
    y>=min(lat)-0.02 & y<=max(lat)+0.02;
role=zeros(height(M),1); representative=(1:height(M)).';
if ismember('LaneMergedInto',M.Properties.VariableNames)
    representative=double(M.LaneMergedInto);
end
if isfield(R,'selectedSections')
    role(ismember(representative,R.selectedSections))=1;
    if isfield(R,'sectionValidationMask')
        v=R.selectedSections(logical(R.sectionValidationMask));
        role(ismember(representative,v))=2;
    end
end
colors=[0.48 0.48 0.48;0.05 0.38 0.75;0.9 0.33 0.08];
names={'Other saved sensors','Estimation sensors / lanes','Validation sensors / lanes'};
% Omit unselected sensor records (grey markers), including their labels.
displayedSensors=inMap & role>0;
for k=1:2
    z=displayedSensors & role==k;
    scatter(ax,x(z),y(z),30,colors(k+1,:),'filled', ...
        'MarkerEdgeColor','w','DisplayName',names{k+1});
end
arrowLength=o.SensorArrowScale*max((max(lon)-min(lon))*cosd(mean(lat)),max(lat)-min(lat));
nDirected=0;
% Draw direction arrows only for selected estimation/validation sensors.
for i=find(displayedSensors).'
    r=representative(i);
    if r<1||r>height(M)||~M.Accepted(r), continue; end
    ids=double(M.DirectedEdgeIDs{r}(:));
    ids=ids(isfinite(ids)&ids>=1&ids<=numedges(G)&ids==fix(ids));
    if isempty(ids), continue; end
    nDirected=nDirected+1;
    for e=ids.'
        dx=(lon(t(e))-lon(s(e)))*cosd(mean(lat)); dy=lat(t(e))-lat(s(e));
        d=hypot(dx,dy); if d==0, continue; end
        quiver(ax,x(i),y(i),arrowLength*dx/d/cosd(mean(lat)), ...
            arrowLength*dy/d,0,'Color',colors(role(i)+1,:), ...
            'LineWidth',2.0,'MaxHeadSize',2.5,'HandleVisibility','off');
    end
end
if o.LabelSensors
    text(ax,x(displayedSensors),y(displayedSensors),"  "+string(M.SectionID(displayedSensors)), ...
        'FontSize',7,'Interpreter','none','Clipping','on');
end
% Sensor roles are explained in the figure caption; omit the legend.


%% 2. Weekday, Saturday, and Sunday peaks on a common map and color scale.
peakNames={'Weekday','Saturday','Sunday / holiday'};
peakIndices=nan(3,1); dayEdges=zeros(numedges(G),3);
dayPaths=cell(3,1);
for c=1:3
    candidates=find(good & category==c);
    if isempty(candidates), continue; end
    [~,j]=max(total(candidates)); peakIndices(c)=candidates(j);
    [dayEdges(:,c),dayPaths{c}]=allocatePaths(R,rates(:,peakIndices(c)),o.PathWeighting);
end
out.peakStratumIndexByDayType=peakIndices;
out.peakEdgeFlowsByDayType=dayEdges;
out.peakPathAllocationsByDayType=dayPaths;
out.peakDayTypeNames=string(peakNames(:));
% An absent category is represented by NaN, not zero estimated traffic.
for c=1:3
    if isnan(peakIndices(c)), out.peakEdgeFlowsByDayType(:,c)=NaN; end
end
figs(2)=newFigure('Parma peak flows by day type',o,[1800 780]);
flowLayout=tiledlayout(figs(2),1,3,'TileSpacing','compact','Padding','compact');
% Discrete bands with common boundaries for the three maps. The axes use
% the band index as color coordinate (CLim=[0 K]) and the color bar shows the
% band boundaries as tick labels. By default each road link is drawn ONCE,
% colored by the sum of its two directed flows: overlapping or adjacent lines
% of different colors would otherwise blend at print size into colors that
% do not appear on the color bar.
combined=strcmpi(o.FlowDirection,'combined');
if combined
    [segS,segT,linkFlows]=combineDirections(s,t,dayEdges);
    flowLabel='Link flow, both directions [vehicles/hour]';
    offset=0;
else
    segS=s; segT=t; linkFlows=dayEdges;
    flowLabel='Directed-edge flow [vehicles/hour]';
    offset=0;
    if o.DirectionOffset
        span=max((max(lon)-min(lon))*cosd(mean(lat)),max(lat)-min(lat));
        offset=0.006*span;
    end
end
out.peakLinkFlowsByDayType=linkFlows;
out.peakLinkEndNodes=[segS(:) segT(:)];
flowAll=linkFlows(:); flowAll=flowAll(isfinite(flowAll)&flowAll>0);
bandEdges=flowBandEdges(flowAll,o.FlowBands,o.FlowBandScale);
K=numel(bandEdges)-1; cmap=flowColormap(K);
out.peakFlowBandEdges=bandEdges(:);
[x1,y1,x2,y2]=offsetSegments(lon,lat,segS,segT,offset);
for c=1:3
    ax=nexttile(flowLayout); hold(ax,'on'); roadMap(ax,lon,lat,s,t);
    colormap(ax,cmap); caxis(ax,[0 K]);
    % Explicit common extent permits direct visual comparison.
    xlim(ax,[min(lon) max(lon)]); ylim(ax,[min(lat) max(lat)]);
    if isnan(peakIndices(c))
        text(ax,.5,.5,'No fitted flow available','Units','normalized', ...
            'HorizontalAlignment','center','BackgroundColor','w');
        continue
    end
    band=flowBandIndex(linkFlows(:,c),bandEdges);
    % Lower bands first, so that larger flows are drawn on top.
    for b=1:K
        ids=find(band==b); if isempty(ids), continue; end
        xx=[x1(ids) x2(ids) nan(numel(ids),1)].';
        yy=[y1(ids) y2(ids) nan(numel(ids),1)].';
        plot(ax,xx(:),yy(:),'Color',cmap(b,:),'LineWidth',1.5+2.5*b/K);
    end
end
cb=colorbar(ax); cb.Layout.Tile='east'; cb.Label.String=flowLabel;
cb.Ticks=0:K; cb.TickLabels=compose('%g',roundFlowLabels(bandEdges));
writetable(table((1:K).',bandEdges(1:end-1).',bandEdges(2:end).',cmap(:,1),cmap(:,2),cmap(:,3), ...
    'VariableNames',{'Band','LowerVehiclesPerHour','UpperVehiclesPerHour','R','G','B'}), ...
    fullfile(o.OutputDir,'peak_flow_color_bands.csv'));


%% 3. Every stratum, including skipped fits, without overlapping weekdays/years.
keys=string(T.DayType)+" | "+string(T.YearLabel);
rowKeys=cell(3,1); nrows=1;
for c=1:3
    rowKeys{c}=unique(keys(category==c),'stable'); nrows=max(nrows,numel(rowKeys{c}));
end
out.strataPanelKeys=rowKeys;
figs(3)=newFigure('All Parma strata',o,[1500 max(600,330*nrows)]);
tl=tiledlayout(figs(3),nrows,3,'TileSpacing','compact','Padding','compact');
headings={'Weekday','Saturday','Sunday / holiday'}; maxTotal=max(total(good));
if maxTotal<=0, maxTotal=1; end
monthTick=days(datetime(2000,1:12,1)-datetime(2000,1,1))+1;
for r=1:nrows
    for c=1:3
        ax=nexttile(tl); hold(ax,'on'); colormap(ax,parula(128)); caxis(ax,[0 maxTotal]);
        if r>numel(rowKeys{c})
            axis(ax,'off');  continue
        end
        ids=find(category==c & keys==rowKeys{c}(r));
        for q=ids.'
            a=double(T.StartMinute(q))/60; b=double(T.EndMinute(q))/60;
            lo=double(T.SeasonStartDOY(q)); hi=double(T.SeasonEndDOY(q))+1;
            if good(q), color=parula(128); color=color(max(1,ceil(128*total(q)/maxTotal)),:);
            else, color=[0.8 0.8 0.8]; end
            patch(ax,[a b b a],[lo lo hi hi],color,'EdgeColor',[.25 .25 .25],'LineWidth',0.6);
            if o.LabelStrata
                text(ax,(a+b)/2,(lo+hi)/2,string(T.Stratum(q)), ...
                    'HorizontalAlignment','center','FontSize',7,'Clipping','on');
            end
        end
        xlim(ax,[0 24]); ylim(ax,[1 367]); xticks(ax,0:4:24);
        yticks(ax,monthTick); yticklabels(ax,{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
        xlabel(ax,'Time of day [h]'); ylabel(ax,'Calendar date');
        
        caxis(ax,[0 maxTotal]); box(ax,'on');
    end
end
cb=colorbar(ax); cb.Layout.Tile='east'; cb.Label.String='Total latent flow [vehicles/hour]';

%% 4. Validation mean-error distributions and overall predicted flow.
figs(4)=newFigure('Validation mean errors versus predicted flow',o,[1550 580]);
tl=tiledlayout(figs(4),1,3,'TileSpacing','compact','Padding','compact');
hasNNLS=o.ShowNNLS && ismember('ErrorNNLSPerHour',V.Properties.VariableNames);
ax=nexttile(tl);
if hasNNLS, violinGroups(ax,V.ErrorPerHour,V.Category,V.ErrorNNLSPerHour,o.DisplayPercentiles,o.DisplayFence);
else, violinGroups(ax,V.ErrorPerHour,V.Category,[],o.DisplayPercentiles,o.DisplayFence); end
ylabel(ax,'Predicted - observed mean [vehicles/hour]');
ax=nexttile(tl);
if hasNNLS, violinGroups(ax,V.PercentOfObserved,V.Category,V.PercentOfObservedNNLS,o.DisplayPercentiles,o.DisplayFence);
else, violinGroups(ax,V.PercentOfObserved,V.Category,[],o.DisplayPercentiles,o.DisplayFence); end
ylabel(ax,'Local relative error [%]');
ax=nexttile(tl); hold(ax,'on');
for c=1:3
    qq=unique(V.StratumIndex(V.Category==c)); mae=nan(numel(qq),1);
    for j=1:numel(qq), mae(j)=mean(abs(V.ErrorPerHour(V.StratumIndex==qq(j)))); end
    scatter(ax,total(qq),mae,28,colors(c,:),'filled','DisplayName',headings{c});
end
xlabel(ax,'Overall predicted latent flow [vehicles/hour]'); ylabel(ax,'Validation MAE [vehicles/hour]');
 grid(ax,'on'); legend(ax,'Location','best');
if isempty(V), text(ax,.5,.5,'No valid validation pairs','Units','normalized','HorizontalAlignment','center'); end


%% 5. Weekday mean profile for a significant sensor.
[profile,sensor]=weekdayProfile(R,T,good,category,o.ProfileSensorRow);
out.weekdayProfile=profile; out.profileSensor=sensor;
figs(5)=newFigure('Weekday predicted and observed mean flow',o,[1200 650]);
ax=axes(figs(5)); hold(ax,'on');
if isempty(profile)
    text(ax,.5,.5,'No fitted weekday observations available', ...
        'Units','normalized','HorizontalAlignment','center');
    
else
    % Independent horizontal segments avoid interpolating across missing data.
    xx=[profile.StartHour profile.EndHour nan(height(profile),1)].';
    realY=[profile.ObservedPerHour profile.ObservedPerHour nan(height(profile),1)].';
    predY=[profile.PredictedPerHour profile.PredictedPerHour nan(height(profile),1)].';
    plot(ax,xx(:),realY(:),'-','Color',[.05 .38 .75],'LineWidth',2.4, ...
        'DisplayName','Observed mean flow');
    plot(ax,xx(:),predY(:),'--','Color',[.9 .33 .08],'LineWidth',2.4, ...
        'DisplayName','Predicted mean flow');
    
    legend(ax,'Location','best');
end
xlabel(ax,'Time of day [h]'); ylabel(ax,'Mean flow [vehicles/hour]');
xlim(ax,[0 24]); xticks(ax,0:2:24); grid(ax,'on'); box(ax,'on');

%% Export figures and the numerical values underlying them.
stems={'01_sensors_directions','02_peak_flows_by_day_type','03_all_strata','04_validation_violins','05_weekday_sensor_profile'};
files=strings(5,3);
for k=1:5
    requireFigure(figs(k),sprintf('before exporting figure %d',k));
    removeFigureTitles(figs(k));
    blackFigureText(figs(k),o);
    base=fullfile(o.OutputDir,stems{k});
    files(k,:)=[string([base '.png']) string([base '.pdf']) string([base '.fig'])];
    exportParmaFigure(figs(k),files(k,1),files(k,2),o.Resolution);
    savefig(figs(k),files(k,3));
end
writetable(out.edgeFlows,fullfile(o.OutputDir,'peak_edge_flows.csv'));
writetable(pathTable,fullfile(o.OutputDir,'peak_path_allocations.csv'));
writetable(V,fullfile(o.OutputDir,'validation_mean_errors.csv'));
writetable(out.peakStratum,fullfile(o.OutputDir,'peak_stratum.csv'));
writetable(profile,fullfile(o.OutputDir,'weekday_sensor_profile.csv'));
% Keep legacy overall-peak CSVs, and export each plotted day-type peak.
peakTags={'weekday','saturday','sunday'};
for c=1:3
    if isnan(peakIndices(c)), continue; end
    E=out.edgeFlows; E.VehiclesPerHour=dayEdges(:,c);
    writetable(E,fullfile(o.OutputDir,[peakTags{c} '_peak_edge_flows.csv']));
    writetable(dayPaths{c},fullfile(o.OutputDir,[peakTags{c} '_peak_path_allocations.csv']));
    writetable(T(peakIndices(c),:),fullfile(o.OutputDir,[peakTags{c} '_peak_stratum.csv']));
end
for k=1:numel(figs)
    requireFigure(figs(k),'before displaying completed figures');
    set(figs(k),'Visible',char(o.Visible));
end
out.figures=figs; out.files=files; out.options=o;
fprintf('Saved five figures to %s. Peak stratum %d: %.3f vehicles/hour.\n', ...
    char(o.OutputDir),T.Stratum(peak),total(peak));
end

function [s,t]=endpoints(G)
e=G.Edges.EndNodes;
if isnumeric(e), s=double(e(:,1)); t=double(e(:,2));
else, s=findnode(G,e(:,1)); t=findnode(G,e(:,2)); end
s=s(:); t=t(:);
assert(all(s>0&t>0),'Invalid graph endpoints.');
end

function f=newFigure(name,o,sz)
% Build offscreen so closing a preview cannot invalidate a pending export.
% Explicit normal windows avoid depending on the user's docking defaults.
f=figure('Name',name,'Color','w','Visible','off', ...
    'WindowStyle','normal','Units','pixels','Position',[40 40 sz]);
setappdata(f,'ParmaExportSizePixels',sz);
requireFigure(f,'after creating a figure');
end

function roadMap(ax,lon,lat,s,t,roadWidth,roadColor)
if nargin<7, roadWidth=.45; end
if nargin<8, roadColor=[.83 .85 .87]; end
x=[lon(s) lon(t) nan(numel(s),1)].'; y=[lat(s) lat(t) nan(numel(s),1)].';
plot(ax,x(:),y(:),'Color',roadColor,'LineWidth',roadWidth,'HandleVisibility','off');
daspect(ax,[1/cosd(mean(lat)) 1 1]); axis(ax,'tight'); box(ax,'on');
xlabel(ax,'Longitude [degrees]'); ylabel(ax,'Latitude [degrees]');

end

function [segS,segT,flows]=combineDirections(s,t,dayEdges)
% One segment per unordered node pair; flows are summed over both directions
% (and over parallel edges between the same nodes). NaN columns stay NaN.
key=sort([s(:) t(:)],2);
[u,~,g]=unique(key,'rows');
segS=u(:,1); segT=u(:,2);
flows=zeros(size(u,1),size(dayEdges,2));
for c=1:size(dayEdges,2)
    v=dayEdges(:,c);
    if all(isnan(v)), flows(:,c)=NaN; continue; end
    v(~isfinite(v))=0;
    flows(:,c)=accumarray(g(:),v(:),[size(u,1) 1]);
end
end

function edges=flowBandEdges(flow,K,scale)
% Band boundaries for positive flows: equal counts ('quantile') or equal
% widths ('linear'), from 0 to the maximum. Duplicated quantiles are merged.
flow=sort(flow(:)); mx=max([flow;1]);
if isempty(flow), edges=[0 mx]; return; end
if strcmpi(scale,'linear')
    edges=linspace(0,mx,K+1);
else
    q=(1:K-1)/K; n=numel(flow);
    inner=flow(max(1,min(n,round(q*n)))).';
    edges=[0 inner mx];
end
edges=unique(edges);
if numel(edges)<2, edges=[0 mx]; end
end

function band=flowBandIndex(flow,edges)
% Band of each flow; 0 for zero or missing flow (not drawn). The upper
% boundary of the last band is included.
K=numel(edges)-1; band=zeros(size(flow));
ok=isfinite(flow)&flow>0;
for b=1:K
    if b<K, in=flow>edges(b)&flow<=edges(b+1); else, in=flow>edges(b); end
    band(ok&in)=b;
end
end

function labels=roundFlowLabels(edges)
% Two significant digits for readability; the exact boundaries are saved
% in peak_flow_color_bands.csv and used for the colors.
labels=edges;
nz=edges>0;
labels(nz)=round(edges(nz),2,'significant');
end

function cmap=flowColormap(K)
% Sequential map from light orange to dark purple, readable in print and
% visible against the grey road network.
anchors=[0.99 0.80 0.30; 0.98 0.52 0.20; 0.86 0.22 0.15; 0.62 0.04 0.28; 0.30 0.00 0.36];
t=linspace(0,1,size(anchors,1)); u=linspace(0,1,K);
cmap=[interp1(t,anchors(:,1),u).' interp1(t,anchors(:,2),u).' interp1(t,anchors(:,3),u).'];
end

function [x1,y1,x2,y2]=offsetSegments(lon,lat,s,t,offset)
% Shift each directed edge to the right of its travel direction, so the two
% directions of a road do not overlap. Offsets are isotropic on the map.
x1=lon(s); y1=lat(s); x2=lon(t); y2=lat(t);
if offset<=0, return; end
k=cosd(mean(lat));
dx=(x2-x1)*k; dy=y2-y1; d=hypot(dx,dy); d(d==0)=Inf;
nx=dy./d; ny=-dx./d;                 % unit normal to the right of travel
x1=x1+offset*nx/k; x2=x2+offset*nx/k;
y1=y1+offset*ny;   y2=y2+offset*ny;
end

function c=dayCategory(day)
d=lower(string(day)); c=zeros(numel(d),1);
c(contains(d,'saturday'))=2; c(contains(d,'sunday')|contains(d,'holiday'))=3;
week=ismember(d,["workday","weekday","monday","tuesday","wednesday","thursday","friday"]);
c(week)=1;
assert(all(c>0),'Unrecognized day type in strataTable.');
end

function [edgeRate,P]=allocatePaths(R,mu,weighting)
n=numedges(R.G); edgeRate=zeros(n,1); C=R.classEdgePaths;
assert(numel(C)==numel(mu),'Class paths and rate rows do not agree.');
cl=[]; pa=[]; wt=[]; fl=[];
for j=1:numel(mu)
    paths=C{j}; assert(iscell(paths)&&~isempty(paths),'Class %d has no saved member paths.',j);
    w=ones(numel(paths),1);
    if strcmpi(weighting,'hits')
        assert(isfield(R,'classPathHits'),'Hits weighting requires classPathHits.');
        w=double(R.classPathHits{j}(:));
        assert(numel(w)==numel(paths)&&all(isfinite(w)&w>0),'Invalid classPathHits.');
    end
    w=w/sum(w);
    for k=1:numel(paths)
        ids=double(paths{k}(:));
        assert(~isempty(ids)&&all(isfinite(ids)&ids==fix(ids)&ids>=1&ids<=n), ...
            'Invalid saved edge path for class %d, member %d.',j,k);
        % ACCUMARRAY also correctly counts repeated traversal of an edge.
        edgeRate=edgeRate+accumarray(ids,repmat(mu(j)*w(k),numel(ids),1),[n 1]);
    end
    cl=[cl;repmat(j,numel(w),1)]; pa=[pa;(1:numel(w)).']; %#ok<AGROW>
    wt=[wt;w]; fl=[fl;mu(j)*w]; %#ok<AGROW>
end
assert(abs(sum(fl)-sum(mu))<=1e-8*max(1,sum(mu)),'Path allocation failed conservation.');
P=table(cl,pa,wt,fl,'VariableNames',{'Class','MemberPath','Fraction','VehiclesPerHour'});
end

function [V,label]=validationValues(R,T,total,good,mode,category,minObs,minTestDays)
if nargin<7, minObs=0; end
if nargin<8, minTestDays=1; end
if strcmpi(mode,'section')
    fields={'sectionValidationMeasuredMeanByStratum','sectionValidationPredictedMeanByStratum', ...
        'sectionValidationObservationCountByStratum','sectionValidationMask'};
    assert(all(isfield(R,fields)),'Section-validation arrays are missing.');
    obs=double(R.sectionValidationMeasuredMeanByStratum);
    pred=double(R.sectionValidationPredictedMeanByStratum);
    count=double(R.sectionValidationObservationCountByStratum);
    allowed=repmat(logical(R.sectionValidationMask(:)),1,height(T));
    label='held-out sensors';
else
    obs=double(R.testMeasuredMeanByStratum); pred=double(R.predictedMeanByStratum);
    count=repmat(double(T.NTest(:)).',size(obs,1),1);
    allowed=logical(R.sensorUsedByStratum); label='temporal hold-out on estimation sensors';
    if isfield(R,'testDaysByStratum')
        allowed=allowed & repmat(R.testDaysByStratum(:).'>=minTestDays,size(allowed,1),1);
    end
end
assert(isequal(size(obs),size(pred),size(count),size(allowed))&&size(obs,2)==height(T), ...
    'Validation arrays have inconsistent dimensions.');
minutes=double(T.ObservationMinutes(:)); assert(all(isfinite(minutes)&minutes>0),'Invalid observation intervals.');
valid=allowed & isfinite(obs)&isfinite(pred)&count>0 & repmat(good.',size(obs,1),1);
[i,q]=find(valid); linear=sub2ind(size(obs),i,q); factor=60./minutes(q);
observed=obs(linear).*factor; predicted=pred(linear).*factor; err=predicted-observed;
den=total(q); pct=nan(size(err)); nonzero=den>0; pct(nonzero)=100*err(nonzero)./den(nonzero);
loc=nan(size(err)); pos=observed>0 & observed>=minObs; loc(pos)=100*err(pos)./observed(pos);
V=table(i,q,T.Stratum(q),category(q),count(linear),observed,predicted,err,den,pct,loc, ...
    'VariableNames',{'MeasuredRow','StratumIndex','Stratum','Category','ObservationCount', ...
    'ObservedPerHour','PredictedPerHour','ErrorPerHour','OverallLatentPerHour','PercentOfOverallFlow', ...
    'PercentOfObserved'});
% Baseline: held-out prediction of the NNLS fit of the training means
% (stored by the modified runParmaTrafficTomography, section mode only).
predN=[];
if strcmpi(mode,'section') && isfield(R,'sectionValidationPredictedMeanNNLSByStratum')
    predN=double(R.sectionValidationPredictedMeanNNLSByStratum);
elseif strcmpi(mode,'temporal') && isfield(R,'predictedMeanNNLSByStratum')
    predN=double(R.predictedMeanNNLSByStratum);
end
if ~isempty(predN)
    if isequal(size(predN),size(obs))
        pN=predN(linear).*factor; eN=pN-observed;
        locN=nan(size(eN)); locN(pos)=100*eN(pos)./observed(pos);
        V.PredictedNNLSPerHour=pN; V.ErrorNNLSPerHour=eN; V.PercentOfObservedNNLS=locN;
    end
end
if isempty(V)
    warning('ParmaFigures:NoValidation','No valid %s pairs. Figure 4 will be empty. Use Validation=temporal explicitly if appropriate.',label);
end
end

function violinGroups(ax,y,category,ref,pct,fence)
% ref (optional): baseline values; their medians are drawn as red diamonds.
% pct (optional): display limits as percentiles of all finite y; values beyond
% them are drawn as small markers on the limit and counted, but medians and
% means are computed from all values.
if nargin<4, ref=[]; end
if nargin<5 || isempty(pct), pct=[0 100]; end
if nargin<6 || isempty(fence), fence=0; end
hold(ax,'on'); palette=[.3 .35 .45;.05 .38 .75;.9 .33 .08;.2 .65 .4];
yy=sort(y(isfinite(y)));
if isempty(yy), return; end
if fence>0
    % Far-outlier fences: only genuinely extreme values fall outside
    q1=pctLocal(yy,25); q3=pctLocal(yy,75); iq=q3-q1;
    lim=[max(min(yy),q1-fence*iq) min(max(yy),q3+fence*iq)];
else
    lim=[pctLocal(yy,pct(1)) pctLocal(yy,pct(2))];
end
if lim(2)<=lim(1), lim=[min(yy) max(yy)]; end
pad=0.06*max(lim(2)-lim(1),eps);
nOut=0;
labels={'Overall','Weekday','Saturday','Sunday/holiday'};
for k=1:4
    if k==1, v=y; else, v=y(category==k-1); end
    v=v(isfinite(v));
    if isempty(v), continue; end
    vAll=v;                                   % for the median and mean
    lowOut=v<lim(1); highOut=v>lim(2);
    v=v(~lowOut & ~highOut);                  % shape drawn from the displayed range
    if any(lowOut)
        plot(ax,k+.18,lim(1)-0.5*pad,'v','Color',palette(k,:),'MarkerFaceColor',palette(k,:),'MarkerSize',5);
        text(ax,k+.27,lim(1)-0.5*pad,sprintf('%d',nnz(lowOut)),'FontSize',8,'VerticalAlignment','middle');
    end
    if any(highOut)
        plot(ax,k+.18,lim(2)+0.5*pad,'^','Color',palette(k,:),'MarkerFaceColor',palette(k,:),'MarkerSize',5);
        text(ax,k+.27,lim(2)+0.5*pad,sprintf('%d',nnz(highOut)),'FontSize',8,'VerticalAlignment','middle');
    end
    if k==1, nOut=nnz(lowOut)+nnz(highOut); end
    if isempty(v), v=vAll; end
    if numel(v)>1 && max(v)>min(v)
        % Gaussian KDE implemented locally: no ksdensity/violinplot dependency.
        h=1.06*std(v)*numel(v)^(-.2); h=max(h,eps(max(abs(v))));
        gridY=linspace(min(v)-2.5*h,max(v)+2.5*h,180).'; density=zeros(size(gridY));
        for j=1:1000:numel(v)
            z=(gridY-v(j:min(j+999,numel(v))).')/h;
            density=density+sum(exp(-.5*z.^2),2);
        end
        width=.36*density/max(density);
        patch(ax,[k-width;flipud(k+width)],[gridY;flipud(gridY)],palette(k,:), ...
            'FaceAlpha',.55,'EdgeColor',palette(k,:));
    else
        plot(ax,[k-.22 k+.22],[v(1) v(1)],'Color',palette(k,:),'LineWidth',3);
    end
    md=median(vAll); mn=mean(vAll);
    plot(ax,[k-.12 k+.12],[md md],'k-','LineWidth',1.7);
    if mn>=lim(1) && mn<=lim(2), plot(ax,k,mn,'ko','MarkerFaceColor','w','MarkerSize',5); end
    if ~isempty(ref)
        if k==1, r=ref; else, r=ref(category==k-1); end
        r=r(isfinite(r));
        if ~isempty(r)
            plot(ax,k,median(r),'d','Color',[.8 .1 .1],'MarkerFaceColor',[.8 .1 .1],'MarkerSize',8);
        end
    end
end
yline(ax,0,'k:'); xlim(ax,[.5 4.5]); xticks(ax,1:4); xticklabels(ax,labels);
grid(ax,'on'); box(ax,'on');
if fence>0 || any(pct~=[0 100])
    ylim(ax,[lim(1)-pad,lim(2)+pad]);
end
end


function blackFigureText(f,o)
% Explicit colors override MATLAB theme/defaults, including layout titles.
set(findall(f,'Type','text'),'Color','k','FontUnits','points','FontSize',o.FontSize);
a=findall(f,'Type','axes');
for k=1:numel(a)
    set(a(k),'XColor','k','YColor','k','ZColor','k','Color','w','FontUnits','points','FontSize',o.FontSize);
    set([a(k).XLabel a(k).YLabel a(k).ZLabel],'Color','k','FontUnits','points','FontSize',o.LabelFontSize);
end
set(findall(f,'Type','legend'),'TextColor','k','Color','w','FontSize',o.FontSize);
b=findall(f,'Type','colorbar');
for k=1:numel(b)
    b(k).Color='k'; b(k).FontSize=o.FontSize;
    b(k).Label.Color='k'; b(k).Label.FontSize=o.LabelFontSize;
end
layouts=findall(f,'Type','tiledlayout');
for k=1:numel(layouts), layouts(k).Title.Color='k'; end
end

function [P,sensor]=weekdayProfile(R,T,good,category,chosen)
% The producer saves observation dates but not each row's native time slot.
% Therefore reconstruct a PIECEWISE STRATUM-MEAN profile, not an exact raw
% hourly average. For each shared time segment, weight each covering stratum
% by valid observation count / stratum duration. This avoids overweighting
% broad time strata. Both curves use identical strata and weights.
% Missing observations may not be uniform within a stratum; exact hourly
% empirical means would require the original archive or saved row timestamps.
P=table('Size',[0 7],'VariableTypes',repmat({'double'},1,7), ...
    'VariableNames',{'MeasuredRow','StartHour','EndHour','ObservedPerHour', ...
    'PredictedPerHour','ObservationDensityPerMinute','ContributingStrata'});
sensor=struct();
assert(isfield(R,'strata')&&isfield(R,'predictedMeanByStratum'), ...
    'Weekday profile requires saved strata and predictedMeanByStratum.');
q=find(good & category==1); if isempty(q), return; end
m=size(R.predictedMeanByStratum,1); n=numel(q);
obs=nan(m,n); pred=nan(m,n); counts=zeros(m,n);
for j=1:n
    Y=double(R.strata(q(j)).YFull); Y(~isfinite(Y)|Y<0)=NaN;
    counts(:,j)=sum(isfinite(Y),1).';
    obs(:,j)=mean(Y,1,'omitnan').'*60/double(T.ObservationMinutes(q(j)));
    pred(:,j)=double(R.predictedMeanByStratum(:,q(j)))*60/double(T.ObservationMinutes(q(j)));
end
valid=isfinite(obs)&isfinite(pred)&counts>0;
counts(~valid)=0; tmp=obs; tmp(~valid)=0;
score=sum(tmp.*counts,2)./max(sum(counts,2),1);
score(sum(counts,2)==0)=-Inf;
if isempty(chosen)
    candidates=find(isfinite(score));
    if isfield(R,'sectionValidationMask')
        v=candidates(logical(R.sectionValidationMask(candidates)));
        if ~isempty(v), candidates=v; end
    end
    if isempty(candidates), return; end
    [~,j]=max(score(candidates)); chosen=candidates(j);
else
    assert(chosen<=m&&isfinite(score(chosen)), ...
        'ProfileSensorRow has no valid fitted weekday data.');
end
section=R.selectedSections(chosen); role="Estimation sensor (descriptive fit comparison)";
if isfield(R,'sectionValidationMask')&&R.sectionValidationMask(chosen)
    role="Held-out validation sensor";
end
sensor=struct('MeasuredRow',chosen,'SectionID',string(R.sectionMatch.SectionID(section)), ...
    'SectionName',string(R.sectionMatch.SectionName(section)),'Role',role, ...
    'MeanObservedPerHour',score(chosen));
start=double(T.StartMinute(q)); stop=double(T.EndMinute(q));
assert(all(stop>start),'Invalid stratum time boundaries.');
edges=unique([start;stop]);
for b=1:numel(edges)-1
    a=edges(b); z=edges(b+1);
    use=find(start<=a & stop>=z & valid(chosen,:).');
    if isempty(use)
        observed=NaN; predicted=NaN; weightSum=0;
    else
        w=counts(chosen,use).'./(stop(use)-start(use)); weightSum=sum(w);
        observed=sum(w.*obs(chosen,use).')/weightSum;
        predicted=sum(w.*pred(chosen,use).')/weightSum;
    end
    P(end+1,:)={chosen,a/60,z/60,observed,predicted,weightSum,numel(use)}; %#ok<AGROW>
end
end


function exportParmaFigure(f,pngFile,pdfFile,resolution)
% Export complete layouts, retaining every panel and shared colorbar.
requireFigure(f,'at export entry');
drawnow nocallbacks;
requireFigure(f,'after rendering');
layouts=findall(f,'Type','tiledlayout');
if ~isempty(layouts)
    target=layouts(1);
else
    ax=findall(f,'Type','axes');
    if numel(ax)==1, target=ax(1); else, target=[]; end
end
paths={char(pngFile),char(pdfFile)};
for formatIndex=1:2
    if ~isempty(target)
        try
            if formatIndex==1
                exportgraphics(target,paths{formatIndex}, ...
                    'Resolution',resolution,'BackgroundColor','white');
            else
                exportgraphics(target,paths{formatIndex}, ...
                    'ContentType','vector','BackgroundColor','white');
            end
            continue
        catch ME
            if ~isscalar(f)||~isgraphics(f,'figure'), rethrow(ME); end
            warning('ParmaFigures:ExportFallback', ...
                'Using whole-figure print for %s: %s',paths{formatIndex},ME.message);
        end
    end
    % Whole-figure fallback preserves all panels and shared colorbars.
    % Match PDF paper size to figure aspect ratio to avoid page clipping.
    requireFigure(f,'before the whole-figure export fallback');
    pixels=getappdata(f,'ParmaExportSizePixels');
    if isempty(pixels)
        oldUnits=get(f,'Units'); set(f,'Units','pixels');
        position=get(f,'Position'); set(f,'Units',oldUnits);
        pixels=position(3:4);
    end
    inches=pixels/get(groot,'ScreenPixelsPerInch');
    set(f,'PaperUnits','inches','PaperPositionMode','manual', ...
        'PaperPosition',[0 0 inches],'PaperSize',inches,'InvertHardcopy','off');
    if formatIndex==1
        print(f,paths{formatIndex},'-dpng',sprintf('-r%d',round(resolution)));
    else
        print(f,paths{formatIndex},'-dpdf','-painters');
    end
end
end


function requireFigure(f,stage)
if ~isscalar(f)||~isgraphics(f,'figure')
    error('ParmaFigures:InvalidFigureHandle', ...
        ['The figure handle is invalid %s. No export fallback was attempted. ' ...
         'Rerun the complete function with the saved results; ensure that no ' ...
         'timer or callback deletes figures while it is running.'],stage);
end
end


function removeFigureTitles(f)
% Publication figures: put descriptive titles and peak stratum details in
% LaTeX captions. Retain axis labels, ticks, colorbars, and legends.
a=findall(f,'Type','axes');
for k=1:numel(a)
    a(k).Title.String='';
    if isprop(a(k),'Subtitle'), a(k).Subtitle.String=''; end
end
layouts=findall(f,'Type','tiledlayout');
for k=1:numel(layouts)
    layouts(k).Title.String='';
    if isprop(layouts(k),'Subtitle'), layouts(k).Subtitle.String=''; end
end
end

function q=pctLocal(v,pr)
% Percentile of sorted vector v (linear interpolation, no toolbox needed).
n=numel(v); if n==1, q=v; return; end
h=1+(n-1)*pr/100; lo=floor(h); hi=ceil(h);
q=v(lo)+(h-lo)*(v(hi)-v(lo));
end
