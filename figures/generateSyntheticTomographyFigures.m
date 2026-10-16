function out=generateSyntheticTomographyFigures(R,varargin)
%GENERATESYNTHETICTOMOGRAPHYFIGURES Four paper figures from saved simulations.
% out=generateSyntheticTomographyFigures(synthetic)
% out=generateSyntheticTomographyFigures('synthetic_osm_tests/synthetic_results.mat')
% Uses runSyntheticOSMTomographyTests output. No estimator, internet, OSM
% download, or Statistics Toolbox is needed. MATLAB R2020b+ recommended.
% Name/value: OutputDir ('synthetic_paper_figures'), FontSize (22),
% LabelFontSize (26), Resolution (300), Visible ('on'), ShowCounts (false),
% ShowNNLS (false): overlay the median error of the NNLS fit of the means
% (red diamonds) when the results contain the NNLSRateError column.
% ShowMethods (true): overlay, at each level, the median error of every other
% method in the results (columns Err_<tag>: NNLS, cumulant matching, GLS,
% projected GMM, pairwise CL), with the markers of the benchmark figures; the
% violins remain those of the projected composite likelihood. A legend is drawn
% in the first panel.
% Violins (option, default off) show the projected CL's trial-error density;
% by default it is drawn, like every other method, as a line of medians.
% Violin shapes show trial-error density at each condition (equal x spacing).
% Grey points are individual trials, black bars medians, white dots means.
% Failures are excluded and counted; no density is invented for constant data.
% The summary CSV retains P10/P90 as additional numerical diagnostics.
% Saves four plots (three for legacy results) as PDF/PNG/FIG, summary.csv, and a LaTeX settings
% fragment used by synthetic_tomography_section.tex. No titles are drawn.

if nargin==0, R=fullfile('synthetic_osm_tests','synthetic_results.mat'); end
p=inputParser;
addParameter(p,'OutputDir','synthetic_paper_figures',@(x)ischar(x)||isstring(x));
addParameter(p,'FontSize',22,@positiveScalar);
addParameter(p,'LabelFontSize',26,@positiveScalar);
addParameter(p,'Resolution',300,@positiveScalar);
addParameter(p,'Visible','on',@(x)any(strcmpi(x,{'on','off'})));
addParameter(p,'ShowCounts',false,@(x)islogical(x)&&isscalar(x));
addParameter(p,'Violins',false,@(x)islogical(x)&&isscalar(x));   % false: projected CL as a median line
addParameter(p,'LogErrors',true,@(x)islogical(x)&&isscalar(x));  % logarithmic error axis
addParameter(p,'ShowNNLS',false,@(x)islogical(x)&&isscalar(x));
addParameter(p,'ShowMethods',true,@(x)islogical(x)&&isscalar(x));
parse(p,varargin{:}); o=p.Results;
% All output paths lie outside the library (see safeOutputPath).
o=guardOutputFields(o,{'OutputDir'});
if ischar(R)||isstring(R)
    data=load(R);
    assert(isfield(data,'results'),'MAT file must contain variable results.');
    R=data.results;
end
assert(isstruct(R)&&isfield(R,'raw')&&istable(R.raw),'Input must contain raw simulation trials.');
T=R.raw;
required={'Experiment','Level','Replicate','RelativeRateError'};
assert(all(ismember(required,T.Properties.VariableNames)),'Raw trial columns are missing.');
T.Experiment=string(T.Experiment);
% Levels may be +Inf (travel speed v=Inf: no travel time); only NaN is invalid.
assert(all(~isnan(T.Level))&&all(isfinite(T.Replicate)),'Invalid trial levels or replicate IDs.');
assert(all(T.RelativeRateError(isfinite(T.RelativeRateError))>=0),'Negative rate error found.');
assert(height(unique(T(:,{'Experiment','Level','Replicate'})))==height(T), ...
    'Duplicate experiment/level/replicate entries: use unique sweep levels.');
if isfield(R,'completedReplicates')&&isfield(R,'settings')&&isfield(R.settings,'Replicates')
    if R.completedReplicates<R.settings.Replicates
        warning('SyntheticFigures:PartialRun','Plotting %d of %d requested replicates from a checkpoint.', ...
            R.completedReplicates,R.settings.Replicates);
    end
end
names=["observations","outliers","impurity","sensors","traveltime"];
assert(all(ismember(names(1:3),unique(T.Experiment))),'Results must include the three original experiments.');
if ~any(T.Experiment=="sensors")
    warning('SyntheticFigures:MissingSensorSweep','Legacy results contain no sensor sweep; rerun the updated simulation to create it.');
    names=names(1:3);
end
if ~any(T.Experiment=="traveltime"), names(names=="traveltime")=[]; end
S=table('Size',[0 7],'VariableTypes',{'string','double','double','double','double','double','double'}, ...
    'VariableNames',{'Experiment','Level','MedianError','P10','P90','Successful','Failed'});
for name=names
    levels=unique(T.Level(T.Experiment==name));
    for level=levels.'
        use=T.Experiment==name&T.Level==level; values=double(T.RelativeRateError(use));
        valid=isfinite(values);
        if ismember('Failure',T.Properties.VariableNames)
            message=string(T.Failure(use));
            valid=valid & (ismissing(message)|strlength(message)==0);
        end
        v=sort(values(valid)); q=localPercentile(v,[.1 .5 .9]);
        S(end+1,:)={name,level,q(2),q(1),q(3),numel(v),nnz(~valid)}; %#ok<AGROW>
    end
end
S.MedianObservableFlowFraction=nan(height(S),1);
S.MedianExpandedClassRateError=nan(height(S),1);
S.MedianFittedClassCount=nan(height(S),1);
S.MedianActiveRateError=nan(height(S),1);
S.MedianInactiveFlowFraction=nan(height(S),1);
S.MedianNNLSRateError=nan(height(S),1);
for k=1:height(S)
    use=T.Experiment==S.Experiment(k)&T.Level==S.Level(k);
    source={'ObservableFlowFraction','ExpandedClassRateError','FittedClassCount', ...
        'ActiveRateError','InactiveFlowFraction','NNLSRateError'};
    target={'MedianObservableFlowFraction','MedianExpandedClassRateError','MedianFittedClassCount', ...
        'MedianActiveRateError','MedianInactiveFlowFraction','MedianNNLSRateError'};
    for j=1:numel(source)
        if ismember(source{j},T.Properties.VariableNames)
            v=double(T.(source{j})(use)); v=v(isfinite(v));
            if ~isempty(v), S.(target{j})(k)=median(v); end
        end
    end
end
if ~isfolder(o.OutputDir), mkdir(o.OutputDir); end
figures=gobjects(numel(names),1); files=strings(numel(names),3);
labels={'Number of observations','Number of contaminated observations','Overdispersion fraction \rho','Number of sensors','Travel speed [km/h]'};
stems={'01_error_observations','02_error_outliers','03_error_impurity','04_error_sensors','05_error_traveltime'};
for j=1:numel(names)
    U=S(S.Experiment==names(j),:);
    f=figure('Color','w','Visible','off','WindowStyle','normal', ...
        'Units','pixels','Position',[60 60 950 700]); figures(j)=f;
    try, theme(f,'light'); catch, end
    ax=axes(f); hold(ax,'on');
    tops=nan(height(U),1); medPCL=nan(height(U),1);
    for k=1:height(U)
        use=T.Experiment==names(j)&T.Level==U.Level(k);
        v=double(T.RelativeRateError(use)); valid=isfinite(v);
        if ismember('Failure',T.Properties.VariableNames)
            message=string(T.Failure(use));
            valid=valid & (ismissing(message)|strlength(message)==0);
        end
        if o.Violins
            tops(k)=drawErrorViolin(ax,k,v(valid));
        elseif any(valid)
            medPCL(k)=median(v(valid)); tops(k)=medPCL(k);
        end
    end
    legH=gobjects(0); legL={};
    if ~o.Violins
        legH(end+1)=plot(ax,1:height(U),medPCL,'-o','Color','k','MarkerFaceColor','k','MarkerSize',8,'LineWidth',1.8);
        legL{end+1}='PCL';
    end
    specs={'NNLS','NNLS',[0.45 0.45 0.45],'s'; 'Cum','CUM',[0.20 0.50 0.80],'d'; ...
           'GLS','GLS',[0.55 0.30 0.70],'p'; 'GMM','GMM',[0.80 0.15 0.05],'*'; ...
           'Pair','PAIR',[0.20 0.65 0.40],'h'; 'CostCL','CA-PCL',[0.85 0.60 0.05],'x'; ...
           'PriorGLS','GLS-P',[0.55 0.35 0.15],'v'};
    hasMethods=o.ShowMethods && any(startsWith(T.Properties.VariableNames,'Err_'));
    if hasMethods
        if o.Violins, legH(end+1)=plot(ax,NaN,NaN,'k-','LineWidth',2.1); legL{end+1}='PCL'; end
        for q=1:size(specs,1)
            col=['Err_' specs{q,1}];
            if ~ismember(col,T.Properties.VariableNames), continue; end
            med=nan(height(U),1);
            for k=1:height(U)
                use=T.Experiment==names(j)&T.Level==U.Level(k);
                v=double(T.(col)(use)); v=v(isfinite(v));
                if ~isempty(v), med(k)=median(v); end
            end
            if any(isfinite(med))
                legH(end+1)=plot(ax,1:height(U),med,['-' specs{q,4}],'Color',specs{q,3}, ...
                    'MarkerFaceColor',specs{q,3},'MarkerSize',9,'LineWidth',1.2); %#ok<AGROW>
                legL{end+1}=specs{q,2}; %#ok<AGROW>
                tops=max(tops,med);
            end
        end
        if j==1 && numel(legH)>1
            legend(ax,legH,legL,'Location','northeast','FontSize',max(8,o.FontSize-6),'Box','off');
        end
    end
    if ~hasMethods && o.ShowNNLS && ismember('MedianNNLSRateError',S.Properties.VariableNames)
        nn=U.MedianNNLSRateError;
        if any(isfinite(nn))
            plot(ax,1:height(U),nn,'d','Color',[.8 .1 .1],'MarkerFaceColor',[.8 .1 .1],'MarkerSize',9);
            tops=max(tops,nn);
        end
    end
    % Numeric labels identify levels; equally spaced positions avoid making
    % narrow violins unreadable on nonuniform experimental grids.
    xticks(ax,1:height(U)); tl=compose('%g',U.Level); tl(isinf(U.Level))={'\infty'}; xticklabels(ax,tl);
    xlim(ax,[.4 height(U)+.6]);
    ymax=max(tops(isfinite(tops))); if isempty(ymax)||ymax<=0, ymax=1; end
    if o.LogErrors
        % Limits from every plotted median (all methods), with a factor 1.5 margin
        yd=get(findobj(ax,'Type','line'),'YData'); if iscell(yd), yd=[yd{:}]; end
        yd=yd(isfinite(yd)&yd>0); if isempty(yd), yd=[1e-3 1]; end
        set(ax,'YScale','log'); ylim(ax,[min(yd)/1.5 max(yd)*1.5]);
    else
        ylim(ax,[0 1.2*ymax]);
    end
    if o.ShowCounts
        for k=1:height(U)
            top=tops(k); if ~isfinite(top), top=0; end
            if o.LogErrors, ty=max(top,realmin)*1.15; else, ty=top+.035*ymax; end
            text(ax,k,ty, ...
                sprintf('%d/%d',U.Successful(k),U.Successful(k)+U.Failed(k)), ...
                'FontSize',max(8,o.FontSize-2),'Color','k','HorizontalAlignment','center');
        end
    end
    if ~any(isfinite(tops))
        text(ax,.5,.5,'No successful fits','Units','normalized', ...
            'FontSize',o.FontSize,'Color','k','HorizontalAlignment','center');
    end
    xlabel(ax,labels{j},'FontSize',o.LabelFontSize,'Color','k');
    ylabel(ax,'Relative class-intensity error','FontSize',o.LabelFontSize,'Color','k');
    set(ax,'FontSize',o.FontSize,'XColor','k','YColor','k','Color','w');
    set([ax.XLabel ax.YLabel],'FontSize',o.LabelFontSize);
    grid(ax,'on'); box(ax,'on');
    base=fullfile(o.OutputDir,stems{j});
    files(j,:)=[string([base '.pdf']),string([base '.png']),string([base '.fig'])];
    exportFigure(f,ax,base,o.Resolution);
end
writetable(S,fullfile(o.OutputDir,'summary.csv'));
writeSettings(R,T,fullfile(o.OutputDir,'simulation_settings.tex'));
for j=1:numel(figures), set(figures(j),'Visible',o.Visible); end
out=struct('figures',figures,'files',files,'summary',S,'options',o);
fprintf('Saved %d figures and simulation_settings.tex to %s\n',numel(figures),char(o.OutputDir));
end

function q=localPercentile(v,prob)
% Linear interpolation at sample plotting positions (i-0.5)/n, clamped.
% Matches the conventional MATLAB percentile rule without a toolbox call.
q=nan(size(prob)); n=numel(v); if n==0, return; end
for k=1:numel(prob)
    h=max(1,min(n,n*prob(k)+.5)); lo=floor(h); hi=ceil(h);
    q(k)=v(lo)+(h-lo)*(v(hi)-v(lo));
end
end

function exportFigure(f,ax,base,dpi)
assert(isscalar(f)&&isgraphics(f,'figure'),'Figure was deleted before export.');
drawnow nocallbacks;
try
    exportgraphics(ax,[base '.pdf'],'ContentType','vector','BackgroundColor','white');
    exportgraphics(ax,[base '.png'],'Resolution',dpi,'BackgroundColor','white');
catch ME
    assert(isgraphics(f,'figure'),'Figure was deleted during export.');
    warning('SyntheticFigures:ExportFallback','Using whole-figure print: %s',ME.message);
    set(f,'Units','inches'); pos=get(f,'Position');
    set(f,'PaperUnits','inches','PaperPosition',[0 0 pos(3:4)], ...
        'PaperSize',pos(3:4),'InvertHardcopy','off');
    print(f,[base '.pdf'],'-dpdf','-painters');
    print(f,[base '.png'],'-dpng',sprintf('-r%d',round(dpi)));
end
savefig(f,[base '.fig']);
end

function writeSettings(R,T,file)
% Definitions are filled from actual output, never assumed defaults.
fid=fopen(file,'w'); assert(fid>=0,'Cannot write LaTeX settings.');
c=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%% Generated from saved simulation results.\n');
keys={'NumTruePaths','NumCandidatePaths','BaseObservations','OutlierMultiplier','ImpurityCV','ObservationMinutes'};
macros={'SimTruePaths','SimCandidatePaths','SimBaseN','SimOutlierScale','SimImpurityCV','SimIntervalMinutes'};
if isfield(R,'settings')
    for j=1:numel(keys)
        if isfield(R.settings,keys{j})
            define(fid,macros{j},sprintf('%g',R.settings.(keys{j})));
        end
    end
    if isfield(R.settings,'RateRange')
        define(fid,'SimRateLow',sprintf('%g',R.settings.RateRange(1)));
        define(fid,'SimRateHigh',sprintf('%g',R.settings.RateRange(2)));
    end
end
if isfield(R,'A')
    define(fid,'SimSensors',sprintf('%d',size(R.A,1)));
    define(fid,'SimClasses',sprintf('%d',size(R.A,2)));
end
define(fid,'SimReplicates',sprintf('%d',numel(unique(T.Replicate))));
for pair={"observations","SimSampleGrid";"outliers","SimOutlierGrid";"impurity","SimImpurityGrid";"sensors","SimSensorGrid";"traveltime","SimTravelGrid"}.'
    levels=unique(T.Level(T.Experiment==pair{1}));
    if isempty(levels), continue; end
    define(fid,char(pair{2}),strrep(strjoin(cellstr(compose('%g',levels)),', '),'Inf','\infty'));
end
end

function define(fid,name,value)
% csname avoids assumptions about whether a parent document predefines macros.
fprintf(fid,'\\expandafter\\def\\csname %s\\endcsname{%s}\n',name,value);
end
function tf=positiveScalar(x)
tf=isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0;
end


function top=drawErrorViolin(ax,position,v)
% Gaussian KDE with reflection at zero, since relative errors are nonnegative.
% Width normalized separately per condition; it represents density, not n.
% No Statistics Toolbox or release-specific violinplot function is required.
v=sort(v(:)); n=numel(v); top=NaN;
if n==0, return; end
color=[.05 .38 .75];
if n>=2 && max(v)-min(v)>eps(max(1,max(v)))
    quartiles=localPercentile(v,[.25 .75]);
    robustScale=(quartiles(2)-quartiles(1))/1.349;
    scale=std(v);
    if robustScale>0, scale=min(scale,robustScale); end
    h=max(.9*scale*n^(-.2),eps(max(1,max(v))));
    gridY=linspace(max(0,min(v)-3*h),max(v)+3*h,240).';
    density=zeros(size(gridY));
    for first=1:500:n
        values=v(first:min(first+499,n)).';
        density=density+sum(exp(-.5*((gridY-values)/h).^2) ...
            +exp(-.5*((gridY+values)/h).^2),2);
    end
    width=.35*density/max(density);
    patch(ax,[position-width;flipud(position+width)], ...
        [gridY;flipud(gridY)],color,'FaceAlpha',.35, ...
        'EdgeColor',color,'LineWidth',1.2);
    top=max(gridY);
else
    plot(ax,position+[-.25 .25],[v(1) v(1)],'-','Color',color,'LineWidth',3);
    top=max(v);
end
% Deterministic jitter shows every finite trial without changing the RNG.
jitter=.11*sin((1:n).'*2.399963229728653);
scatter(ax,position+jitter,v,15,[.4 .4 .4],'filled');
plot(ax,position+[-.15 .15],[median(v) median(v)],'k-','LineWidth',2.1);
plot(ax,position,mean(v),'ko','MarkerFaceColor','w','MarkerSize',6,'LineWidth',1.2);
end
