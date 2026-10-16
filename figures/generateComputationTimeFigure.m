function out = generateComputationTimeFigure(synthetic,parma,varargin)
%GENERATECOMPUTATIONTIMEFIGURE Computation times of the estimator (paper Fig. 6).
%
%   out = generateComputationTimeFigure(syntheticResults,parmaResults)
%   out = generateComputationTimeFigure('synthetic_osm_tests/synthetic_results.mat', ...
%             'parma_tomography_results.mat','OutputDir','computation_figures')
%
%   Panel (a): synthetic study, sample-size sweep. Median time of projection
%              selection (design) and of the trust-region fit versus N, with
%              10-90% ranges over the Monte Carlo trials.
%   Panel (b): synthetic study, sensor sweep, versus the number of sensors m;
%              the number of observable classes is printed above each level.
%   Panel (c): Parma study, fit time of each reliably estimated stratum versus
%              its number of classes; the color gives the training intervals.
%
%   Times are wall-clock seconds measured by the estimator (info.designRuntimeSeconds,
%   info.runtimeSeconds) and by the Parma pipeline (estimationTimeByStratum). They
%   depend on the computer: results.machine records it when available.
%
%   Name/value: OutputDir ('computation_figures'), FontSize (16),
%   LabelFontSize (18), Resolution (300), Visible ('off').
%   Writes computation_times.pdf/.png/.fig and computation_times.csv.

p = inputParser;
addParameter(p,'OutputDir','computation_figures',@(x)ischar(x)||isstring(x));
addParameter(p,'FontSize',16,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'LabelFontSize',18,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'Resolution',300,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'Visible','off',@(x)any(strcmpi(x,{'on','off'})));
parse(p,varargin{:}); o = p.Results;
% All output paths lie outside the library (see safeOutputPath).
o = guardOutputFields(o,{'OutputDir'});
S = readResults(synthetic); P = readResults(parma);
T = S.raw;
need = {'Experiment','Level','N','DesignSeconds','FitSeconds','FittedClassCount'};
assert(all(ismember(need,T.Properties.VariableNames)), ...
    'generateComputationTimeFigure:OldResults', ...
    'The synthetic results lack timing columns: rerun the modified runSyntheticOSMTomographyTests.');
need = {'estimationTimeByStratum','reducedClassCountFitted'};
assert(all(isfield(P,need)),'generateComputationTimeFigure:OldResults', ...
    'The Parma results lack timing fields: rerun the modified runParmaTrafficTomography.');
ok = true(height(T),1);
if ismember('Failure',T.Properties.VariableNames)
    ok = ismissing(T.Failure) | strlength(string(T.Failure))==0;
end
T.Experiment = string(T.Experiment);

% Summaries by level
A = levelSummary(T(ok & T.Experiment=="observations",:),'N');
B = levelSummary(T(ok & T.Experiment=="sensors",:),'Level');

% Parma: reliable (finite-estimate) strata only
valid = isfield(P,'validFitByStratum') && ~isempty(P.validFitByStratum);
if valid, use = logical(P.validFitByStratum(:)); else, use = isfinite(P.estimationTimeByStratum(:)); end
tP = double(P.estimationTimeByStratum(:)); nP = double(P.reducedClassCountFitted(:));
if isfield(P,'nTrainByStratum'), NP = double(P.nTrainByStratum(:)); else, NP = nan(size(tP)); end
use = use & isfinite(tP) & isfinite(nP) & tP>0;

f = figure('Color','w','Visible','off','Units','pixels','Position',[40 40 1700 520]);
% Force a white, light-themed figure even when MATLAB uses a dark theme.
try, theme(f,'light'); catch, end
% Three panels of 8 tiles each and one narrow tile for the color bar of panel (c)
tl = tiledlayout(f,1,25,'TileSpacing','compact','Padding','compact');
panelTile = [1 9 17];
designColor = [.05 .38 .75]; fitColor = [.85 .33 .10];

ax = nexttile(tl,panelTile(1),[1 8]); hold(ax,'on');
drawSeries(ax,A.x,A.design,designColor,'o','Projection selection');
drawSeries(ax,A.x,A.fit,fitColor,'s','Trust-region fit');
set(ax,'XScale','log','YScale','log'); grid(ax,'on'); box(ax,'on');
xlabel(ax,'Number of observations N'); ylabel(ax,'Time [s]');
xticks(ax,A.x); xticklabels(ax,compose('%g',A.x));

ax = nexttile(tl,panelTile(2),[1 8]); hold(ax,'on');
drawSeries(ax,B.x,B.design,designColor,'o','Projection selection');
drawSeries(ax,B.x,B.fit,fitColor,'s','Trust-region fit');
set(ax,'YScale','log'); grid(ax,'on'); box(ax,'on');
xlabel(ax,'Number of sensors m'); ylabel(ax,'Time [s]');
xticks(ax,B.x);
yl = ylim(ax);
for k = 1:numel(B.x)
    text(ax,B.x(k),yl(2),sprintf('n=%g',B.classes(k)),'HorizontalAlignment','center', ...
        'VerticalAlignment','top','FontSize',max(8,o.FontSize-4),'Color','k');
end

ax = nexttile(tl,panelTile(3),[1 8]); hold(ax,'on');
if any(use)
    if any(isfinite(NP(use)))
        % Explicit RGB colors and a hand-drawn color bar: independent of the
        % figure/axes colormap and of MATLAB's color theme.
        c = NP(use);
        lo = min(c); hi = max(c); if hi <= lo, hi = lo+1; end
        cmap = sequentialMap(256);
        rgb = cmap(1+round(255*min(max((c-lo)/(hi-lo),0),1)),:);
        scatter(ax,nP(use),tP(use),40,rgb,'filled','MarkerEdgeColor',[.25 .25 .25]);
        drawColorBar(nexttile(tl,25),cmap,lo,hi,'Training intervals in the stratum',o.FontSize);
    else
        scatter(ax,nP(use),tP(use),36,fitColor,'filled','MarkerEdgeColor',[.2 .2 .2]);
    end
    set(ax,'YScale','log');
else
    text(ax,.5,.5,'No estimated strata','Units','normalized','HorizontalAlignment','center');
end
grid(ax,'on'); box(ax,'on');
xlabel(ax,'Number of classes n'); ylabel(ax,'Fit time per stratum [s]');

set(findall(f,'Type','axes'),'FontSize',o.FontSize,'XColor','k','YColor','k', ...
    'Color','w','GridColor',[.15 .15 .15],'GridAlpha',.15,'MinorGridColor',[.15 .15 .15]);
set(findall(f,'Type','text'),'Color','k');
labels = [findall(f,'Type','axes'); findall(f,'Type','colorbar')];
for k = 1:numel(labels)
    try, labels(k).Label.FontSize = o.LabelFontSize; catch, end
    try, labels(k).XLabel.FontSize = o.LabelFontSize; labels(k).YLabel.FontSize = o.LabelFontSize; catch, end
end

folder = char(o.OutputDir); if ~isfolder(folder), mkdir(folder); end
base = fullfile(folder,'computation_times');
drawnow;
exportgraphics(tl,[base '.pdf'],'ContentType','vector','BackgroundColor','white');
exportgraphics(tl,[base '.png'],'Resolution',o.Resolution,'BackgroundColor','white');
savefig(f,[base '.fig']);

csvRows = [table(repmat("synthetic N",numel(A.x),1),A.x(:),A.design(:,2),A.fit(:,2),nan(numel(A.x),1), ...
               'VariableNames',{'Panel','Level','MedianDesignSeconds','MedianFitSeconds','Classes'}); ...
           table(repmat("synthetic m",numel(B.x),1),B.x(:),B.design(:,2),B.fit(:,2),B.classes(:), ...
               'VariableNames',{'Panel','Level','MedianDesignSeconds','MedianFitSeconds','Classes'}); ...
           table(repmat("parma stratum",nnz(use),1),find(use),nan(nnz(use),1),tP(use),nP(use), ...
               'VariableNames',{'Panel','Level','MedianDesignSeconds','MedianFitSeconds','Classes'})];
writetable(csvRows,[base '.csv']);
set(f,'Visible',o.Visible);
out = struct('figure',f,'files',string({[base '.pdf'],[base '.png'],[base '.fig'],[base '.csv']}), ...
    'syntheticN',A,'syntheticSensors',B);
fprintf('Saved computation-time figure to %s\n',folder);
end

% -------------------------------------------------------------------------
function R = readResults(x)
if ischar(x) || (isstring(x) && isscalar(x))
    d = load(x,'results'); assert(isfield(d,'results'),'MAT file must contain results.');
    R = d.results;
else
    R = x;
end
end

function U = levelSummary(T,levelColumn)
% Median and 10/90 percentiles of design and fit times at each level.
x = unique(double(T.(levelColumn)));
U.x = x(:).';
U.design = nan(numel(x),3); U.fit = nan(numel(x),3); U.classes = nan(numel(x),1);
for k = 1:numel(x)
    r = double(T.(levelColumn))==x(k);
    U.design(k,:) = pct(double(T.DesignSeconds(r)));
    U.fit(k,:) = pct(double(T.FitSeconds(r)));
    c = double(T.FittedClassCount(r)); c = c(isfinite(c));
    if ~isempty(c), U.classes(k) = median(c); end
end
end

function q = pct(v)
v = sort(v(isfinite(v) & v>0)); q = nan(1,3);
if isempty(v), return; end
n = numel(v); pr = [.1 .5 .9];
for k = 1:3
    h = max(1,min(n,n*pr(k)+.5)); lo = floor(h); hi = ceil(h);
    q(k) = v(lo)+(h-lo)*(v(hi)-v(lo));
end
end

function drawSeries(ax,x,q,color,marker,name)
good = all(isfinite(q),2).';
if ~any(good), return; end
x = x(good); q = q(good,:);
patch(ax,[x fliplr(x)],[q(:,1).' fliplr(q(:,3).')],color,'FaceAlpha',.15, ...
    'EdgeColor','none','HandleVisibility','off');
plot(ax,x,q(:,2),['-' marker],'Color',color,'MarkerFaceColor',color, ...
    'LineWidth',2,'MarkerSize',7,'DisplayName',name);
end

function cmap = sequentialMap(n)
% Light yellow -> orange -> red -> dark purple; monotone in lightness, prints well.
anchors = [1.00 0.93 0.60; 0.99 0.70 0.30; 0.90 0.38 0.20; 0.65 0.12 0.30; 0.28 0.05 0.35];
t = linspace(0,1,size(anchors,1)); u = linspace(0,1,n);
cmap = [interp1(t,anchors(:,1),u).' interp1(t,anchors(:,2),u).' interp1(t,anchors(:,3),u).'];
end

function drawColorBar(cax,cmap,lo,hi,label,fontSize)
% Continuous color bar drawn as an RGB image in its own (narrow) axes.
img = reshape(cmap,[size(cmap,1) 1 3]);
image(cax,[0 1],[lo hi],img);
set(cax,'YDir','normal','XTick',[],'YAxisLocation','right','FontSize',fontSize, ...
    'XColor','k','YColor','k','Box','on','Color','w');
set(cax,'YTick',unique(round(linspace(lo,hi,5)))); ylim(cax,[lo hi]);
ylabel(cax,label,'Color','k');
end
