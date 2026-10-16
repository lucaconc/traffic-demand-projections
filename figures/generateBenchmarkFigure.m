function out = generateBenchmarkFigure(R,varargin)
%GENERATEBENCHMARKFIGURE Error versus sample size for the identification benchmark.
%
%   out = generateBenchmarkFigure(R)                 % R from runIdentificationBenchmark
%   out = generateBenchmarkFigure('benchmark_results/benchmark_results.mat','RateScale',10)
%
%   One panel per routing matrix (easy, intermediate, difficult, complex): median relative
%   error over the paired replicates against N, both axes logarithmic, one line
%   per method, labelled at its right end (no legend). Methods whose identifying
%   system has rank below the number of paths are drawn dashed: their error
%   floor reflects missing information. Shaded bands are the interquartile ranges
%   of the projected composite likelihood.
%
%   Name/value: RateScale (middle of R.rateScales), OutputDir ('benchmark_figures'),
%   FileName ('benchmark_error'), Visible ('off'), FontSize (11).

p = inputParser;
addParameter(p,'RateScale',[],@(x)isempty(x)||isnumeric(x));
addParameter(p,'OutputDir','benchmark_figures',@(x)ischar(x)||isstring(x));
addParameter(p,'FileName','benchmark_error',@(x)ischar(x)||isstring(x));
addParameter(p,'Visible','off',@(x)any(strcmpi(x,{'on','off'})));
addParameter(p,'FontSize',11,@(x)isnumeric(x)&&isscalar(x));
parse(p,varargin{:}); o = p.Results;
% All output paths lie outside the library (see safeOutputPath).
o = guardOutputFields(o,{'OutputDir'});
if ischar(R) || (isstring(R) && isscalar(R)), d = load(R,'R'); R = d.R; end
scales = R.rateScales(:).';
if isempty(o.RateScale), sc = scales(ceil(numel(scales)/2)); else, sc = o.RateScale; end
S = R.summary;
colors = [0.45 0.45 0.45; 0.40 0.70 0.90; 0.20 0.50 0.80; 0.05 0.30 0.62; 0.02 0.15 0.40; ...
          0.55 0.30 0.70; 0.95 0.55 0.25; 0.80 0.15 0.05; 0.20 0.65 0.40; 0.85 0.60 0.05];
markers = {'s','^','v','d','<','p','o','*','h','x'};
short = {'NNLS','CUM_2','CUM_3','CUM_4','CUM_5','GLS','PCL','GMM','PAIR','CA-PCL'};

profiles = {'balanced'};
if isfield(R,'rateProfiles') && ismember('Profile',S.Properties.VariableNames), profiles = R.rateProfiles; end
profileLabel = containers.Map({'balanced','heavyPairs'},{'balanced rates','heavy two-sensor paths'});
f = figure('Color','w','Visible','off','Units','pixels','Position',[40 40 380*numel(R.cases)+120 420*numel(profiles)+50]);
try, theme(f,'light'); catch, end
tl = tiledlayout(f,numel(profiles),numel(R.cases),'TileSpacing','compact','Padding','compact');
for pr = 1:numel(profiles)
for c = 1:numel(R.cases)
    ax = nexttile(tl); hold(ax,'on');
    n = size(R.cases(c).A,2);
    yAll = []; lab = struct('x',[],'y',[],'t',{{}},'c',zeros(0,3));
    for m = 1:numel(R.methods)
        if ~R.applicable(c,m), continue; end
        sel = S.Case==string(R.cases(c).name) & S.RateScale==sc & S.Method==string(R.methods{m});
        if ismember('Profile',S.Properties.VariableNames), sel = sel & S.Profile==string(profiles{pr}); end
        T = sortrows(S(sel,:),'N');
        if isempty(T) || all(~isfinite(T.MedianError)), continue; end
        deficient = isfinite(R.rank(c,m)) && R.rank(c,m) < n;
        ls = '-'; if deficient, ls = '--'; end
        plot(ax,T.N,T.MedianError,ls,'Color',colors(m,:),'LineWidth',1.6, ...
            'Marker',markers{m},'MarkerFaceColor',colors(m,:),'MarkerSize',5);
        ok = isfinite(T.MedianError);
        last = find(ok,1,'last');
        lab.x(end+1) = T.N(last)*1.12; lab.y(end+1) = T.MedianError(last);
        lab.t{end+1} = short{m}; lab.c(end+1,:) = colors(m,:);
        yAll = [yAll; T.MedianError(ok)]; %#ok<AGROW>
    end
    set(ax,'XScale','log','YScale','log','FontSize',o.FontSize,'Color','w','XColor','k','YColor','k', ...
        'GridColor',[.15 .15 .15],'GridAlpha',.15,'Box','on');
    grid(ax,'on');
    xlim(ax,[min(R.sampleSizes)/1.3 max(R.sampleSizes)*2.6]);
    yAll = yAll(isfinite(yAll) & yAll>0);
    if ~isempty(yAll), ylim(ax,[min(yAll)/1.5 max(yAll)*1.5]); end
    xticks(ax,R.sampleSizes);
    xlabel(ax,'Number of observations N','Color','k');
    if c == 1
        pl = profiles{pr}; if isKey(profileLabel,pl), pl = profileLabel(pl); end
        ylabel(ax,{'Median relative error',['(' pl ')']},'Color','k');
    end
    if pr == 1, title(ax,R.cases(c).label,'FontWeight','normal','Color','k'); end
    placeLabels(ax,lab.x,lab.y,lab.t,lab.c,o.FontSize-2);
end
end
folder = char(o.OutputDir); if ~isfolder(folder), mkdir(folder); end
base = fullfile(folder,char(o.FileName));
drawnow;
exportgraphics(tl,[base '.pdf'],'ContentType','vector','BackgroundColor','white');
exportgraphics(tl,[base '.png'],'Resolution',300,'BackgroundColor','white');
savefig(f,[base '.fig']);
set(f,'Visible',o.Visible);
out = struct('figure',f,'rateScale',sc,'files',string({[base '.pdf'],[base '.png'],[base '.fig']}));
end

function placeLabels(ax,x,y,txt,col,fontSize)
% Line-end labels spread vertically (logarithmic axis) so that they do not
% overlap; a thin connector links a moved label to its line end.
if isempty(x), return; end
yl = get(ax,'YLim'); ly = log10(max(y,realmin)); span = diff(log10(yl));
gap = 0.055*span;
[ly0,ord] = sort(ly); lyA = ly0;
for i = 2:numel(lyA), lyA(i) = max(lyA(i),lyA(i-1)+gap); end
over = lyA(end) - (log10(yl(2)) - 0.02*span);
if over > 0, lyA = lyA - over; end
for i = numel(lyA)-1:-1:1, lyA(i) = min(lyA(i),lyA(i+1)-gap); end
for i = 1:numel(ord)
    k = ord(i); yy = 10^lyA(i);
    if abs(lyA(i)-ly0(i)) > 1e-9
        plot(ax,[x(k) x(k)],[y(k) yy],'-','Color',[col(k,:) 0.5],'LineWidth',0.5,'HandleVisibility','off');
    end
    text(ax,x(k),yy,txt{k},'Color',col(k,:),'FontSize',fontSize,'VerticalAlignment','middle');
end
end
