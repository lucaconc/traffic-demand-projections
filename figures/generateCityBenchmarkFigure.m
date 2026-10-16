function out = generateCityBenchmarkFigure(R,varargin)
%GENERATECITYBENCHMARKFIGURE Error and solution time on city-like networks.
%
%   out = generateCityBenchmarkFigure(R)      % R from runCityBenchmark (or its MAT file)
%
%   Two rows of four panels, one column per sweep (sensors, paths, observations,
%   contamination): median relative error of the class rates (top) and median
%   seconds per fit (bottom), both logarithmic, one line per method labelled at
%   its right end (no legend). A method is dashed at a level where it is rank
%   deficient in more than half of the replicates. The horizontal axis shows the
%   levels as categories; the baseline is marked in bold.
%
%   Name/value: OutputDir ('benchmark_figures'), FileName ('city_benchmark'),
%   Visible ('off'), FontSize (11).
p = inputParser;
addParameter(p,'OutputDir','benchmark_figures'); addParameter(p,'FileName','city_benchmark');
addParameter(p,'Visible','off'); addParameter(p,'FontSize',11);
parse(p,varargin{:}); o = p.Results;
% All output paths lie outside the library (see safeOutputPath).
o = guardOutputFields(o,{'OutputDir'});
if ischar(R) || (isstring(R) && isscalar(R)), d = load(R,'R'); R = d.R; end
S = R.summary; methods = R.methods; sweeps = R.sweeps; base = R.settings.Baseline;
% Style of each method by name (the run may contain any subset of them).
allNames = {'NNLS (means)','Cumulant matching','Bootstrap GLS','Projected CL', ...
            'Projected CL, robust','Projected GMM','Projected GMM, robust','Pairwise CL'};
allColors = [0.45 0.45 0.45; 0.20 0.50 0.80; 0.55 0.30 0.70; 0.95 0.55 0.25; 0.60 0.35 0.10; ...
             0.80 0.15 0.05; 0.45 0.05 0.05; 0.20 0.65 0.40];
allMarkers = {'s','d','p','o','^','*','x','h'};
allShort = {'NNLS','Cum.','GLS','Proj. CL','Proj. CL rob.','Proj. GMM','Proj. GMM rob.','Pair. CL'};
[tf,pos] = ismember(methods,allNames); pos(~tf) = 1;
colors = allColors(pos,:); markers = allMarkers(pos); short = allShort(pos);
xlab = {'Number of sensors','Number of paths','Number of observations N','Contaminated share of intervals'};
baseVal = [base.Sensors base.Paths base.N base.Contamination];
f = figure('Color','w','Visible','off','Units','pixels','Position',[40 40 1700 760]);
try, theme(f,'light'); catch, end
tl = tiledlayout(f,2,size(sweeps,1),'TileSpacing','compact','Padding','compact');
for row = 1:2
    for s = 1:size(sweeps,1)
        ax = nexttile(tl,(row-1)*size(sweeps,1)+s); hold(ax,'on');
        levels = sweeps{s,2}; L = numel(levels); lab = struct('x',[],'y',[],'t',{{}},'c',zeros(0,3));
        for k = 1:numel(methods)
            sel = S.Sweep==string(sweeps{s,1}) & S.Method==string(methods{k});
            T = sortrows(S(sel,:),'Level');
            if isempty(T), continue; end
            if row == 1, v = T.MedianError; else, v = T.MedianSeconds; end
            x = arrayfun(@(lv) find(levels==lv,1),T.Level);
            if all(~isfinite(v)), continue; end
            plot(ax,x,v,'-','Color',colors(k,:),'LineWidth',1.4,'Marker',markers{k}, ...
                'MarkerFaceColor',colors(k,:),'MarkerSize',5);
            def = T.RankDeficientShare > 0.5;             % dashed overlay where deficient
            if any(def) && row == 1
                plot(ax,x(def),v(def),'--','Color',[1 1 1],'LineWidth',0.8);
            end
            last = find(isfinite(v),1,'last');
            lab.x(end+1) = x(last)+0.08; lab.y(end+1) = v(last);
            lab.t{end+1} = short{k}; lab.c(end+1,:) = colors(k,:);
        end
        set(ax,'YScale','log','FontSize',o.FontSize,'Color','w','XColor','k','YColor','k', ...
            'GridColor',[.15 .15 .15],'GridAlpha',.15,'Box','on');
        grid(ax,'on'); xlim(ax,[0.7 L+0.9]); xticks(ax,1:L);
        tickLab = arrayfun(@(lv) sprintf('%g',lv),levels,'UniformOutput',false);
        b = find(levels==baseVal(s),1); if ~isempty(b), tickLab{b} = ['\bf' tickLab{b}]; end
        xticklabels(ax,tickLab);
        if row == 2, xlabel(ax,xlab{s},'Color','k'); end
        drawnow; placeLabels(ax,lab.x,lab.y,lab.t,lab.c,o.FontSize-3);
        if s == 1
            if row == 1, ylabel(ax,'Median relative error','Color','k');
            else, ylabel(ax,'Median seconds per fit','Color','k'); end
        end
    end
end
folder = char(o.OutputDir); if ~isfolder(folder), mkdir(folder); end
base_ = fullfile(folder,char(o.FileName)); drawnow;
exportgraphics(tl,[base_ '.pdf'],'ContentType','vector','BackgroundColor','white');
exportgraphics(tl,[base_ '.png'],'Resolution',300,'BackgroundColor','white');
savefig(f,[base_ '.fig']); set(f,'Visible',o.Visible);
out = struct('figure',f,'files',string({[base_ '.pdf'],[base_ '.png'],[base_ '.fig']}));
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
