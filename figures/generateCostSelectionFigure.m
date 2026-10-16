function file = generateCostSelectionFigure(files,outDir)
%GENERATECOSTSELECTIONFIGURE Error against time of the cost-aware selection.
%
%   file = generateCostSelectionFigure({corridorMat, syntheticMat}, outDir)
%
%   One panel per design whose results file exists (runCostSelectionStudy):
%   median relative error against median seconds per fit, both logarithmic,
%   with thin bars for the interquartile range of the error. CA-PCL designs are
%   joined in order of decreasing alpha and labelled by it; PCL (default design)
%   and GMM are reference points. Drawn at its final size for the full text
%   width (16 x 6.2 cm, 9 pt text), with a legend above the panels and axis
%   limits that leave room for every point, bar and label.
outDir = char(safeOutputPath(outDir)); if ~isfolder(outDir), mkdir(outDir); end
files = files(cellfun(@(f) ~isempty(f) && isfile(f), files));
file = '';
if isempty(files), return; end
names = struct('corridor','Corridor','synthetic','Realistic network');
cCost = [0.85 0.45 0.05]; cPCL = [0.10 0.10 0.10]; cGMM = [0.75 0.10 0.10];
light = @(c) c + 0.55*(1-c);                       % lighter shade for the bars
fs = 9;
fig = figure('Visible','off','Color','w','Units','centimeters','Position',[2 2 16 6.2]);
try, theme(fig,'light'); catch, end
tl = tiledlayout(fig,1,numel(files),'TileSpacing','compact','Padding','compact');
hLeg = gobjects(1,3);
for f = 1:numel(files)
    S = load(files{f},'R'); R = S.R;
    nM = numel(R.methods); med = nan(nM,1); lo = med; hi = med; tm = med;
    for k = 1:nM
        med(k) = q(R.err(:,k),0.5); lo(k) = q(R.err(:,k),0.25); hi(k) = q(R.err(:,k),0.75);
        tm(k) = q(R.secs(:,k),0.5);
    end
    isCost = startsWith(R.methods,'alpha=');
    a = cellfun(@(s) sscanf(s,'alpha=%g'), R.methods(isCost));
    [a,ord] = sort(a,'descend'); idx = find(isCost); idx = idx(ord);
    iD = find(strcmp(R.methods,'default')); iG = find(strcmp(R.methods,'GMM'));
    ax = nexttile(tl); hold(ax,'on');
    % interquartile ranges first, so the markers sit on top
    all_ = [idx(:); iD(:); iG(:)];
    cols = [repmat(cCost,numel(idx),1); repmat(cPCL,numel(iD),1); repmat(cGMM,numel(iG),1)];
    for k = 1:numel(all_)
        j = all_(k);
        if isfinite(lo(j)) && isfinite(hi(j))
            plot(ax,[tm(j) tm(j)],[lo(j) hi(j)],'-','Color',light(cols(k,:)),'LineWidth',1.6,'HandleVisibility','off');
        end
    end
    h1 = plot(ax,tm(idx),med(idx),'-o','Color',cCost,'MarkerFaceColor',cCost,'MarkerEdgeColor','w', ...
        'MarkerSize',6.5,'LineWidth',1.3);
    % One label per group of (nearly) coincident points, e.g. 'alpha=1, 0.5'
    done = false(numel(idx),1);
    for k = 1:numel(idx)
        if done(k), continue; end
        near = ~done & abs(log(tm(idx)./tm(idx(k)))) < 0.25 & abs(log(med(idx)./med(idx(k)))) < 0.12;
        done = done | near;
        txt = ['\\alpha=' strjoin(arrayfun(@(v) sprintf('%g',v),a(near),'UniformOutput',false),', ')];
        text(ax,max(tm(idx(near)))*1.15,mean(med(idx(near))),txt,'FontSize',fs-1,'Color',cCost, ...
            'VerticalAlignment','middle','HorizontalAlignment','left');
    end
    h2 = gobjects(1); h3 = gobjects(1);
    if ~isempty(iD), h2 = plot(ax,tm(iD),med(iD),'s','MarkerFaceColor',cPCL,'MarkerEdgeColor','w','MarkerSize',7.5); end
    if ~isempty(iG), h3 = plot(ax,tm(iG),med(iG),'d','MarkerFaceColor',cGMM,'MarkerEdgeColor','w','MarkerSize',7.5); end
    if f == 1, hLeg = [h1 h2 h3]; end
    % limits with margins: bars included, room on the right for the labels
    xs = tm(all_); xs = xs(isfinite(xs) & xs > 0);
    ys = [lo(all_); hi(all_); med(all_)]; ys = ys(isfinite(ys) & ys > 0);
    xl = [min(xs)/1.7, max(xs)*2.6]; yl = [min(ys)/1.2, max(ys)*1.2];
    set(ax,'XScale','log','YScale','log','XLim',xl,'YLim',yl,'FontSize',fs,'Box','on','Color','w', ...
        'XColor','k','YColor','k','XMinorTick','off','YMinorTick','off','XMinorGrid','off','YMinorGrid','off', ...
        'GridColor',[0.15 0.15 0.15],'GridAlpha',0.12,'TickLength',[0.012 0.012],'LineWidth',0.6);
    grid(ax,'on');
    xt = [0.01 0.02 0.05 0.1 0.2 0.5 1 2 5 10 20 50 100 200]; xt = xt(xt >= xl(1) & xt <= xl(2));
    yt = [0.02 0.03 0.05 0.07 0.1 0.15 0.2 0.3 0.5 0.7 1]; yt = yt(yt >= yl(1) & yt <= yl(2));
    set(ax,'XTick',xt,'XTickLabel',arrayfun(@(v) sprintf('%g',v),xt,'UniformOutput',false), ...
        'YTick',yt,'YTickLabel',arrayfun(@(v) sprintf('%g',v),yt,'UniformOutput',false));
    xlabel(ax,'Median time per fit [s]','FontSize',fs);
    if f == 1, ylabel(ax,'Median relative error','FontSize',fs); end
    d = lower(char(R.design)); if isfield(names,d), nm = names.(d); else, nm = d; end
    title(ax,sprintf('(%c) %s',char('a'+f-1),nm),'FontWeight','normal','FontSize',fs);
end
ok = isgraphics(hLeg);
labels = {'CA-PCL (labelled by \alpha)','PCL','GMM'};
lg = legend(hLeg(ok),labels(ok),'Orientation','horizontal','Box','off','FontSize',fs);
try
    lg.Layout.Tile = 'north';                      % above both panels (R2020b and later)
catch
    lg.Location = 'northoutside';                  % older releases: above the first panel
end
file = fullfile(outDir,'cost_selection.pdf');
if exist('exportgraphics','file')
    exportgraphics(fig,file,'ContentType','vector','BackgroundColor','white');
else                                               % releases before R2020a
    set(fig,'PaperUnits','centimeters','PaperSize',[16 6.2],'PaperPosition',[0 0 16 6.2]);
    print(fig,file,'-dpdf','-painters');
end
fprintf('Cost-selection figure written: %s\n',file);
close(fig);
end

function v = q(x,p)
% Quantile by linear interpolation, ignoring non-finite values.
x = sort(x(isfinite(x))); n = numel(x);
if n == 0, v = NaN; return; end
h = 1+(n-1)*p; a = floor(h); b = ceil(h); v = x(a)+(h-a)*(x(b)-x(a));
end
