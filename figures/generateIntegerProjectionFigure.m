function file = generateIntegerProjectionFigure(R,outDir)
%GENERATEINTEGERPROJECTIONFIGURE Distance to the ML estimate and time against L.
%
%   file = generateIntegerProjectionFigure(R,outDir)
%
%   R from runIntegerProjectionStudy. Left: median distance of the projected
%   estimate to muML, in units of the ML error, against L (log scales); L*, the
%   largest observed count plus one (median over replicates), is filled.
%   Right: median seconds per fit against L, with the time of the exact ML fit
%   (convolution) as a dashed line of the same colour.
outDir = char(safeOutputPath(outDir)); if ~isfolder(outDir), mkdir(outDir); end
S = R.summary; names = R.cases;
cols = [0.20 0.45 0.75; 0.85 0.35 0.10; 0.30 0.60 0.30];
mk = {'o','s','^'};
fig = figure('Visible','off','Color','w','Units','centimeters','Position',[2 2 17 6.5]);
try, theme(fig,'light'); catch, end            % white background whatever the desktop theme
tl = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
axD = nexttile(tl); hold(axD,'on'); axT = nexttile(tl); hold(axT,'on');
h = gobjects(0); lab = {};
for c = 1:numel(names)
    T = S(strcmp(S.Case,names{c}),:);
    fin = isfinite(T.L); star = ~fin;
    x = T.LValue; dmed = max(T.MedianDistanceToML,1e-8);
    [xs,ord] = sort(x); ds = dmed(ord); ts = T.MedianSeconds(ord);
    h(end+1) = plot(axD,xs,ds,['-' mk{c}],'Color',cols(c,:),'LineWidth',1.2,'MarkerSize',5); %#ok<AGROW>
    lab{end+1} = names{c}; %#ok<AGROW>
    plot(axD,x(star),dmed(star),mk{c},'Color',cols(c,:),'MarkerFaceColor',cols(c,:),'MarkerSize',7);
    plot(axT,xs,ts,['-' mk{c}],'Color',cols(c,:),'LineWidth',1.2,'MarkerSize',5);
    plot(axT,x(star),T.MedianSeconds(star),mk{c},'Color',cols(c,:),'MarkerFaceColor',cols(c,:),'MarkerSize',7);
    yline(axT,T.MedianSecondsML(1),'--','Color',cols(c,:),'LineWidth',1);
end
set(axD,'XScale','log','YScale','log','FontSize',9,'Box','on','Color','w','XColor','k','YColor','k', ...
    'GridColor',[.15 .15 .15],'GridAlpha',.15); grid(axD,'on');
xlabel(axD,'L'); ylabel(axD,'Distance to \mu_{ML} / ML error');
title(axD,'(a) Convergence to maximum likelihood','FontWeight','normal','FontSize',9);
legend(axD,h,lab,'Location','southwest','Box','off','FontSize',8);
set(axT,'XScale','log','YScale','log','FontSize',9,'Box','on','Color','w','XColor','k','YColor','k', ...
    'GridColor',[.15 .15 .15],'GridAlpha',.15); grid(axT,'on');
xlabel(axT,'L'); ylabel(axT,'Median seconds per fit');
title(axT,'(b) Time (dashed: exact likelihood)','FontWeight','normal','FontSize',9);
file = fullfile(outDir,'integer_projection.pdf');
exportgraphics(fig,file,'ContentType','vector','BackgroundColor','white');
close(fig);
end
