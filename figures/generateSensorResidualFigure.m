function out = generateSensorResidualFigure(results,varargin)
%GENERATESENSORRESIDUALFIGURE Map of the per-sensor residuals of the NNLS fit.
%
%   out = generateSensorResidualFigure(results)
%   out = generateSensorResidualFigure('parma_tomography_results.mat','OutputDir','diagnostics')
%
%   Diagnostic, not used by the paper. Each sensor is a marker colored by its
%   relative residual (predicted - observed)/observed of the NNLS fit of the
%   sensor means, and sized by its mean count:
%     (a) dictionary-wide NNLS fit of the pooled means, all candidate classes;
%     (b) median over strata of the stratum NNLS fit with the fitted classes.
%   Blue: the paths predict too little traffic; red: too much. Isolated strong
%   colors usually indicate a wrong match or direction, a sensor near the map
%   boundary or a faulty detector; colored stretches along a road indicate
%   traffic entering or leaving between sensors. Circles are estimation
%   sensors, triangles validation sensors (predicted, not fitted). Sensors flagged
%   by diagnoseParmaSensors are ringed in black and labelled with their number.
%
%   Name/value: OutputDir ('diagnostics'), Clip (1: color scale -1..1),
%   Visible ('off'), Resolution (300). Writes sensor_residuals.pdf/.png/.fig.

p = inputParser;
addParameter(p,'OutputDir','diagnostics',@(x)ischar(x)||isstring(x));
addParameter(p,'Clip',1,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'Visible','off',@(x)any(strcmpi(x,{'on','off'})));
addParameter(p,'Resolution',300,@(x)isnumeric(x)&&isscalar(x)&&x>0);
parse(p,varargin{:}); o = p.Results;
% All output paths lie outside the library (see safeOutputPath).
o = guardOutputFields(o,{'OutputDir'});
if ischar(results) || (isstring(results) && isscalar(results))
    d = load(results,'results'); results = d.results;
end
assert(isfield(results,'sensorResiduals') && isfield(results,'G'), ...
    'generateSensorResidualFigure:OldResults', ...
    'Results lack sensorResiduals: rerun the modified runParmaTrafficTomography.');
T = results.sensorResiduals; G = results.G;
lon = double(G.Nodes.Lon); lat = double(G.Nodes.Lat);
en = G.Edges.EndNodes;
if isnumeric(en), s = double(en(:,1)); t = double(en(:,2));
else, s = findnode(G,en(:,1)); t = findnode(G,en(:,2)); end
xs = [lon(s) lon(t) nan(numel(s),1)].'; ys = [lat(s) lat(t) nan(numel(s),1)].';

f = figure('Color','w','Visible','off','Units','pixels','Position',[40 40 1500 700]);
tl = tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');
vals = {T.RelResidualDictionaryNNLS, T.MedianRelResidualStrataNNLS};
titles = {'(a) Dictionary-wide NNLS, pooled means','(b) Median over strata, stratum NNLS'};
cmap = divergingMap(256);
sz = 12+60*sqrt(max(T.MeanCount,0)/max([T.MeanCount;1]));
isVal = T.Role=="validation";
for k = 1:2
    ax = nexttile(tl); hold(ax,'on');
    plot(ax,xs(:),ys(:),'-','Color',[.82 .82 .82],'LineWidth',0.6);
    v = max(-o.Clip,min(o.Clip,vals{k}));
    ok = isfinite(v) & isfinite(T.Lon) & isfinite(T.Lat);
    scatter(ax,T.Lon(ok&~isVal),T.Lat(ok&~isVal),sz(ok&~isVal),v(ok&~isVal),'filled', ...
        'MarkerEdgeColor',[.25 .25 .25]);
    scatter(ax,T.Lon(ok&isVal),T.Lat(ok&isVal),sz(ok&isVal),v(ok&isVal),'^','filled', ...
        'MarkerEdgeColor','k','LineWidth',1);
    % Sensors flagged by diagnoseParmaSensors: black ring and measurement number
    if isfield(results,'sensorDiagnostics') && istable(results.sensorDiagnostics) && ...
            ismember('Flag',results.sensorDiagnostics.Properties.VariableNames)
        fl = strlength(results.sensorDiagnostics.Flag)>0 & isfinite(T.Lon) & isfinite(T.Lat);
        scatter(ax,T.Lon(fl),T.Lat(fl),sz(fl)+80,'o','MarkerEdgeColor','k','LineWidth',1.6);
        idx = find(fl);
        for k = idx(:).'
            text(ax,T.Lon(k),T.Lat(k),sprintf('  %d',k),'FontSize',8,'Color','k');
        end
    end
    colormap(ax,cmap); caxis(ax,[-o.Clip o.Clip]);
    axis(ax,'equal'); xlim(ax,[min(lon) max(lon)]); ylim(ax,[min(lat) max(lat)]);
    set(ax,'XTick',[],'YTick',[],'Box','on');
    title(ax,titles{k},'FontWeight','normal');
end
cb = colorbar(ax); cb.Layout.Tile = 'east';
cb.Label.String = sprintf('(predicted - observed)/observed, clipped to \\pm%g',o.Clip);
folder = char(o.OutputDir); if ~isfolder(folder), mkdir(folder); end
base = fullfile(folder,'sensor_residuals');
drawnow;
exportgraphics(tl,[base '.pdf'],'ContentType','vector','BackgroundColor','white');
exportgraphics(tl,[base '.png'],'Resolution',o.Resolution,'BackgroundColor','white');
savefig(f,[base '.fig']);
set(f,'Visible',o.Visible);
out = struct('figure',f,'files',string({[base '.pdf'],[base '.png'],[base '.fig']}));
fprintf('Saved sensor-residual map to %s\n',folder);
end

function c = divergingMap(n)
% Blue - white - red, symmetric around zero.
h = floor(n/2); u = linspace(0,1,h).';
blue = [0.13 0.40 0.75]; red = [0.80 0.15 0.15];
lower = (1-u)*blue + u*[1 1 1];
upper = (1-u)*[1 1 1] + u*red;
c = [lower; upper];
if size(c,1) < n, c = [c; red]; end
end
