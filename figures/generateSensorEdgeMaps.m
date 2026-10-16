function out = generateSensorEdgeMaps(results,varargin)
%GENERATESENSOREDGEMAPS Local maps of sensors with dubious edge assignments.
%
%   out = generateSensorEdgeMaps(results)
%   out = generateSensorEdgeMaps('parma_tomography_results.mat','OutputDir','diagnostics/edge_maps')
%   out = generateSensorEdgeMaps(R,'ExtraNames',{'PRM12-S1 | GRAMSCI (LATO EST) >> OVEST /'})
%
%   Diagnostic, not used by the paper. One zoomed map per case:
%     * every group of indistinguishable sensors (identical path signatures,
%       results.indistinguishableSensors) whose counts differ (max/min above
%       MinCountRatio), that contains opposite directions (INGRESSO/USCITA,
%       NORD/SUD, EST/OVEST), or whose sensors share matched edges;
%     * every sensor flagged by diagnoseParmaSensors as possibly reversed,
%       swapped or inconsistent with the other sensors;
%     * any sensor named in ExtraNames.
%   Each map shows the road network around the sensors (grey), the directed
%   edges each sensor is matched to (thick colored lines with an arrow at the
%   end showing the direction of travel, drawn slightly to the right of the
%   road axis), the reverse edges of those edges where they exist (dashed),
%   and the sensor GPS position (circle with the sensor number). The title
%   gives the reason; the panel text lists, per sensor, its name, mean count,
%   NNLS prediction and matched edge IDs, and whether a reverse edge exists.
%
%   Name/value: OutputDir ('diagnostics/edge_maps'), RadiusMeters (350),
%   MinCountRatio (1.5), MaxMaps (60), ExtraNames ({}), Visible ('off').
%   Writes sensor_edge_maps.pdf (one page per map), one PNG per map and
%   sensor_edge_maps.csv (one row per mapped sensor).

p = inputParser;
addParameter(p,'OutputDir',fullfile('diagnostics','edge_maps'),@(x)ischar(x)||isstring(x));
addParameter(p,'RadiusMeters',350,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'MinCountRatio',1.5,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
addParameter(p,'MaxMaps',60,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
addParameter(p,'ExtraNames',{},@(x)iscell(x)||isstring(x));
addParameter(p,'Visible','off',@(x)any(strcmpi(x,{'on','off'})));
parse(p,varargin{:}); o = p.Results;
% All output paths lie outside the library (see safeOutputPath).
o = guardOutputFields(o,{'OutputDir'});
if ischar(results) || (isstring(results) && isscalar(results))
    d = load(results,'results'); results = d.results;
end
G = results.G; M = results.sectionMatch; sel = results.selectedSections(:);
names = strtrim(string(M.SectionName(sel)));
names = regexprep(names,'[\x00-\x1F\x7F\xA0]',' ');
lonN = double(G.Nodes.Lon); latN = double(G.Nodes.Lat);
en = G.Edges.EndNodes;
if ~isnumeric(en), en = [findnode(G,en(:,1)) findnode(G,en(:,2))]; end
en = double(en);
edgeSets = M.DirectedEdgeIDs(sel);
latS = nan(numel(sel),1); lonS = latS;
try, latS = double(M.Lat(sel)); lonS = double(M.Lon(sel)); catch, end
meanC = nan(numel(sel),1); predC = meanC; flag = strings(numel(sel),1);
if isfield(results,'sensorDiagnostics') && istable(results.sensorDiagnostics) && height(results.sensorDiagnostics)==numel(sel)
    D = results.sensorDiagnostics;
    meanC = double(D.MeanCount); predC = double(D.PredictedNNLS);
    if ismember('Flag',D.Properties.VariableNames), flag = string(D.Flag); flag(ismissing(flag)) = ""; end
end

% ---- cases to map -------------------------------------------------------
cases = struct('members',{},'reason',{});
if isfield(results,'indistinguishableSensors') && istable(results.indistinguishableSensors) ...
        && height(results.indistinguishableSensors) > 0
    T = results.indistinguishableSensors;
    for g = unique(T.Group).'
        mem = T.Measurement(T.Group==g); mem = mem(:).';
        mc = meanC(mem); ratio = max(mc)/max(min(mc),0.1);
        sharedEdges = numel(unique(cellfun(@(e)mat2str(sort(round(double(e(:)')))),edgeSets(mem),'UniformOutput',false))) < numel(mem);
        opposite = hasOppositeDirections(names(mem));
        if ratio > o.MinCountRatio || sharedEdges || opposite
            why = strings(0,1);
            if opposite, why(end+1) = "opposite directions share all paths"; end %#ok<AGROW>
            if sharedEdges, why(end+1) = "sensors matched to the same edges"; end %#ok<AGROW>
            if ratio > o.MinCountRatio, why(end+1) = sprintf("counts differ by a factor %.1f",ratio); end %#ok<AGROW>
            cases(end+1) = struct('members',mem,'reason',"Indistinguishable group: "+strjoin(why,"; ")); %#ok<AGROW>
        end
    end
end
mapped = [cases.members];
flagged = find(contains(flag,["reversed","swapped","inconsistent"]));
for k = flagged(:).'
    if ~ismember(k,mapped)
        cases(end+1) = struct('members',k,'reason',"Flagged by the sensor diagnostics: "+flag(k)); %#ok<AGROW>
    end
end
if ~isempty(o.ExtraNames)
    extra = find(ismember(lower(names),lower(strtrim(string(o.ExtraNames(:))))));
    for k = extra(:).'
        if ~ismember(k,[cases.members])
            cases(end+1) = struct('members',k,'reason',"Requested"); %#ok<AGROW>
        end
    end
end
cases = cases(1:min(numel(cases),o.MaxMaps));
folder = char(o.OutputDir); if ~isfolder(folder), mkdir(folder); end
pdfFile = fullfile(folder,'sensor_edge_maps.pdf');
if isfile(pdfFile), delete(pdfFile); end
rows = {};
colors = lines(8);

for c = 1:numel(cases)
    mem = cases(c).members;
    % Center and extent: sensor positions and matched edges
    px = []; py = [];
    for k = mem
        ids = validEdges(edgeSets{k},size(en,1));
        px = [px; lonN(en(ids,1)); lonN(en(ids,2)); lonS(k)]; %#ok<AGROW>
        py = [py; latN(en(ids,1)); latN(en(ids,2)); latS(k)]; %#ok<AGROW>
    end
    ok = isfinite(px) & isfinite(py);
    if ~any(ok), continue; end
    cx = mean(px(ok)); cy = mean(py(ok)); kx = cosd(cy);
    rDeg = o.RadiusMeters/111320;
    span = max([rDeg, (max(px(ok))-min(px(ok)))*kx/2*1.3, (max(py(ok))-min(py(ok)))/2*1.3]);
    xl = cx+[-1 1]*span/kx; yl = cy+[-1 1]*span;

    f = figure('Color','w','Visible','off','Units','pixels','Position',[40 40 1200 820]);
    try, theme(f,'light'); catch, end
    ax = axes(f,'Position',[0.05 0.30 0.62 0.62]); hold(ax,'on');
    % Background network within the window
    inWin = (lonN(en(:,1))>=xl(1)-span & lonN(en(:,1))<=xl(2)+span & latN(en(:,1))>=yl(1)-span & latN(en(:,1))<=yl(2)+span) | ...
            (lonN(en(:,2))>=xl(1)-span & lonN(en(:,2))<=xl(2)+span & latN(en(:,2))>=yl(1)-span & latN(en(:,2))<=yl(2)+span);
    e = find(inWin);
    xs = [lonN(en(e,1)) lonN(en(e,2)) nan(numel(e),1)].'; ys = [latN(en(e,1)) latN(en(e,2)) nan(numel(e),1)].';
    plot(ax,xs(:),ys(:),'-','Color',[.75 .75 .75],'LineWidth',1);
    off = 0.06*span;                                   % draw directions to the right of the axis
    % Edge IDs of the roads in the window (to the right of each direction), for
    % the correction table of runParmaTrafficTomography (EdgeOverrideCSV).
    mx = (lonN(en(e,1))+lonN(en(e,2)))/2; my = (latN(en(e,1))+latN(en(e,2)))/2;
    inner = find(mx>=xl(1) & mx<=xl(2) & my>=yl(1) & my<=yl(2));
    for t = inner(:).'
        [x1,y1,x2,y2] = shiftRight(lonN(en(e(t),1)),latN(en(e(t),1)),lonN(en(e(t),2)),latN(en(e(t),2)),1.6*off,kx);
        text(ax,(x1+x2)/2,(y1+y2)/2,sprintf('%d',e(t)),'FontSize',6,'Color',[.40 .40 .40], ...
            'HorizontalAlignment','center','Clipping','on');
    end
    txt = strings(0,1);
    for j = 1:numel(mem)
        k = mem(j); col = colors(mod(j-1,size(colors,1))+1,:);
        ids = validEdges(edgeSets{k},size(en,1));
        hasRev = false(numel(ids),1);
        for t = 1:numel(ids)
            a = en(ids(t),1); b = en(ids(t),2);
            [x1,y1,x2,y2] = shiftRight(lonN(a),latN(a),lonN(b),latN(b),off,kx);
            plot(ax,[x1 x2],[y1 y2],'-','Color',col,'LineWidth',4);
            drawArrowHead(ax,x1,y1,x2,y2,col,0.12*span,kx);
            r = find(en(:,1)==b & en(:,2)==a,1);
            if ~isempty(r)
                hasRev(t) = true;
                [u1,v1,u2,v2] = shiftRight(lonN(b),latN(b),lonN(a),latN(a),off,kx);
                plot(ax,[u1 u2],[v1 v2],'--','Color',col,'LineWidth',1.5);
            end
        end
        if isfinite(lonS(k))
            plot(ax,lonS(k),latS(k),'o','MarkerSize',11,'MarkerFaceColor',col,'MarkerEdgeColor','k','LineWidth',1);
            text(ax,lonS(k),latS(k),sprintf('  %d',k),'FontSize',11,'FontWeight','bold','Color','k');
        end
        revText = "no reverse edge";
        if all(hasRev), revText = "reverse edges exist"; elseif any(hasRev), revText = "reverse edge for some"; end
        txt(end+1) = sprintf('%d  %s\n     mean %.1f, NNLS %.1f; edges %s (%s)%s', k, names(k), meanC(k), predC(k), ...
            mat2str(ids(:).'), revText, ternary(flag(k)~="", "; flag: "+flag(k), "")); %#ok<AGROW>
        rows(end+1,:) = {c, cases(c).reason, k, names(k), meanC(k), predC(k), string(mat2str(ids(:).')), revText, flag(k)}; %#ok<AGROW>
    end
    axis(ax,'equal'); xlim(ax,xl); ylim(ax,yl);
    set(ax,'XTick',[],'YTick',[],'Box','on','Color','w','XColor','k','YColor','k');
    title(ax,sprintf('Map %d. %s',c,cases(c).reason),'FontWeight','normal','Color','k','Interpreter','none');
    % Scale bar (100 m)
    sb = 100/111320/kx;
    plot(ax,xl(1)+0.05*diff(xl)+[0 sb],yl(1)+0.05*diff(yl)*[1 1],'k-','LineWidth',2);
    text(ax,xl(1)+0.05*diff(xl),yl(1)+0.09*diff(yl),'100 m','Color','k','FontSize',9);
    % Text panel
    tax = axes(f,'Position',[0.05 0.02 0.9 0.24],'Visible','off');
    text(tax,0,1,strjoin(txt,newline),'VerticalAlignment','top','FontSize',10,'Interpreter','none','Color','k');
    lax = axes(f,'Position',[0.70 0.30 0.27 0.62],'Visible','off');
    text(lax,0,1,sprintf(['How to read this map\n\n' ...
        'Grey: road network.\nThick colored line: edge matched to the sensor of\nthe same color; the arrowhead shows the direction\nof travel (drawn right of the road axis).\n' ...
        'Dashed: reverse edge of a matched edge.\nCircle: sensor GPS position and number.\nSmall grey numbers: edge IDs, for the\ncorrection table (EdgeOverrideCSV).\n\n' ...
        'Opposite directions on one edge, or a missing\nreverse edge for a two-way road, point to a\nmatching or network problem.']), ...
        'VerticalAlignment','top','FontSize',10,'Color','k');
    base = fullfile(folder,sprintf('edge_map_%02d',c));
    exportgraphics(f,[base '.png'],'Resolution',150,'BackgroundColor','white');
    if c == 1, exportgraphics(f,pdfFile,'ContentType','vector','BackgroundColor','white');
    else, exportgraphics(f,pdfFile,'ContentType','vector','BackgroundColor','white','Append',true); end
    close(f);
end
if isempty(rows)
    T = table(); fprintf('No dubious sensor assignments to map.\n');
else
    T = cell2table(rows,'VariableNames',{'Map','Reason','Measurement','SectionName','MeanCount', ...
        'PredictedNNLS','MatchedEdges','ReverseEdges','Flag'});
    writetable(T,fullfile(folder,'sensor_edge_maps.csv'));
    fprintf('Saved %d local maps of dubious edge assignments to %s\n',numel(cases),pdfFile);
end
out = struct('maps',numel(cases),'table',T,'pdf',string(pdfFile));
end

% -------------------------------------------------------------------------
function ids = validEdges(e,nE)
ids = round(double(e(:))); ids = ids(isfinite(ids) & ids>=1 & ids<=nE);
end

function tf = hasOppositeDirections(nm)
% True if the names contain opposite direction words.
u = upper(nm);
pairs = {["INGRESSO","USCITA"],[">> NORD",">> SUD"],[">> EST",">> OVEST"]};
tf = false;
for k = 1:numel(pairs)
    if any(contains(u,pairs{k}(1))) && any(contains(u,pairs{k}(2))), tf = true; return; end
end
end

function [x1,y1,x2,y2] = shiftRight(x1,y1,x2,y2,off,kx)
% Shift a segment to the right of its direction of travel (isotropic offset).
dx = (x2-x1)*kx; dy = y2-y1; d = hypot(dx,dy); if d == 0, return; end
nx = dy/d; ny = -dx/d;
x1 = x1+off*nx/kx; x2 = x2+off*nx/kx; y1 = y1+off*ny; y2 = y2+off*ny;
end

function drawArrowHead(ax,x1,y1,x2,y2,col,len,kx)
% Filled triangle at (x2,y2) pointing along the segment.
dx = (x2-x1)*kx; dy = y2-y1; d = hypot(dx,dy); if d == 0, return; end
ux = dx/d; uy = dy/d; L = min(len,0.5*d); W = 0.45*L;
bx = x2*kx-L*ux; by = y2-L*uy;
px = [x2*kx, bx+W*uy, bx-W*uy]/kx; py = [y2, by-W*ux, by+W*ux];
patch(ax,px,py,col,'EdgeColor','k','LineWidth',0.5);
end

function v = ternary(c,a,b)
if c, v = a; else, v = b; end
end
