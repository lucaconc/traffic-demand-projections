function T = listEdgeCandidates(results,names,varargin)
%LISTEDGECANDIDATES Road edges around sensors, for the EdgeOverrideCSV table.
%
%   T = listEdgeCandidates(results,names)
%   T = listEdgeCandidates('parma_tomography_results.mat','parma_edge_overrides.csv', ...
%                          'OutputFile','parma_edge_candidates.csv')
%
%   For every sensor named in NAMES (a string array, or a CSV with a SectionName
%   column such as parma_edge_overrides.csv), lists the directed road edges whose
%   midpoint lies within Radius metres of the sensor's GPS position: edge ID,
%   OpenStreetMap street name, road class, compass direction of travel, length,
%   distance from the sensor, and whether it is the sensor's current match; for
%   the sensor, whether it is used (SensorUsed: false if excluded) and how its
%   edges were assigned (AssignmentSource: "manual override" after a correction). The
%   correct edge for a detector is normally the one whose street name matches the
%   street in the detector's name and whose direction matches its direction
%   words (">> NORD", "INGRESSO" towards the centre, ...).
%
%   Uses the results of runParmaTrafficTomography (struct or MAT file); no rerun.
%   Name/value: Radius (150 m), MinLength (0 m: shorter edges are listed but
%   marked), OutputFile ('parma_edge_candidates.csv'; written outside the library).
p = inputParser;
addParameter(p,'Radius',150); addParameter(p,'MinLength',0);
addParameter(p,'OutputFile','parma_edge_candidates.csv');
parse(p,varargin{:}); o = p.Results;
o.OutputFile = safeOutputPath(o.OutputFile);
if ischar(results) || (isstring(results) && isscalar(results))
    d = load(results,'results'); results = d.results;
end
if ischar(names), names = string(names); end            % a char path or name is ONE text, not characters
if isstring(names) && isscalar(names) && isfile(names)
    t = readtable(names,'TextType','string','Delimiter',',');
    assert(ismember('SectionName',t.Properties.VariableNames), ...
        'listEdgeCandidates:Columns','The file %s has no SectionName column.',names);
    names = t.SectionName;
elseif isstring(names) && isscalar(names) && (endsWith(names,".csv",'IgnoreCase',true) || contains(names,filesep))
    error('listEdgeCandidates:File','File not found: %s',names);
end
names = strtrim(string(names(:)));
G = results.G; M = results.sectionMatch;
en = G.Edges.EndNodes; if ~isnumeric(en), en = [findnode(G,en(:,1)) findnode(G,en(:,2))]; end
en = double(en);
lat = double(G.Nodes.Lat); lon = double(G.Nodes.Lon);
road = strings(size(en,1),1); cls = road;
if ismember('RoadName',G.Edges.Properties.VariableNames), road = string(G.Edges.RoadName); end
if ismember('Highway',G.Edges.Properties.VariableNames), cls = string(G.Edges.Highway); end
road(ismissing(road)) = ""; cls(ismissing(cls)) = "";
mName = regexprep(strtrim(string(M.SectionName)),'[\x00-\x1F\x7F\xA0]',' ');
rows = {};
for k = 1:numel(names)
    r = find(strcmpi(strtrim(mName),strtrim(regexprep(names(k),'[\x00-\x1F\x7F\xA0]',' '))),1);
    if isempty(r)
        warning('listEdgeCandidates:Name','Sensor "%s" not found.',names(k)); continue
    end
    la = double(M.Lat(r)); lo = double(M.Lon(r)); kx = cosd(la);
    mx = (lon(en(:,1))+lon(en(:,2)))/2; my = (lat(en(:,1))+lat(en(:,2)))/2;
    dist = 111320*hypot((mx-lo)*kx, my-la);
    len = 111320*hypot((lon(en(:,2))-lon(en(:,1)))*kx, lat(en(:,2))-lat(en(:,1)));
    br = mod(atan2d((lon(en(:,2))-lon(en(:,1)))*kx, lat(en(:,2))-lat(en(:,1))),360);   % 0 = north, 90 = east
    cur = double(M.DirectedEdgeIDs{r}(:));
    used = true; try, used = logical(M.Accepted(r)); catch, end
    src = ""; try, src = string(M.AssignmentSource(r)); catch, end
    near = find(dist <= o.Radius);
    [~,ord] = sort(dist(near)); near = near(ord);
    for e = near(:).'
        rows(end+1,:) = {mName(r),used,src,e,road(e),cls(e),compass(br(e)),round(br(e)),round(len(e),1), ...
            round(dist(e),1),ismember(e,cur),len(e) < o.MinLength}; %#ok<AGROW>
    end
end
vars = {'SectionName','SensorUsed','AssignmentSource','EdgeID','StreetName','RoadClass','Direction', ...
    'BearingDeg','LengthM','DistanceM','CurrentMatch','Short'};
if isempty(rows)
    T = cell2table(cell(0,numel(vars)),'VariableNames',vars);
    warning('listEdgeCandidates:Empty','No sensor of the list was found in the results; the table is empty.');
else
    T = cell2table(rows,'VariableNames',vars);
end
writetable(T,o.OutputFile);
fprintf('Edge candidates for %d sensors written to %s (%d rows).\n',numel(unique(T.SectionName)),o.OutputFile,height(T));
end

function c = compass(b)
d = ["N","NE","E","SE","S","SW","W","NW"];
c = d(mod(round(b/45),8)+1);
end
