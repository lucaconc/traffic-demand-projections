function file = exportSyntheticGraph(cfg)
%EXPORTSYNTHETICGRAPH  Bundle the road graph of the synthetic study with the library.
%
%   file = exportSyntheticGraph(cfg)
%
%   Saves the road graph of central Parma used by the synthetic study (from
%   the synthetic results, else from the Parma results) to
%   <library>/data/central_parma_osm_graph.mat. The synthetic stage then uses it
%   instead of downloading OpenStreetMap roads, so that others can reproduce the
%   synthetic results without network access. The graph is OpenStreetMap data,
%   distributed under the Open Database License (ODbL): keep the attribution.
if nargin < 1 || isempty(cfg), cfg = paperConfig(); end
src = {cfg.paths.syntheticMat, cfg.paths.parmaMat}; G = [];
for k = 1:numel(src)
    f = src{k}; if ~isfile(f), f = fullfile(cfg.outputRoot,f); end
    if isfile(f)
        d = load(f,'results');
        if isfield(d,'results') && isfield(d.results,'G') && ~isempty(d.results.G)
            G = d.results.G; fprintf('Road graph taken from %s\n',f); break
        end
    end
end
assert(~isempty(G),'exportSyntheticGraph:NoGraph','No saved road graph found: run the synthetic stage first.');
attribution = 'Road data (c) OpenStreetMap contributors, available under the Open Database License (ODbL), https://www.openstreetmap.org/copyright';
file = fullfile(fileparts(fileparts(mfilename('fullpath'))),'data','central_parma_osm_graph.mat');
if ~isfolder(fileparts(file)), mkdir(fileparts(file)); end
save(file,'G','attribution');
fprintf('Saved %s (%d nodes, %d edges).\n',file,numnodes(G),numedges(G));
end
