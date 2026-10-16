function inst = cityRoutingInstance(nPaths,nSensors,seed,varargin)
%CITYROUTINGINSTANCE Random city-like routing instance.
%
%   inst = cityRoutingInstance(nPaths,nSensors,seed)
%   inst = cityRoutingInstance(40,8,1,'GridSize',10,'CostNoise',0.3)
%
%   Network: a GridSize x GridSize grid of intersections with two-way links
%   (each direction a separate directed link) and random base costs in
%   [1, 1.5] (heterogeneous street speeds). Paths: nPaths distinct shortest
%   paths between random origin-destination pairs, each under link costs
%   multiplied by exp(CostNoise*xi), xi standard normal (stochastic route
%   choice). Sensors: nSensors directed links among those used by the paths,
%   drawn without replacement with probability proportional to the number of
%   paths using them (busy roads are more often instrumented).
%
%   Paths crossing the same sensors form a class; paths crossing no sensor are
%   unobservable and dropped. Fields:
%     A           sensors x classes routing matrix (distinct nonzero columns)
%     pathToClass class of each path (0: unobservable)
%     classSize   number of paths per class
%     Apath       sensors x paths matrix before merging
%     sensorLinks, linkUse, nNodes, nLinks
%   Path rates mu_path give class rates A-merged as accumarray(pathToClass, mu_path).
%   No toolbox is needed (own Dijkstra), so instances are reproducible anywhere.

p = inputParser;
addParameter(p,'GridSize',10); addParameter(p,'CostNoise',0.3);
parse(p,varargin{:}); o = p.Results;
st = seed;                                   % portable linear congruential stream
    function u = U()
        st = mod(1103515245*st + 12345, 2^31); u = st/2^31;
    end
    function z = Z()
        z = sqrt(-2*log(max(U(),1e-12)))*cos(2*pi*U());
    end
g = o.GridSize; nNodes = g*g; id = @(r,c) (r-1)*g + c;
E = zeros(0,2);
for r = 1:g
    for c = 1:g
        if c < g, E(end+1,:) = [id(r,c) id(r,c+1)]; E(end+1,:) = [id(r,c+1) id(r,c)]; end %#ok<AGROW>
        if r < g, E(end+1,:) = [id(r,c) id(r+1,c)]; E(end+1,:) = [id(r+1,c) id(r,c)]; end %#ok<AGROW>
    end
end
nLinks = size(E,1);
base = zeros(nLinks,1); for e = 1:nLinks, base(e) = 1 + 0.5*U(); end
out = cell(nNodes,1); for e = 1:nLinks, out{E(e,1)}(end+1) = e; end
paths = {}; keys = {}; tries = 0;
while numel(paths) < nPaths && tries < 50*nPaths
    tries = tries + 1;
    a = 1 + floor(U()*nNodes); b = 1 + floor(U()*nNodes);
    if a == b, continue; end
    w = base; for e = 1:nLinks, w(e) = w(e)*exp(o.CostNoise*Z()); end
    pe = dijkstraPath(a,b,w,E,out,nNodes);
    if isempty(pe), continue; end
    k = sprintf('%d,',pe);
    if any(strcmp(keys,k)), continue; end
    paths{end+1} = pe; keys{end+1} = k; %#ok<AGROW>
end
nP = numel(paths);
linkUse = zeros(nLinks,1);
for j = 1:nP, linkUse(paths{j}) = linkUse(paths{j}) + 1; end
cand = find(linkUse > 0); sensors = zeros(0,1);
wts = linkUse(cand);
for s = 1:min(nSensors,numel(cand))
    cw = cumsum(wts)/sum(wts); k = find(cw >= U(),1);
    sensors(end+1,1) = cand(k); wts(k) = 0; %#ok<AGROW>
    if sum(wts) == 0, break; end
end
Apath = zeros(numel(sensors),nP);
for j = 1:nP, Apath(:,j) = double(ismember(sensors,paths{j})); end
obs = any(Apath,1);
% Classes in order of first appearance (portable stable unique of rows)
[u,ia,ic] = unique(Apath(:,obs).','rows');
[~,ord] = sort(ia); rk = zeros(numel(ord),1); rk(ord) = 1:numel(ord);
cols = u(ord,:); cls = rk(ic);
pathToClass = zeros(nP,1); pathToClass(obs) = cls;
inst = struct('A',double(cols.'),'pathToClass',pathToClass, ...
    'classSize',accumarray(cls(:),1).','Apath',Apath,'sensorLinks',sensors, ...
    'linkUse',linkUse,'nNodes',nNodes,'nLinks',nLinks,'nPaths',nP);
end

function pe = dijkstraPath(a,b,w,E,out,nNodes)
% Shortest path from a to b (edge list), positive weights.
D = inf(nNodes,1); prevE = zeros(nNodes,1); done = false(nNodes,1); D(a) = 0;
while true
    dd = D; dd(done) = inf; [dv,u] = min(dd);
    if ~isfinite(dv) || u == b, break; end
    done(u) = true;
    for e = out{u}
        v = E(e,2); nd = dv + w(e);
        if nd < D(v), D(v) = nd; prevE(v) = e; end
    end
end
if ~isfinite(D(b)), pe = []; return; end
pe = zeros(0,1); v = b;
while v ~= a, e = prevE(v); pe(end+1,1) = e; v = E(e,1); end %#ok<AGROW>
pe = flipud(pe);
end
