function T = runActivePathsSweep(varargin)
%RUNACTIVEPATHSSWEEP  How the estimators compare as more paths carry traffic.
%
%   T = runActivePathsSweep()
%   T = runActivePathsSweep('ActiveCounts',[6 12 24 36 48 60],'Replicates',20)
%
%   Runs the synthetic study of Section 5.3 at its baseline only (N = 500, clean
%   counts, all sensors, synchronous counts) for an increasing number of active
%   paths among the candidates, and returns and prints the median relative error
%   of each estimator. The other experiments are reduced to their baseline to
%   save time. Outputs go to <outputRoot>/active_paths_sweep, so the paper's
%   results are not touched.
%
%   Options: 'ActiveCounts' (default [6 12 24 36 48 60]), 'Replicates' (20),
%   'OutputDir' (default <paperConfig().outputRoot>/active_paths_sweep), 'Graph'
%   (default: the road graph saved by the synthetic run, else by the Parma run, so
%   that nothing is downloaded and the paper's graph is used), and any
%   extra name/value pairs of runSyntheticOSMTomographyTests in 'StudyOptions'.
ip = inputParser;
addParameter(ip,'ActiveCounts',[6 12 24 36 48 60],@(x)isnumeric(x)&&all(x>0));
addParameter(ip,'Replicates',20,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(ip,'OutputDir','');
addParameter(ip,'StudyOptions',{},@iscell);
addParameter(ip,'Graph',[]);
parse(ip,varargin{:}); o = ip.Results;
cfg = paperConfig();
if isempty(o.OutputDir), o.OutputDir = fullfile(cfg.outputRoot,'active_paths_sweep'); end
G = o.Graph;
if isempty(G)
    % Reuse the road graph of an earlier run instead of downloading it again.
    if isfield(cfg.synthetic,'graphFile') && isfile(cfg.synthetic.graphFile)
        g = load(cfg.synthetic.graphFile,'G'); G = g.G; fprintf('Road graph bundled with the library: %s\n',cfg.synthetic.graphFile);
    end
    candidates = {cfg.paths.syntheticMat, cfg.paths.parmaMat};
    for k = 1:numel(candidates)*isempty(G)
        f = candidates{k};
        if ~isfile(f), f = fullfile(cfg.outputRoot,f); end
        if isfile(f)
            data = load(f,'results');
            if isfield(data,'results') && isfield(data.results,'G') && ~isempty(data.results.G)
                G = data.results.G; fprintf('Road graph loaded from %s\n',f); break
            end
        end
    end
    assert(~isempty(G),'runActivePathsSweep:NoGraph', ...
        ['No saved road graph found (looked in %s and %s). Run the synthetic or parma stage ' ...
         'once, or pass the graph: runActivePathsSweep(''Graph'',G).'],candidates{:});
end
tags  = {'RelativeRateError','Err_GMM','Err_Pair','Err_CostCL','Err_GLS','Err_Cum','Err_NNLS'};
names = {'PCL','GMM','PAIR','CA_PCL','GLS','CUM','NNLS'};
E = nan(numel(o.ActiveCounts),numel(tags));
for k = 1:numel(o.ActiveCounts)
    R = runSyntheticOSMTomographyTests('NumTruePaths',o.ActiveCounts(k),'Replicates',o.Replicates, ...
        'ObservationCounts',500,'OutlierCounts',0,'ImpurityFractions',0,'SensorCounts',30, ...
        'TravelSpeeds',Inf,'FullLikelihoodReplicates',0,'Visible','off','Graph',G, ...
        'OutputDir',fullfile(o.OutputDir,sprintf('active_%02d',o.ActiveCounts(k))),o.StudyOptions{:});
    raw = R.raw;
    base = strcmp(string(raw.Experiment),"observations") & raw.N == 500;
    for t = 1:numel(tags)
        if ismember(tags{t},raw.Properties.VariableNames)
            E(k,t) = median(raw.(tags{t})(base),'omitnan');
        end
    end
    fprintf('done: %d active paths\n',o.ActiveCounts(k));
end
T = array2table([o.ActiveCounts(:) E],'VariableNames',[{'Active'} names]);
fprintf('\nMedian relative error at N = 500 (%d replicates)\n',o.Replicates); disp(T);
if ~isfolder(o.OutputDir), mkdir(o.OutputDir); end
save(fullfile(o.OutputDir,'active_paths_sweep.mat'),'T');
end
