function zipFile = makeShareableLibrary(zipFile,varargin)
%MAKESHAREABLELIBRARY  Zip the library for sharing, without restricted data.
%
%   zipFile = makeShareableLibrary('tomography_paper_lib.zip')
%   zipFile = makeShareableLibrary(zipFile,'IncludeEdgeOverrides',true)
%
%   Copies the library to a temporary folder and zips it, leaving out
%     - data/parma_edge_overrides.csv, the corrections to the matching of the
%       Municipality's counting sections (metadata of the data provider; include
%       it only if the data agreement allows: 'IncludeEdgeOverrides',true);
%     - compiled MEX binaries (platform specific; rebuilt by buildPaperMex);
%     - every .mat file except data/central_parma_osm_graph.mat, so that no
%       result computed from the Parma counts is shared.
%   By default ('PaperText',false) the package produces the experiments' results
%   and figures only: the code that writes the paper's sentences (the paperdata
%   folder, its tests and runParmaScreeningSensitivity) is left out, and
%   cfg.paperText is set to false in the package's paperConfig.
%   The bundled road graph (exportSyntheticGraph) is kept: it is OpenStreetMap
%   data under the ODbL, with its attribution. The code is scanned for
%   absolute personal paths, which are reported. The Parma counts are never part
%   of the library.
ip = inputParser;
addParameter(ip,'IncludeEdgeOverrides',false,@(x)islogical(x)&&isscalar(x));
addParameter(ip,'PaperText',false,@(x)islogical(x)&&isscalar(x));
parse(ip,varargin{:}); o = ip.Results;
if nargin < 1 || isempty(zipFile), zipFile = fullfile(pwd,'tomography_paper_lib.zip'); end
libRoot = fileparts(fileparts(mfilename('fullpath')));
tmp = tempname; dst = fullfile(tmp,'tomography_paper_lib'); mkdir(dst);
cleanup = onCleanup(@() rmdir(tmp,'s'));
copyfile(libRoot,dst);
removed = {};
% restricted or unneeded files
if ~o.IncludeEdgeOverrides
    f = fullfile(dst,'data','parma_edge_overrides.csv');
    if isfile(f), delete(f); removed{end+1} = 'data/parma_edge_overrides.csv'; end
end
if ~o.PaperText
    % results and figures only: no code that writes the paper's text
    if isfolder(fullfile(dst,'paperdata')), rmdir(fullfile(dst,'paperdata'),'s'); removed{end+1} = 'paperdata/'; end
    for f = {fullfile('tests','testPaperData.m'), fullfile('parma','runParmaScreeningSensitivity.m')}
        if isfile(fullfile(dst,f{1})), delete(fullfile(dst,f{1})); removed{end+1} = strrep(f{1},'\','/'); end %#ok<AGROW>
    end
    cfgFile = fullfile(dst,'paperConfig.m'); txt = fileread(cfgFile);
    assert(contains(txt,'cfg.paperText = true;'),'makeShareableLibrary:Config','cfg.paperText not found in paperConfig.m.');
    txt = strrep(txt,'cfg.paperText = true;','cfg.paperText = false;');
    fid = fopen(cfgFile,'w'); fwrite(fid,txt); fclose(fid);
    readme = fullfile(dst,'README.md');
    if isfile(readme)
        fid = fopen(readme,'a');
        fprintf(fid,['\n## About this version\n\nThis version reproduces the experiments of the paper and ' ...
            'their figures. The code that wrote the paper''s numerical sentences and captions is not ' ...
            'included: cfg.paperText is false, and the stages ''captions'' and ''macros'' are reported ' ...
            'as skipped. Numerical results are in the results files (.mat) of each study.\n']);
        fclose(fid);
    end
end
files = dir(fullfile(dst,'**','*'));
files = files(~[files.isdir]);
for k = 1:numel(files)
    f = fullfile(files(k).folder,files(k).name); rel = strrep(f,[dst filesep],'');
    [~,~,ext] = fileparts(f);
    isMex = startsWith(ext,'.mex');
    isMat = strcmpi(ext,'.mat') && ~strcmp(strrep(rel,'\','/'),'data/central_parma_osm_graph.mat');
    isJunk = any(strcmpi(ext,{'.asv','.DS_Store'})) || strcmp(files(k).name,'.DS_Store');
    if isMex || isMat || isJunk, delete(f); removed{end+1} = rel; end %#ok<AGROW>
end
% personal paths in the code
code = dir(fullfile(dst,'**','*.m')); flagged = {};
code = code(~strcmp({code.name},'makeShareableLibrary.m'));   % its own search pattern is not a path
for k = 1:numel(code)
    txt = fileread(fullfile(code(k).folder,code(k).name));
    if ~isempty(regexp(txt,'(/Users/|/home/[a-z]|[A-Z]:\\Users\\)','once'))
        flagged{end+1} = strrep(fullfile(code(k).folder,code(k).name),[dst filesep],''); %#ok<AGROW>
    end
end
if isfile(zipFile), delete(zipFile); end
zip(zipFile,'tomography_paper_lib',tmp);
fprintf('Shareable library: %s\n',zipFile);
if isempty(removed), fprintf('Nothing removed.\n'); else, fprintf('Left out: %s\n',strjoin(removed,', ')); end
if ~isfile(fullfile(dst,'data','central_parma_osm_graph.mat'))
    fprintf(['No bundled road graph: others will need internet access for the synthetic study.\n' ...
        'Run exportSyntheticGraph(cfg) first to include it.\n']);
end
if ~isempty(flagged)
    warning('makeShareableLibrary:PersonalPaths','Absolute personal paths found in: %s',strjoin(flagged,', '));
end
end
