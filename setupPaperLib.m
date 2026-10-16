function info = setupPaperLib(varargin)
%SETUPPAPERLIB Add the library to the MATLAB path and check its requirements.
%
%   info = setupPaperLib()          add folders, warn about missing requirements
%   info = setupPaperLib('Strict',true)   error instead of warning
%
%   The library needs MATLAB R2020b or newer (tiledlayout, exportgraphics),
%   the Optimization Toolbox (lsqnonneg, fmincon) and the Statistics and
%   Machine Learning Toolbox (poissrnd, gamrnd, prctile in the synthetic
%   study).  No other toolbox is used.  Internet access is needed only to
%   download OpenStreetMap roads the first time (results are cached).
%
%   info.root       library folder
%   info.folders    folders added to the path
%   info.missing    names of missing requirements (empty when complete)
%   info.notes      warnings about files and toolboxes needed by the Parma stage

p = inputParser;
addParameter(p,'Strict',false,@(x)islogical(x)&&isscalar(x));
parse(p,varargin{:});

root = fileparts(mfilename('fullpath'));
sub = {'estimator','synthetic','benchmark','parma','figures','paperdata','util'};
folders = [{root}, fullfile(root,sub)];
folders = folders(cellfun(@isfolder,folders));   % the public package has no paperdata folder
addpath(folders{:});
% Compiled MEX files live OUTSIDE the library (buildPaperMex), in a sibling folder;
% they are optional: without them the MATLAB code is used.
mexDir = fullfile(fileparts(root),'tomography_paper_mex');
if isfolder(mexDir), addpath(mexDir); end

missing = {};
if exist('OCTAVE_VERSION','builtin') ~= 0
    missing{end+1} = 'MATLAB R2020b or newer (GNU Octave is not supported)';
elseif verLessThan('matlab','9.9')      % R2020b
    missing{end+1} = 'MATLAB R2020b or newer';
end
required = {'lsqnonneg','Optimization Toolbox'; 'fmincon','Optimization Toolbox'; ...
            'poissrnd','Statistics and Machine Learning Toolbox'; ...
            'gamrnd','Statistics and Machine Learning Toolbox'};
for k = 1:size(required,1)
    if exist(required{k,1},'file') == 0
        missing{end+1} = sprintf('%s (%s)',required{k,2},required{k,1}); %#ok<AGROW>
    end
end
missing = unique(missing,'stable');

% Parma stage: files from the original project and optional toolboxes.
notes = {};
if exist('findStableParmaInterval','file') == 0
    notes{end+1} = ['findStableParmaInterval.m is missing: the parma stage needs it ' ...
        '(AutoStableDateRange=true selects the calibration period). Copy it to parma/.'];
end
if exist('diagnoseParmaFanoDecomposition','file') == 0
    notes{end+1} = ['diagnoseParmaFanoDecomposition.m is missing: optional diagnostic only ' ...
        '(set ''FanoDecomposition'',false in cfg.parma.options to silence the warning).'];
end
if exist('shaperead','file') == 0 || exist('projcrs','file') == 0 || exist('projinv','file') == 0
    notes{end+1} = ['Mapping Toolbox (shaperead/projcrs/projinv) is missing: sensor matching ' ...
        'falls back to OpenStreetMap only, without the official street layer of the Municipality.'];
end
for k = 1:numel(notes)
    warning('setupPaperLib:Parma','%s',notes{k});
end

if ~isempty(missing)
    msg = sprintf('Missing requirements: %s.',strjoin(missing,'; '));
    if p.Results.Strict
        error('setupPaperLib:Requirements','%s',msg);
    else
        warning('setupPaperLib:Requirements','%s',msg);
    end
end
% Interactive use: other copies of library functions (e.g. old files in the
% current folder) would be called instead of the library versions.
% buildPaperArtifacts avoids this by itself; here we only warn.
shadow = checkLibraryShadowing();
if ~isempty(shadow)
    warning('setupPaperLib:Shadowed', ['These files shadow the library when called from ' ...
        'the current folder or path (buildPaperArtifacts is not affected):\n  %s'], ...
        strjoin(shadow,'\n  '));
end
info = struct('root',root,'folders',{folders},'missing',{missing},'notes',{notes},'shadowed',{shadow});
if nargout == 0
    fprintf('Paper library ready: %s\n',root);
end
end
