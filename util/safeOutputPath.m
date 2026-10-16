function p = safeOutputPath(p)
%SAFEOUTPUTPATH Absolute output path guaranteed to lie outside the library.
%
%   p = safeOutputPath(p)
%
%   Every function of the library that writes a file or creates a folder passes
%   its output path through this function, so that nothing is ever written
%   inside the library folder tree:
%     * an empty path is returned unchanged (meaning "do not write");
%     * a relative path is made absolute with respect to the current folder,
%       unless the current folder is inside the library: then it is resolved
%       in the sibling folder <parent of the library>/tomography_paper_output
%       (a warning says so, once per session);
%     * an absolute path inside the library raises an error.
%   See also outputBaseFolder.
if isempty(p) || (isstring(p) && all(strlength(p) == 0)), return; end
wasString = isstring(p);
p = char(p);
if ~isAbsolutePath(p), p = fullfile(outputBaseFolder(), p); end
p = normalizePath(p);
root = normalizePath(libraryRootFolder());
if isInsideFolder(p,root)
    error('tomography:OutputInsideLibrary', ...
        ['Output path %s lies inside the library folder %s. Write outputs outside the ' ...
         'library (for example in the paper folder).'],p,root);
end
if wasString, p = string(p); end
end

function b = outputBaseFolder()
% Folder for relative output paths: the current folder, unless it lies inside
% the library tree.
persistent warned
root = normalizePath(libraryRootFolder()); here = normalizePath(pwd);
if isInsideFolder(here,root)
    b = fullfile(fileparts(root),'tomography_paper_output');
    if isempty(warned)
        warning('tomography:OutputBase', ['The current folder is inside the library; ' ...
            'relative outputs are written to %s.'],b);
        warned = true;
    end
else
    b = here;
end
end

function r = libraryRootFolder()
% The library root is the parent of this file's folder (util/).
r = fileparts(fileparts(mfilename('fullpath')));
end

function tf = isAbsolutePath(p)
tf = startsWith(p,'/') || startsWith(p,'\\') || ~isempty(regexp(p,'^[A-Za-z]:[\\/]','once'));
end

function p = normalizePath(p)
% Resolve '.' and '..' and unify separators, without requiring the path to exist.
p = strrep(p,'\',filesep); p = strrep(p,'/',filesep);
lead = ''; drive = regexp(p,'^[A-Za-z]:','match','once');
if ~isempty(drive), lead = drive; p = p(numel(drive)+1:end); end
isAbs = startsWith(p,filesep);
parts = strsplit(p,filesep); out = {};
for k = 1:numel(parts)
    s = parts{k};
    if isempty(s) || strcmp(s,'.'), continue; end
    if strcmp(s,'..')
        if ~isempty(out), out(end) = []; end
    else
        out{end+1} = s; %#ok<AGROW>
    end
end
p = strjoin(out,filesep);
if isAbs, p = [filesep p]; end
p = [lead p];
end

function tf = isInsideFolder(p,root)
% True if p equals root or lies below it (case-insensitive on Windows).
if ispc, p = lower(p); root = lower(root); end
tf = strcmp(p,root) || startsWith(p,[root filesep]);
end
