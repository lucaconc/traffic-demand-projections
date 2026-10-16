function conflicts = checkLibraryShadowing()
%CHECKLIBRARYSHADOWING Find other copies of library functions that MATLAB
%would call instead of the library versions.
%
%   conflicts = checkLibraryShadowing()
%
%   MATLAB resolves a function in the current folder before the path, so an
%   old copy of, e.g., generateParmaRequestedFigures.m left in the paper's
%   matlab/ folder silently replaces the library version. Returns a cell array
%   of 'name: file' strings, empty when every library function resolves to
%   the library.
root = fileparts(fileparts(mfilename('fullpath')));
sub = {'estimator','synthetic','benchmark','parma','figures','paperdata','util'};
conflicts = {};
for k = 1:numel(sub)
    files = dir(fullfile(root,sub{k},'*.m'));
    for j = 1:numel(files)
        [~,name] = fileparts(files(j).name);
        w = which(name);
        mine = fullfile(root,sub{k},files(j).name);
        if ~isempty(w) && ~strcmp(canonical(w),canonical(mine))
            conflicts{end+1} = sprintf('%s: %s',name,w); %#ok<AGROW>
        end
    end
end
end

function p = canonical(p)
p = lower(strrep(char(p),'\','/'));
end
