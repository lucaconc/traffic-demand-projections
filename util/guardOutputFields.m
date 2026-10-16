function o = guardOutputFields(o,fields)
%GUARDOUTPUTFIELDS Pass the named output options through safeOutputPath.
%   o = guardOutputFields(o,{'OutputDir','SaveFile',...}) makes each existing,
%   nonempty field an absolute path outside the library (error if inside).
for k = 1:numel(fields)
    f = fields{k};
    if isfield(o,f) && ~isempty(o.(f)), o.(f) = safeOutputPath(o.(f)); end
end
end
