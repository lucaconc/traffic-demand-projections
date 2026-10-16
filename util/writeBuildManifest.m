function manifestFile = writeBuildManifest(cfg,status,report)
%WRITEBUILDMANIFEST Record how the paper artifacts were produced.
%
%   manifestFile = writeBuildManifest(cfg,status,report)
%
%   Writes build_manifest.json in cfg.outputRoot with the date, the MATLAB
%   release and installed toolboxes, the configuration, the stages run in
%   this call, and size and SHA-256 hash of every file read by the paper.
m = struct();
m.created = datestr(now,'yyyy-mm-ddTHH:MM:SS'); %#ok<TNOW1,DATST>
v = ver('MATLAB');
m.matlab = sprintf('%s %s',v.Version,v.Release);
tb = ver;
m.toolboxes = arrayfun(@(t)sprintf('%s %s',t.Name,t.Version),tb,'UniformOutput',false);
m.library = fileparts(fileparts(mfilename('fullpath')));
c = cfg;
c.synthetic.options = cellfun(@toText,c.synthetic.options,'UniformOutput',false);
c.parma.options = cellfun(@toText,c.parma.options,'UniformOutput',false);
c.figures.parmaOptions = cellfun(@toText,c.figures.parmaOptions,'UniformOutput',false);
m.config = c;
if nargin >= 3 && ~isempty(report)
    m.stages = rmfield(report,'message');
end
files = cell(height(status),1);
for k = 1:height(status)
    f = struct('file',char(status.File(k)),'role',char(status.Role(k)), ...
        'present',status.Present(k),'bytes',status.Bytes(k),'sha256','');
    if status.Present(k)
        fp = f.file;
        if ~(startsWith(fp,'/') || ~isempty(regexp(fp,'^[A-Za-z]:[\\/]','once')))
            fp = fullfile(cfg.outputRoot,fp);
        end
        f.sha256 = sha256File(fp);
    end
    files{k} = f;
end
m.files = [files{:}];
manifestFile = safeOutputPath(fullfile(cfg.outputRoot,'build_manifest.json'));
fid = fopen(manifestFile,'w');
assert(fid >= 0,'writeBuildManifest:Open','Cannot write %s.',manifestFile);
closer = onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(m,'PrettyPrint',true));
fprintf('Manifest written to %s\n',manifestFile);
end

function h = sha256File(file)
% SHA-256 through the Java runtime shipped with MATLAB; '' if unavailable.
h = '';
try
    fid = fopen(file,'r'); bytes = fread(fid,inf,'*uint8'); fclose(fid);
    md = java.security.MessageDigest.getInstance('SHA-256');
    d = typecast(md.digest(bytes),'uint8');
    h = lower(reshape(dec2hex(d,2).',1,[]));
catch
end
end

function s = toText(x)
if ischar(x) || isstring(x), s = char(x);
elseif isnumeric(x) || islogical(x), s = mat2str(x);
else, s = class(x);
end
end
