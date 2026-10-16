function info = buildPaperMex(varargin)
%BUILDPAPERMEX Compile the optional MEX files of the library, outside the library.
%
%   info = buildPaperMex()
%   info = buildPaperMex('OutputDir','/path/to/tomography_paper_mex')
%
%   Compiles estimator/mex/compoundLogPmfMex.c (the compound Poisson log-probability
%   recursion, about half of the GMM's run time) into OutputDir, by default the
%   sibling folder <parent of the library>/tomography_paper_mex, which setupPaperLib
%   adds to the path when it exists. Nothing is written inside the library.
%
%   After compiling, the MEX result is compared with the MATLAB code on random
%   cases; the build fails if they differ by more than 1e-12 (relative).
%   Needs a C compiler configured for MATLAB (mex -setup C; on macOS, Xcode's
%   command line tools). Run "clear functions" afterwards so that already loaded
%   functions see the new MEX file.
libRoot = fileparts(fileparts(mfilename('fullpath')));
p = inputParser;
addParameter(p,'OutputDir',fullfile(fileparts(libRoot),'tomography_paper_mex'));
parse(p,varargin{:});
outDir = char(safeOutputPath(p.Results.OutputDir));
if ~isfolder(outDir), mkdir(outDir); end
src = fullfile(libRoot,'estimator','mex','compoundLogPmfMex.c');
fprintf('Compiling %s\n      -> %s\n',src,outDir);
mex('-O','-outdir',outDir,src);
addpath(outDir); clear('compoundLogPmfMex');
% Check against the MATLAB recursion (vector path, independent of the MEX switch)
rs = RandStream('mt19937ar','Seed',7); worst = 0; bad = 0;
for r = 1:300
    K = randi(rs,5); jumps = sort(randperm(rs,6,K)); lambda = 30*rand(rs,1,K);
    if rand(rs) < 0.2, lambda(randi(rs,K)) = 0; end
    zmax = randi(rs,600)-1;
    a = compoundLogPmfMex(lambda(:),jumps(:),zmax);
    b = referencePmf(lambda,jumps,zmax);
    if ~isequal(isfinite(a),isfinite(b)), bad = bad+1; continue; end
    f = isfinite(a);
    if any(f), worst = max(worst, max(abs(a(f)-b(f))./max(1,abs(b(f))))); end
end
assert(bad == 0 && worst <= 1e-12, 'buildPaperMex:Mismatch', ...
    'MEX and MATLAB results differ (%d pattern mismatches, max relative difference %.2e).',bad,worst);
fprintf('MEX check passed on 300 random cases (max relative difference %.1e).\n',worst);
fprintf('Run "clear functions" so that loaded functions use the MEX file.\n');
info = struct('outputDir',outDir,'maxRelDiff',worst);
end

function lp = referencePmf(lambda,jumps,zmax)
% The MATLAB recursion of estimateMuEfficientGMM (small-vector version).
lp = -inf(zmax+1,1); lp(1) = -sum(lambda);
K = max(jumps); coef = zeros(1,K); coef(jumps) = jumps(:).'.*lambda(:).';
for z = 1:zmax
    k = 1:min(z,numel(coef)); k = k(coef(k)>0);
    if isempty(k), continue; end
    t = log(coef(k)) + lp(z-k+1).';
    mx = max(t); if ~isfinite(mx), continue; end
    lp(z+1) = mx + log(sum(exp(t-mx))) - log(z);
end
end
