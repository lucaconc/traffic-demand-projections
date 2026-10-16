function out = generateNRMSEvsSizeFigure(results,varargin)
%GENERATENRMSEVSSIZEFIGURE Per-stratum test NRMSE against stratum size.
%
%   out = generateNRMSEvsSizeFigure(results)
%   out = generateNRMSEvsSizeFigure('parma_tomography_results.mat','OutputDir','diagnostics')
%
%   Diagnostic, not used by the paper. One point per reliably estimated stratum:
%     (a) test NRMSE against the number of training intervals (data for the fit);
%     (b) test NRMSE against the number of test intervals (noise of the test means).
%   Colors give the day type; black squares are the median NRMSE in size bins
%   (quantiles of the size). Each panel reports Spearman's rank correlation rho
%   with a two-sided permutation p-value (no toolbox needed). A clearly negative
%   rho in (a) means that small strata are fitted worse; in (b), that small test
%   sets inflate the error through noise in the test means.
%
%   Name/value: OutputDir ('diagnostics'), Bins (6), Permutations (5000),
%   Visible ('off'), Resolution (300), Seed (1).
%   Writes nrmse_vs_stratum_size.pdf/.png/.fig and nrmse_vs_stratum_size.csv.

p = inputParser;
addParameter(p,'OutputDir','diagnostics',@(x)ischar(x)||isstring(x));
addParameter(p,'Bins',6,@(x)isnumeric(x)&&isscalar(x)&&x>=2);
addParameter(p,'Permutations',5000,@(x)isnumeric(x)&&isscalar(x)&&x>=100);
addParameter(p,'Visible','off',@(x)any(strcmpi(x,{'on','off'})));
addParameter(p,'Resolution',300,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'Seed',1,@(x)isnumeric(x)&&isscalar(x));
parse(p,varargin{:}); o = p.Results;
% All output paths lie outside the library (see safeOutputPath).
o = guardOutputFields(o,{'OutputDir'});
if ischar(results) || (isstring(results) && isscalar(results))
    d = load(results,'results'); results = d.results;
end
T = results.strataTable;
need = {'TestNRMSE','nTrainVec','nTestVec','DayType'};
assert(all(ismember(need,T.Properties.VariableNames)), ...
    'generateNRMSEvsSizeFigure:Fields','strataTable lacks %s.',strjoin(need,', '));
valid = results.validFitByStratum(:);
y = double(T.TestNRMSE); ntr = double(T.nTrainVec); nte = double(T.nTestVec);
use = valid & isfinite(y) & ntr>0 & nte>0;
dayType = string(T.DayType);
types = unique(dayType(use),'stable');
palette = [.3 .35 .45; .05 .38 .75; .9 .33 .08; .2 .65 .4];
rs = RandStream('mt19937ar','Seed',o.Seed);

f = figure('Color','w','Visible','off','Units','pixels','Position',[40 40 1300 520]);
tl = tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');
xs = {ntr,nte};
xlabels = {'Training intervals in the stratum','Test intervals in the stratum'};
panel = {'(a)','(b)'};
stats = zeros(2,2);
for k = 1:2
    ax = nexttile(tl); hold(ax,'on');
    x = xs{k};
    for t = 1:numel(types)
        sel = use & dayType==types(t);
        scatter(ax,x(sel),y(sel),22,palette(min(t,4),:),'filled','MarkerFaceAlpha',.6, ...
            'DisplayName',char(types(t)));
    end
    [bx,by] = binnedMedians(x(use),y(use),o.Bins);
    plot(ax,bx,by,'ks-','MarkerFaceColor','k','LineWidth',1.5,'DisplayName','median in size bins');
    [rho,pval] = spearmanPermutation(x(use),y(use),o.Permutations,rs);
    stats(k,:) = [rho pval];
    set(ax,'XScale','log'); grid(ax,'on'); box(ax,'on');
    xlabel(ax,xlabels{k}); ylabel(ax,'Test NRMSE');
    title(ax,sprintf('%s Spearman \\rho = %.2f, p = %s (%d strata)',panel{k},rho,pformat(pval,o.Permutations),nnz(use)), ...
        'FontWeight','normal');
    if k == 1, legend(ax,'Location','northeast'); end
end
folder = char(o.OutputDir); if ~isfolder(folder), mkdir(folder); end
base = fullfile(folder,'nrmse_vs_stratum_size');
drawnow;
exportgraphics(tl,[base '.pdf'],'ContentType','vector','BackgroundColor','white');
exportgraphics(tl,[base '.png'],'Resolution',o.Resolution,'BackgroundColor','white');
savefig(f,[base '.fig']);
idx = find(use);
writetable(table(idx,dayType(idx),ntr(idx),nte(idx),y(idx), ...
    'VariableNames',{'Stratum','DayType','TrainingIntervals','TestIntervals','TestNRMSE'}),[base '.csv']);
set(f,'Visible',o.Visible);
fprintf(['NRMSE vs stratum size (%d strata): Spearman rho %.2f (p %s) with training intervals, ' ...
    '%.2f (p %s) with test intervals.\n'],nnz(use),stats(1,1),pformat(stats(1,2),o.Permutations), ...
    stats(2,1),pformat(stats(2,2),o.Permutations));
out = struct('figure',f,'rhoTrain',stats(1,1),'pTrain',stats(1,2),'rhoTest',stats(2,1),'pTest',stats(2,2), ...
    'files',string({[base '.pdf'],[base '.png'],[base '.fig'],[base '.csv']}));
end

function [bx,by] = binnedMedians(x,y,nb)
% Medians of x and y in bins with equal numbers of points (by x quantiles).
[xs,ord] = sort(x(:)); ys = y(ord); n = numel(xs);
nb = max(1,min(nb,floor(n/5)));
edges = round(linspace(0,n,nb+1));
bx = nan(nb,1); by = nan(nb,1);
for b = 1:nb
    r = edges(b)+1:edges(b+1);
    if isempty(r), continue; end
    bx(b) = median(xs(r)); by(b) = median(ys(r));
end
end

function [rho,p] = spearmanPermutation(x,y,nperm,rs)
% Spearman rank correlation (average ranks for ties) and two-sided permutation p.
rx = tiedRanks(x(:)); ry = tiedRanks(y(:));
rho = corrPearson(rx,ry);
if ~isfinite(rho), p = NaN; return; end
count = 0;
for k = 1:nperm
    if abs(corrPearson(rx,ry(randperm(rs,numel(ry))))) >= abs(rho), count = count+1; end
end
p = (count+1)/(nperm+1);
end

function r = tiedRanks(v)
[s,ord] = sort(v); r = zeros(size(v)); n = numel(v); k = 1;
while k <= n
    j = k; while j < n && s(j+1) == s(k), j = j+1; end
    r(ord(k:j)) = (k+j)/2; k = j+1;
end
end

function c = corrPearson(a,b)
a = a-mean(a); b = b-mean(b);
c = (a.'*b)/max(sqrt((a.'*a)*(b.'*b)),eps);
end

function s = pformat(p,nperm)
if ~isfinite(p), s = 'n/a';
elseif p <= 1/(nperm+1)*1.0001, s = sprintf('< %.1g',2/(nperm+1));
else, s = sprintf('%.3f',p);
end
end
