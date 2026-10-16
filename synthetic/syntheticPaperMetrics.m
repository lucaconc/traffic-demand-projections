function M=syntheticPaperMetrics(muHat,muTrue,A,Y)
%SYNTHETICPAPERMETRICS Error measures of one synthetic fit, as in the paper.
%
%   M = syntheticPaperMetrics(muHat,muTrue)
%   M = syntheticPaperMetrics(muHat,muTrue,A,Y)
%
%   muHat, muTrue : estimated and true class intensities (same class order).
%   A, Y          : optional class matrix and the observations (N-by-m) used
%                   for the fit; if given, the NNLS baseline is computed from
%                   the sample means, as in the initialization of the paper.
%
%   M.E     relative error ||muHat-muTrue||/||muTrue||          (eq. synthetic-error)
%   M.Eact  same error restricted to the active classes (muTrue>0)
%   M.F0    fraction of the estimated flow assigned to inactive classes
%   M.ENNLS relative error of the NNLS fit of the means (NaN if A, Y absent)
%   M.muNNLS the NNLS fit itself
muHat=muHat(:); muTrue=muTrue(:);
act=muTrue>0;
M.E=norm(muHat-muTrue)/max(norm(muTrue),eps);
M.Eact=norm(muHat(act)-muTrue(act))/max(norm(muTrue(act)),eps);
M.F0=sum(muHat(~act))/max(sum(muHat),eps);
M.ENNLS=NaN;
M.muNNLS=[];
if nargin>=4 && ~isempty(A)
    muN=lsqnonneg(A,mean(Y,1).',optimset('Display','off'));
    M.muNNLS=muN;
    M.ENNLS=norm(muN-muTrue)/max(norm(muTrue),eps);
end
end
