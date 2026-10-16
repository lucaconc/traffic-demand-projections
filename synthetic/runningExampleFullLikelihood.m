function out=runningExampleFullLikelihood(varargin)
%RUNNINGEXAMPLEFULLLIKELIHOOD Composite versus exact likelihood on the
%running example of the paper (two sensors, three paths).
%
%   out = runningExampleFullLikelihood('N',500,'R',200,'Mu',[4;3;5],'Seed',1)
%
%   For the running example, Y1=X1+X3 and Y2=X2+X3, so the exact probability
%   of (y1,y2) is a finite sum over x3=0..min(y1,y2).  Each Monte Carlo
%   replicate fits both the exact maximum-likelihood estimator (fminsearch on
%   log mu, started from the composite estimate) and the composite estimator
%   of estimateMuOptimizedProjections.  out.ratio = MSE(composite)/MSE(exact)
%   measures the efficiency lost by the composite formulation; out.macro is
%   the text used for \ResSynFullLik.
ip=inputParser;
addParameter(ip,'N',500);
addParameter(ip,'R',200);
addParameter(ip,'Mu',[4;3;5]);
addParameter(ip,'Seed',1);
parse(ip,varargin{:});
o=ip.Results;
A=[1 0 1;0 1 1];
mu=o.Mu(:);
rng(o.Seed);
errC=nan(o.R,3); errE=nan(o.R,3);
for r=1:o.R
    X=[poissonSample(mu(1),o.N) poissonSample(mu(2),o.N) poissonSample(mu(3),o.N)];
    Y=X*A.';
    muC=estimateMuOptimizedProjections(A,Y,'Display','off','ComputeSandwich',false);
    nll=@(th) -exactLogLik(exp(th),Y);
    th=fminsearch(nll,log(max(muC(:),1e-3)),optimset('TolX',1e-8,'TolFun',1e-10, ...
        'MaxFunEvals',4000,'MaxIter',4000,'Display','off'));
    errC(r,:)=(muC(:)-mu).';
    errE(r,:)=(exp(th(:))-mu).';
end
out.mseComposite=mean(sum(errC.^2,2));
out.mseExact=mean(sum(errE.^2,2));
out.ratio=out.mseComposite/out.mseExact;
out.N=o.N; out.R=o.R; out.mu=mu;
out.macro=sprintf(['the ratio %.2f between the mean squared errors of the composite and of the exact ' ...
    'maximum-likelihood estimators ($N=%d$, %d replicates, $\\mu=(%s)^T$)'], ...
    out.ratio,o.N,o.R,strjoin(arrayfun(@(x)num2str(x),mu.','UniformOutput',false),','));
end

function ll=exactLogLik(mu,Y)
% log p(y1,y2) = log sum_{x3} Pois(y1-x3;mu1) Pois(y2-x3;mu2) Pois(x3;mu3)
[U,~,iu]=unique(Y,'rows');
cnt=accumarray(iu,1);
lp=zeros(size(U,1),1);
for k=1:size(U,1)
    x3=(0:min(U(k,:))).';
    t=lpois(U(k,1)-x3,mu(1))+lpois(U(k,2)-x3,mu(2))+lpois(x3,mu(3));
    mx=max(t);
    lp(k)=mx+log(sum(exp(t-mx)));
end
ll=sum(cnt.*lp);
end

function v=lpois(x,lam)
v=x.*log(lam)-lam-gammaln(x+1);
end

function x=poissonSample(lam,n)
% Inversion sampler (no Statistics Toolbox needed); adequate for small means.
x=zeros(n,1);
u=rand(n,1);
p0=exp(-lam);
for i=1:n
    k=0; p=p0; F=p0;
    while u(i)>F
        k=k+1; p=p*lam/k; F=F+p;
    end
    x(i)=k;
end
end
