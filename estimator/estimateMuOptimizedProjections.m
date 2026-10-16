function [muHat,info] = estimateMuOptimizedProjections(A,Y,varargin)
%ESTIMATEMUOPTIMIZEDPROJECTIONS Poisson tomography via optimized 1-D projections.
%
%   [muHat,info] = estimateMuOptimizedProjections(A,Y)
%   [muHat,info] = estimateMuOptimizedProjections(A,Y,'Name',value,...)
%
% MODEL
%   X_j ~ Poisson(mu_j), independent, j=1,...,p
%   Y   = A X, where A is m-by-p and Y is N-by-m (observations in rows).
%
% METHOD
%   Choose sparse nonnegative integer combinations
%
%       Z_h = c_h' Y,
%
%   with c_h in {0,1}^m.  If b_h=A'*c_h, then
%
%       Z_h = sum_j b_hj X_j,
%
%   is a one-dimensional compound Poisson random variable.  Its PMF is
%   computed exactly (up to ordinary floating-point arithmetic) by the
%   Panjer-type recursion
%
%       p(0) = exp(-sum_{j:b_j>0} mu_j),
%       p(z) = (1/z) sum_{k>=1} k*lambda_k*p(z-k),
%
%   where lambda_k=sum_{j:b_j=k} mu_j.
%
%   Candidate c-vectors are selected greedily by a D-optimal Fisher-
%   information criterion evaluated at a pilot mu.  The final estimator
%   maximizes the resulting composite log-likelihood using a constrained
%   trust-region method with mu>=0.
%
% IMPORTANT
%   No multidimensional Fourier inversion is used.  Every selected marginal
%   PMF is nonnegative by construction.  The criterion is a composite
%   likelihood because the projected variables are dependent.
%
% OPTIONS
%   Mu0                  Initial mu ([]: NNLS from A*mu ~= mean(Y)).
%   PilotMu              Pilot for projection design ([]: floored NNLS).
%   PilotFloor           Positive floor for automatic pilot (0.05).
%   MaxCombinationSize   Max number of measured components in c (3).
%   MaxCandidates        Max candidate projections retained/generated (500).
%   NumProjections       Number selected ([]: min(max(p,2*p),candidates)).
%   CandidateSeed        Seed for random candidate subsampling (1).
%   RequireGCDOne        Keep only projections with gcd(b_j>0)=1 (true).
%   DesignRidge          Ridge in D-optimal logdet criterion (1e-6).
%   DesignTailSigma      PMF support used in Fisher design, mean+s*sd (8).
%   DesignProbFloor      Ignore smaller probabilities in design FI (1e-14).
%   InitialRadius        Initial trust-region radius (1).
%   MaxRadius            Maximum trust-region radius (100).
%   MaxIter              Maximum trust-region iterations (100).
%   MaxSubproblemIter    SQP iterations for quadratic subproblem (100).
%   TolGradient          Projected-gradient tolerance (1e-6).
%   TolStep              Accepted-step tolerance (1e-10).
%   TolRadius            Minimum trust-region radius (1e-10).
%   EtaAccept            Minimum ratio for accepting a trial step (0.1).
%   EtaExpand            Ratio for trust-region expansion (0.75).
%   ShrinkFactor         Trust-region shrink factor (0.25).
%   ExpandFactor         Trust-region expansion factor (2).
%   BoundaryFraction     Expand when ||s|| >= fraction*Delta (0.8).
%   Display              'iter', 'final', or 'off' ('iter').
%   CheckDerivatives     Finite-difference derivative check (false).
%   ProjectionMatrix     [] or fixed H-by-m projection matrix. When supplied,
%                        candidate generation/Fisher D-optimal design are skipped.
%   DesignOnly           false. If true, stop after projection design and return.
%   TrimToIdentifiable   false. Stop greedy selection once both the jump-class
%                        and Fisher designs have full rank.
%   IncludeSignatureProjections true. Add the signature projections c=a_j'
%                        (one per column of A) to the candidates, even when
%                        they involve more than MaxCombinationSize sensors.
%                        With distinct nonzero columns, C=A' satisfies the
%                        jump-class rank condition (Proposition 4.2 of the paper).
%   SignatureRankFallback true. If the greedy selection does not reach
%                        jump-class rank p, append surviving signature
%                        projections until the rank is reached or none is left.
%   Contamination        0. Weight eps of a fixed background law in each scalar
%                        likelihood, p~(z)=(1-eps)p(z|mu)+eps*g(z), with g a
%                        geometric law with the sample mean of the projection.
%                        eps>0 makes counts that are impossible under Y=AX finite
%                        (weight zero); eps=0 is the exact composite likelihood. Use a
%                        tiny eps (e.g. 1e-10): with overdispersed data a large eps also
%                        downweights tail counts and biases the estimates downward.
%   IncludeSingleSensors false. If true, the single-sensor projections c=e_i of all
%                        observed sensors are added to the design (selected or given),
%                        with the same weight as every other projection; each sensor's
%                        own law then enters the criterion, which keeps the fitted
%                        sensor means close to the observed ones. Adding rows cannot
%                        lower the rank of the design.
%   MeanAnchorWeight     0. With w>0, the single-sensor projections c=e_i of all sensors
%                        are added to the selected ones (unless already selected), and the
%                        composite weights become 1/H_o for the H_o other projections and
%                        w/m_s for the m_s single-sensor ones. Each single-sensor term is a
%                        Poisson likelihood maximized at the sensor mean, so w>0 anchors the
%                        fit on the sensor means (as NNLS does, in Poisson deviance), while
%                        the multi-sensor terms still separate the classes. Useful when the
%                        path set cannot reproduce the joint count distributions exactly.
%   SubproblemSolver     'fast' (default): the trust-region subproblem (quadratic model,
%                        ball of radius Delta, bounds mu+s>=0) is solved by an active-set
%                        method with an exact ball-constrained solve on the free variables
%                        (eigendecomposition and secular equation), compared with the
%                        projected Cauchy step. 'fmincon': SQP via fmincon (slower; the
%                        previous behaviour). The acceptance test is the same, so both
%                        converge to the same maximum; individual steps may differ.
%   UseParallel          false. true evaluates the candidate projections and the
%                        multiple starts with parfor (Parallel Computing Toolbox;
%                        without it, or with false, the same loops run serially).
%                        Results are identical either way.
%   NumStarts            1. Number of initial points. With NumStarts>1 the projections
%                        are selected once and the trust-region method is run from:
%                        the NNLS fit (or Mu0), the floored pilot, and NumStarts-2
%                        random points pilot.*exp(MultiStartSpread*randn). The
%                        solution with the largest criterion is returned, and
%                        info.multistart reports how the solutions depend on the
%                        start (a check of nonconvexity).
%   MultiStartSpread     1 (log-scale spread of the random starts)
%   MultiStartSeed       1 (seed of the random starts; the global RNG is not used)
%   MultiStartTol        1e-2 (relative distance under which two solutions coincide)
%   ComputeSandwich      true. Return the Godambe (sandwich) covariance of the
%                        interior components of muHat.
%   ComputeCandidateFisherRanks true. Candidate-rank diagnostics; may be disabled
%                        for faster projection design without changing selection.
%
% OUTPUT HIGHLIGHTS
%   info.C                       H-by-m selected projection matrix.
%   info.B                       H-by-p jump/signature matrix C*A.
%   info.projectionWeights       Composite weights (sum to one).
%   info.designGains             Greedy D-optimal gains.
%   info.projectionDesignRank    Rank of stacked jump-class design.
%   info.fisherDesignRank        Rank of selected pilot Fisher information.
%   info.rankConditionHolds      projectionDesignRank==p.
%   info.nSignatureCandidates    Signature projections among the candidates.
%   info.nFallbackProjections    Signature projections appended by the fallback.
%   info.sandwichCovariance      Godambe covariance (NaN outside interior).
%   info.standardErrors          sqrt(diag(sandwichCovariance)).
%   info.interiorMask            Components treated as interior (muHat>tol).
%   info.sensitivityCondition    Condition number of -Hessian of the criterion at
%                                muHat, on the interior components, after scaling
%                                to unit diagonal (Inf if not positive definite):
%                                large values flag a poorly determined estimate.
%   info.nBoundary               Number of components of muHat at the bound 0.
%   info.contamination           eps used.
%   info.downweightedFraction    Share of (observation, projection) pairs whose
%                                weight (1-eps)p/p~ is below 0.5 at muHat.
%   info.multistart              (NumStarts>1) struct with the starts, the solutions,
%                                their criteria, relativeDistance=||mu_k-muHat||/||muHat||,
%                                fractionAtBest, nDistinctMaxima, maxRelativeDistance,
%                                criterionGap and criterionGapPerObservation.
%   info.criterion               Final composite log-likelihood.
%   info.firstOrderOptimality    Projected-gradient residual.
%   info.history                 Trust-region history table.
%
% Requires Optimization Toolbox: lsqnonneg and fmincon.

%% Validate model/data
if nargin<2
    error('estimateMuOptimizedProjections:Inputs','A and Y are required.');
end
validateattributes(A,{'numeric','logical'},{'2d','real','finite','nonempty'});
validateattributes(Y,{'numeric'},{'2d','real','finite','nonempty'});
A=full(double(A));
Y=full(double(Y));
[m,p]=size(A);
N=size(Y,1);
if size(Y,2)~=m
    error('estimateMuOptimizedProjections:Dimensions', ...
        'Y must have m=%d columns.',m);
end
if any(A(:)~=0 & A(:)~=1)
    error('estimateMuOptimizedProjections:BinaryA','A must be binary.');
end
if any(Y(:)<0) || any(abs(Y(:)-round(Y(:)))>1e-10)
    error('estimateMuOptimizedProjections:Counts', ...
        'Y must contain nonnegative integer counts.');
end
Y=round(Y);
if any(sum(A,1)==0)
    error('estimateMuOptimizedProjections:ZeroColumn', ...
        'A contains a zero column, so the corresponding mu is unobservable.');
end
if size(unique(A.','rows'),1)<p
    error('estimateMuOptimizedProjections:DuplicateColumns', ...
        'A has duplicate columns; individual rates are not identifiable.');
end
if exist('lsqnonneg','file')~=2 || exist('fmincon','file')~=2
    error('estimateMuOptimizedProjections:Toolbox', ...
        'Optimization Toolbox (lsqnonneg and fmincon) is required.');
end

%% Options
pos=@(x)isnumeric(x)&&isreal(x)&&isscalar(x)&&isfinite(x)&&x>0;
nonneg=@(x)isnumeric(x)&&isreal(x)&&isscalar(x)&&isfinite(x)&&x>=0;
posint=@(x)pos(x)&&x==floor(x);
bool=@(x)(islogical(x)||isnumeric(x))&&isscalar(x)&&any(x==[0 1]);
ip=inputParser;
ip.FunctionName=mfilename;
addParameter(ip,'Mu0',[],@(x)isempty(x)||isnumeric(x));
addParameter(ip,'PilotMu',[],@(x)isempty(x)||isnumeric(x));
addParameter(ip,'PilotFloor',0.05,pos);
addParameter(ip,'MaxCombinationSize',3,posint);
addParameter(ip,'MaxCandidates',500,posint);
addParameter(ip,'NumProjections',[],@(x)isempty(x)||posint(x));
addParameter(ip,'CandidateSeed',1,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=0&&x==floor(x));
addParameter(ip,'RequireGCDOne',true,bool);
addParameter(ip,'DesignRidge',1e-6,pos);
addParameter(ip,'DesignTailSigma',8,pos);
addParameter(ip,'DesignProbFloor',1e-14,pos);
addParameter(ip,'InitialRadius',1,pos);
addParameter(ip,'MaxRadius',100,pos);
addParameter(ip,'MaxIter',100,posint);
addParameter(ip,'MaxSubproblemIter',100,posint);
addParameter(ip,'TolGradient',1e-6,pos);
addParameter(ip,'TolStep',1e-10,pos);
addParameter(ip,'TolRadius',1e-10,pos);
addParameter(ip,'EtaAccept',0.1,pos);
addParameter(ip,'EtaExpand',0.75,pos);
addParameter(ip,'ShrinkFactor',0.25,pos);
addParameter(ip,'ExpandFactor',2,pos);
addParameter(ip,'BoundaryFraction',0.8,pos);
addParameter(ip,'Display','iter',@(x)ischar(x)||isstring(x));
addParameter(ip,'CheckDerivatives',false,bool);
addParameter(ip,'ProjectionMatrix',[],@(x)isempty(x)||isnumeric(x)||islogical(x));
addParameter(ip,'ProjectionWeights',[],@(x)isempty(x)||(isnumeric(x)&&all(x(:)>=0)));
addParameter(ip,'DesignOnly',false,bool);
addParameter(ip,'TrimToIdentifiable',false,bool);
addParameter(ip,'ComputeCandidateFisherRanks',true,bool);
addParameter(ip,'IncludeSignatureProjections',true,bool);
addParameter(ip,'SignatureRankFallback',true,bool);
addParameter(ip,'ComputeSandwich',true,bool);
addParameter(ip,'UseParallel',false,bool);
addParameter(ip,'SubproblemSolver','fast',@(x)any(strcmpi(x,{'fast','fmincon'})));
addParameter(ip,'MeanAnchorWeight',0,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=0);
addParameter(ip,'IncludeSingleSensors',false,@(x)islogical(x)||isnumeric(x));
addParameter(ip,'NumStarts',1,@(x)isnumeric(x)&&isscalar(x)&&x>=1&&x==fix(x));
addParameter(ip,'MultiStartSpread',1,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(ip,'MultiStartSeed',1,@(x)isnumeric(x)&&isscalar(x)&&x>=0&&x==fix(x));
addParameter(ip,'MultiStartTol',1e-2,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(ip,'Contamination',0,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=0&&x<1);
parse(ip,varargin{:});
o=ip.Results;
o.Display=validatestring(o.Display,{'off','final','iter'});
o.RequireGCDOne=logical(o.RequireGCDOne);
o.CheckDerivatives=logical(o.CheckDerivatives);
o.DesignOnly=logical(o.DesignOnly);
o.TrimToIdentifiable=logical(o.TrimToIdentifiable);
o.ComputeCandidateFisherRanks=logical(o.ComputeCandidateFisherRanks);
o.IncludeSignatureProjections=logical(o.IncludeSignatureProjections);
o.SignatureRankFallback=logical(o.SignatureRankFallback);
o.ComputeSandwich=logical(o.ComputeSandwich);
if o.NumStarts>1 && ~o.DesignOnly
    % Multi-start: same projections, several initial points, best criterion kept.
    % (A design-only call selects projections and fits nothing: NumStarts is
    % irrelevant there and is ignored.)
    [muHat,info]=multiStartEstimate(A,Y,varargin,o);
    return
end
if o.MaxCombinationSize>m
    o.MaxCombinationSize=m;
end
if isempty(o.ProjectionMatrix) && o.MaxCandidates<m
    error('estimateMuOptimizedProjections:MaxCandidates', ...
        'MaxCandidates must be at least m=%d so singleton projections are included.',m);
end
if ~(o.EtaAccept<o.EtaExpand && o.EtaExpand<1 && o.ShrinkFactor<1 && ...
        o.ExpandFactor>1 && o.BoundaryFraction<=1 && ...
        o.InitialRadius<=o.MaxRadius && o.TolRadius<o.InitialRadius)
    error('estimateMuOptimizedProjections:Options','Inconsistent trust-region options.');
end

%% Pilot and candidate projection design
meanY=mean(Y,1).';
% Avoid a second NNLS solve when the caller already supplies both the pilot
% and the optimization starting point (the Parma pipeline does this per stratum).
needNNLS=isempty(o.PilotMu) || (~o.DesignOnly && isempty(o.Mu0));
if needNNLS
    muNNLS=lsqnonneg(A,meanY,optimset('Display','off'));
else
    muNNLS=nan(p,1);
end
if isempty(o.PilotMu)
    muPilot=max(muNNLS,o.PilotFloor);
else
    validateattributes(o.PilotMu,{'numeric'},{'vector','real','finite','nonnegative','numel',p});
    muPilot=full(double(o.PilotMu(:)));
end

designTic=tic;
fixedProjection=~isempty(o.ProjectionMatrix);
if ~fixedProjection
    Ccand=generateCandidateMatrix(m,o.MaxCombinationSize,o.MaxCandidates,o.CandidateSeed);
    % Signature projections c=a_j' (Proposition 4.2 of the paper).  They may
    % involve more than MaxCombinationSize sensors and are added on top of
    % the MaxCandidates budget.
    Csig=double(A.'>0);
    if o.IncludeSignatureProjections
        Ccand=unique([double(Ccand);Csig],'rows','stable');
    end
    Bcand=Ccand*A; % candidates-by-p
    keep=any(Bcand>0,2);
    if o.RequireGCDOne
        for i=1:size(Bcand,1)
            if keep(i)
                bp=Bcand(i,Bcand(i,:)>0);
                keep(i)=vectorGCD(bp)==1;
            end
        end
    end
    Ccand=Ccand(keep,:);
    Bcand=Bcand(keep,:);
    if isempty(Ccand)
        error('estimateMuOptimizedProjections:Candidates','No valid projection candidates remain.');
    end
    % Same b gives the same scalar model distribution; keep one representative.
    [~,ia]=unique(Bcand,'rows','stable');
    Ccand=Ccand(ia,:);
    Bcand=Bcand(ia,:);
    Q=size(Ccand,1);

    Fcand=zeros(p,p,Q);
    fiTrace=zeros(Q,1);
    fiRank=zeros(Q,1);
    % Compute observed projected maxima in BLAS-backed blocks.  This is exact
    % and much faster than Q separate Y*c multiplications.
    zDataMaxAll=zeros(Q,1);
    designBlock=64;
    for q0=1:designBlock:Q
        qq=q0:min(Q,q0+designBlock-1);
        Zblk=Y*Ccand(qq,:).';
        zDataMaxAll(qq)=max(Zblk,[],1).';
    end
    parWorkers=0; if o.UseParallel, parWorkers=Inf; end
    tailSigma=o.DesignTailSigma; probFloor=o.DesignProbFloor; wantRanks=o.ComputeCandidateFisherRanks;
    parfor (q=1:Q, parWorkers)
        F=projectionFisherInformation(Bcand(q,:).',muPilot,zDataMaxAll(q),tailSigma,probFloor);
        Fcand(:,:,q)=F;
        fiTrace(q)=trace(F);
        if wantRanks
            fiRank(q)=rank(F,1e-10*max(1,norm(F,2)));
        else
            fiRank(q)=NaN;
        end
    end
    informative=fiTrace>1e-14;
    Ccand=Ccand(informative,:);
    Bcand=Bcand(informative,:);
    Fcand=Fcand(:,:,informative);
    fiTrace=fiTrace(informative);
    fiRank=fiRank(informative);
    Q=size(Ccand,1);
    if Q==0
        error('estimateMuOptimizedProjections:Information','All projection candidates have negligible pilot information.');
    end

    isSignature=ismember(Ccand,Csig,'rows');
    if isempty(o.NumProjections)
        Hwanted=min(Q,max(p,min(2*p,20)));
    else
        Hwanted=min(Q,o.NumProjections);
    end
    [selected,designGains,designScores,Fselected]=greedyDOptimal( ...
        Fcand,Hwanted,o.DesignRidge,Bcand,p,o.TrimToIdentifiable);
    selected=selected(:);
    nFallback=0;
    if o.SignatureRankFallback
        rankNow=rank(jumpClassDesign(Bcand(selected,:)));
        if rankNow<p
            cand=find(isSignature & ~ismember((1:Q).',selected));
            for q=cand.'
                rankTrial=rank(jumpClassDesign(Bcand([selected;q],:)));
                if rankTrial>rankNow
                    selected=[selected;q]; %#ok<AGROW>
                    designGains=[designGains(:);NaN];
                    designScores=[designScores(:);NaN];
                    Fselected=Fselected+Fcand(:,:,q);
                    rankNow=rankTrial;
                    nFallback=nFallback+1;
                    if rankNow>=p, break; end
                end
            end
        end
    end
    C=Ccand(selected,:);
    B=Bcand(selected,:);
    nSignatureCandidates=nnz(isSignature);
else
    C=full(double(o.ProjectionMatrix));
    if size(C,2)~=m || isempty(C) || any(~isfinite(C(:))) || ...
            any(C(:)<0) || any(abs(C(:)-round(C(:)))>1e-12) || any(sum(C,2)==0)
        error('estimateMuOptimizedProjections:ProjectionMatrix', ...
            'ProjectionMatrix must be a nonempty H-by-m nonnegative integer matrix with nonzero rows.');
    end
    C=round(C);
    B=C*A;
    keep=any(B>0,2);
    C=C(keep,:); B=B(keep,:);
    if isempty(C)
        error('estimateMuOptimizedProjections:ProjectionMatrix','ProjectionMatrix has no informative rows for A.');
    end
    Hproj=size(C,1);
    Q=Hproj;
    Ccand=C; Bcand=B;
    selected=(1:Hproj).';
    designGains=nan(Hproj,1);
    designScores=nan(Hproj,1);
    fiTrace=nan(Hproj,1);
    fiRank=nan(Hproj,1);
    Fcand=[];
    Fselected=nan(p);
    nFallback=0;
    nSignatureCandidates=0;
end
if logical(o.IncludeSingleSensors) && ~o.DesignOnly
    % Every observed sensor's own projection, with ordinary (equal) weight.
    isSingle0=sum(C~=0,2)==1 & sum(C,2)==1;
    have0=false(m,1);
    if any(isSingle0), [~,col0]=max(C(isSingle0,:),[],2); have0(col0)=true; end
    observed0=find(any(A>0,2)); add0=observed0(~have0(observed0));
    if ~isempty(add0)
        Cadd0=zeros(numel(add0),m); Cadd0(sub2ind(size(Cadd0),(1:numel(add0)).',add0(:)))=1;
        C=[C;Cadd0]; B=[B;Cadd0*A];
    end
end
Hproj=size(C,1);
projectionWeights=ones(Hproj,1)/Hproj;
if ~isempty(o.ProjectionWeights)
    % Weights of the projections (with ProjectionMatrix), e.g. from estimateMuWeightedProjections.
    assert(numel(o.ProjectionWeights)==Hproj,'ProjectionWeights must have one entry per projection.');
    projectionWeights=double(o.ProjectionWeights(:))/max(sum(o.ProjectionWeights(:)),realmin);
end
if o.MeanAnchorWeight>0 && ~o.DesignOnly
    % Mean anchoring: append the single-sensor projections and reweight.
    isSingle=sum(C~=0,2)==1 & sum(C,2)==1;
    observed=find(any(A>0,2));
    have=false(m,1);
    if any(isSingle), [~,col]=max(C(isSingle,:),[],2); have(col)=true; end
    add=observed(~have(observed));
    if ~isempty(add)
        Cadd=zeros(numel(add),m); Cadd(sub2ind(size(Cadd),(1:numel(add)).',add(:)))=1;
        C=[C;Cadd]; B=[B;Cadd*A];
        isSingle=[isSingle;true(numel(add),1)];
        Hproj=size(C,1);
    end
    nS=nnz(isSingle); nO=Hproj-nS;
    projectionWeights=zeros(Hproj,1);
    if nO>0, projectionWeights(~isSingle)=1/nO; end
    projectionWeights(isSingle)=o.MeanAnchorWeight/max(nS,1);
    projectionWeights=projectionWeights/sum(projectionWeights);
end
designRuntimeSeconds=toc(designTic);

%% Projection identifiability diagnostics and sufficient statistics
R=jumpClassDesign(B);
projectionDesignRank=rank(R);
projectionDesignCondition=matrixCondition(R,p);
if fixedProjection
    fisherDesignRank=NaN;
    fisherDesignCondition=NaN;
else
    fisherDesignRank=rank(Fselected,1e-10*max(1,norm(Fselected,2)));
    fisherDesignCondition=matrixCondition(Fselected,p);
end
if projectionDesignRank<p
    warning('estimateMuOptimizedProjections:ProjectionRank', ...
        ['Selected scalar projections have jump-class design rank %d < p=%d. ', ...
         'Increase NumProjections, MaxCombinationSize, or MaxCandidates.'], ...
        projectionDesignRank,p);
end
if ~isnan(fisherDesignRank) && fisherDesignRank<p
    warning('estimateMuOptimizedProjections:FisherRank', ...
        ['Selected pilot Fisher information has rank %d < p=%d. ', ...
         'The composite criterion may not identify every rate locally.'], ...
        fisherDesignRank,p);
end

if o.DesignOnly
    muHat=muPilot;
    info=struct();
    info.C=C;
    info.B=B;
    info.projectionWeights=projectionWeights;
    info.nCandidates=Q;
    info.nSelectedProjections=Hproj;
    info.selectedCandidateIndices=selected(:);
    info.designGains=designGains(:);
    info.designScores=designScores(:);
    info.candidateC=Ccand;
    info.candidateB=Bcand;
    info.candidateFisherTrace=fiTrace;
    info.candidateFisherRank=fiRank;
    info.selectedPilotFisher=Fselected;
    info.projectionDesignMatrix=R;
    info.projectionDesignRank=projectionDesignRank;
    info.rankConditionHolds=projectionDesignRank>=p;
    info.nSignatureCandidates=nSignatureCandidates;
    info.nFallbackProjections=nFallback;
    info.projectionDesignCondition=projectionDesignCondition;
    info.fisherDesignRank=fisherDesignRank;
    info.fisherDesignCondition=fisherDesignCondition;
    info.muNNLS=muNNLS;
    info.muPilot=muPilot;
    info.options=o;
    info.designRuntimeSeconds=designRuntimeSeconds;
    info.designOnly=true;
    info.fixedProjectionMatrix=fixedProjection;
    return
end

proj=repmat(struct('c',[],'b',[],'zValues',[],'counts',[],'maxZ',0,'eps',0,'logG',[]),Hproj,1);
Zselected=Y*C.';
for h=1:Hproj
    z=Zselected(:,h);
    [zv,~,iz]=unique(z);
    cnt=accumarray(iz(:),1);
    proj(h).c=C(h,:).';
    proj(h).b=B(h,:).';
    proj(h).zValues=zv(:);
    proj(h).counts=cnt(:);
    proj(h).maxZ=max(zv);
    proj(h).eps=o.Contamination;
    % Background law: geometric with the sample mean (fixed, independent of mu).
    mbar=max(mean(z),1); th=mbar/(1+mbar);
    proj(h).logG=log(1-th)+double(zv(:))*log(th);
end

%% Initial point
if isempty(o.Mu0)
    if any(~isfinite(muNNLS)), muNNLS=lsqnonneg(A,meanY,optimset('Display','off')); end
    mu=muNNLS;
    startName='nonnegative least squares';
else
    validateattributes(o.Mu0,{'numeric'},{'vector','real','finite','nonnegative','numel',p});
    mu=full(double(o.Mu0(:)));
    startName='user supplied';
end
v=evaluateComposite(mu,proj,projectionWeights,p,true);
if ~v.valid
    % The floored pilot is strictly positive and is typically numerically safer.
    mu=max(muPilot,o.PilotFloor);
    startName='positive pilot fallback';
    v=evaluateComposite(mu,proj,projectionWeights,p,true);
end
if ~v.valid
    error('estimateMuOptimizedProjections:InitialPoint', ...
        ['Could not evaluate a finite projected composite likelihood at the initial point. ', ...
         'The PMF is evaluated in log scale, so this usually indicates a structurally ', ...
         'impossible observed projected count or a numerically extreme derivative ratio.']);
end
mu0=mu;

%% Trust-region optimization
innerOptions=optimoptions('fmincon','Algorithm','sqp','Display','off', ...
    'SpecifyObjectiveGradient',true, ...
    'SpecifyConstraintGradient',true, ...
    'MaxIterations',o.MaxSubproblemIter, ...
    'MaxFunctionEvaluations',max(1000,150*p), ...
    'OptimalityTolerance',1e-10,'ConstraintTolerance',1e-10,'StepTolerance',1e-12);

Delta=o.InitialRadius;
exitflag=0;
message='Maximum number of trust-region iterations reached.';
iterDone=0;
subproblemCalls=0;
derivativeCheck=struct();

% iteration, criterion, projectedGradient, radiusBefore, radiusAfter,
% stepNorm, predictedIncrease, actualIncrease, rho, accepted
history=nan(o.MaxIter+1,10);
muHistory=nan(p,o.MaxIter+1);
muHistory(:,1)=mu;
pg=projectedGradientResidual(mu,v.g);
history(1,:)=[0,v.ell,pg,Delta,Delta,0,NaN,NaN,NaN,false];

if o.CheckDerivatives
    derivativeCheck=checkCompositeDerivatives(mu,proj,projectionWeights,p);
end

if strcmp(o.Display,'iter')
    fprintf('Optimized-projection compound-Poisson trust-region estimator\n');
    fprintf('Candidates: %d; selected projections: %d; parameters: %d\n',Q,Hproj,p);
    fprintf('Projection design rank: %d/%d; Fisher design rank: %d/%d\n', ...
        projectionDesignRank,p,fisherDesignRank,p);
    fprintf('Start: %s; criterion %.10g\n',startName,v.ell);
    fprintf(' iter       criterion          projPG       radius       rho    accepted\n');
    fprintf('%5d  %16.8g  %12.4g  %10.4g  %8s  %8d\n',0,v.ell,pg,Delta,'-',0);
end

t0=tic;
for k=1:o.MaxIter
    iterDone=k;
    pg=projectedGradientResidual(mu,v.g);
    if pg<=o.TolGradient
        exitflag=1;
        message='Projected first-order optimality tolerance satisfied.';
        iterDone=k-1;
        break
    end
    if Delta<=o.TolRadius
        exitflag=2;
        message='Trust-region radius tolerance satisfied.';
        iterDone=k-1;
        break
    end

    oldMu=mu;
    oldEll=v.ell;
    oldDelta=Delta;
    if strcmpi(o.SubproblemSolver,'fast')
        [s,subflag,nCalls]=solveQuadraticModelFast(mu,v.g,v.H,Delta);
    else
        [s,subflag,nCalls]=solveQuadraticModel(mu,v.g,v.H,Delta,innerOptions);
    end
    subproblemCalls=subproblemCalls+nCalls;
    stepNorm=norm(s);
    pred=v.g.'*s+0.5*s.'*v.H*s;
    accepted=false;
    actual=NaN;
    rho=-Inf;

    if isempty(s) || any(~isfinite(s)) || subflag<=-2 || ...
            ~(isfinite(pred) && pred>100*eps*max(1,abs(oldEll)))
        Delta=max(o.ShrinkFactor*Delta,o.TolRadius);
    else
        trialMu=max(0,mu+s);
        % The criterion alone decides acceptance; gradient and Hessian are
        % computed only for accepted points (rejected steps skip that work).
        trial=evaluateComposite(trialMu,proj,projectionWeights,p,false);
        if trial.valid
            actual=trial.ell-oldEll;
            rho=actual/pred;
            if actual>0 && rho>=o.EtaAccept
                trialFull=evaluateComposite(trialMu,proj,projectionWeights,p,true);
                if trialFull.valid
                    mu=trialMu;
                    v=trialFull;
                    accepted=true;
                else
                    % Derivatives not finite: treat as an invalid trial, as before.
                    actual=NaN; rho=-Inf;
                end
            end
        end

        if ~isfinite(rho) || rho<o.EtaAccept
            Delta=max(o.ShrinkFactor*oldDelta,o.TolRadius);
        elseif rho>=o.EtaExpand && stepNorm>=o.BoundaryFraction*oldDelta
            Delta=min(o.ExpandFactor*oldDelta,o.MaxRadius);
        else
            Delta=oldDelta;
        end
    end

    pgNew=projectedGradientResidual(mu,v.g);
    history(k+1,:)=[k,v.ell,pgNew,oldDelta,Delta,stepNorm,pred,actual,rho,accepted];
    muHistory(:,k+1)=mu;

    if strcmp(o.Display,'iter')
        fprintf('%5d  %16.8g  %12.4g  %10.4g  %8.3g  %8d\n', ...
            k,v.ell,pgNew,Delta,rho,accepted);
    end

    if accepted && norm(mu-oldMu)<=o.TolStep*(1+norm(oldMu))
        exitflag=3;
        message='Accepted-step tolerance satisfied.';
        break
    end
end
runtimeSeconds=toc(t0);
if iterDone==o.MaxIter && exitflag==0
    message='Maximum number of trust-region iterations reached.';
end

% Trim history. If convergence was detected before executing iteration k,
% the last filled row is iterDone+1 with iteration number iterDone.
lastRow=find(isfinite(history(:,1)),1,'last');
history=history(1:lastRow,:);
muHistory=muHistory(:,1:lastRow);
muHat=mu;
finalPG=projectedGradientResidual(muHat,v.g);

histTable=array2table(history,'VariableNames', ...
    {'Iteration','Criterion','ProjectedGradient','RadiusBefore','RadiusAfter', ...
     'StepNorm','PredictedIncrease','ActualIncrease','Rho','Accepted'});

%% Output
info=struct();
info.criterion=v.ell;
info.meanCriterion=v.ell/N;
info.gradient=v.g;
info.hessian=v.H;
info.firstOrderOptimality=finalPG;
info.exitflag=exitflag;
info.converged=exitflag>0;
info.message=message;
info.nIterations=max(history(:,1));
info.subproblemCalls=subproblemCalls;
info.runtimeSeconds=runtimeSeconds;
info.designRuntimeSeconds=designRuntimeSeconds;
info.fixedProjectionMatrix=fixedProjection;
info.history=histTable;
info.muHistory=muHistory;
info.mu0=mu0;
info.muNNLS=muNNLS;
info.muPilot=muPilot;
info.startName=startName;
info.C=C;
info.B=B;
info.projectionWeights=projectionWeights;
info.nCandidates=Q;
info.nSelectedProjections=Hproj;
info.selectedCandidateIndices=selected(:);
info.designGains=designGains(:);
info.designScores=designScores(:);
info.candidateC=Ccand;
info.candidateB=Bcand;
info.candidateFisherTrace=fiTrace;
info.candidateFisherRank=fiRank;
info.selectedPilotFisher=Fselected;
info.projectionDesignMatrix=R;
info.projectionDesignRank=projectionDesignRank;
info.rankConditionHolds=projectionDesignRank>=p;
info.nSignatureCandidates=nSignatureCandidates;
info.nFallbackProjections=nFallback;
info.projectionDesignCondition=projectionDesignCondition;
info.fisherDesignRank=fisherDesignRank;
info.fisherDesignCondition=fisherDesignCondition;
info.projections=proj;
if o.ComputeSandwich
    [info.sandwichCovariance,info.standardErrors,info.interiorMask]= ...
        compositeSandwich(muHat,proj,projectionWeights,Zselected,v.H);
else
    info.sandwichCovariance=[]; info.standardErrors=[]; info.interiorMask=[];
end
[info.sensitivityCondition,info.nBoundary]=sensitivityConditioning(muHat,v.H);
info.contamination=o.Contamination;
info.downweightedFraction=downweightedShare(muHat,proj);
info.derivativeCheck=derivativeCheck;
info.options=o;
info.pmfEvaluation='log-domain Panjer recursion';
info.derivativeImplementation='vectorized exact shift identities';
info.minLogProbability=v.minLogProbability;

if strcmp(o.Display,'final') || strcmp(o.Display,'iter')
    fprintf('%s\n',message);
    fprintf('Final composite criterion = %.10g; projected gradient = %.4g; elapsed = %.3f s\n', ...
        info.criterion,info.firstOrderOptimality,info.runtimeSeconds);
end
end

% =========================================================================
function C=generateCandidateMatrix(m,sMax,maxCandidates,seed)
% Always include singleton measured components. Add all small supports when
% possible; otherwise fill the remaining budget by reproducible random supports.
C=eye(m);
if sMax<=1 || size(C,1)>=maxCandidates
    % (fixed: the original indexing failed when sMax<=1 and m<maxCandidates)
    C=C(1:min(maxCandidates,size(C,1)),:);
    return
end
rs=RandStream('mt19937ar','Seed',seed);
for s=2:sMax
    remaining=maxCandidates-size(C,1);
    if remaining<=0, break; end
    logN=sum(log((m-s+1):m))-sum(log(1:s));
    nComb=round(exp(logN));
    if isfinite(nComb) && nComb<=remaining && nComb<=2e5
        idx=nchoosek(1:m,s);
        block=zeros(size(idx,1),m);
        for r=1:size(idx,1)
            block(r,idx(r,:))=1;
        end
        C=[C;block]; %#ok<AGROW>
    else
        block=zeros(0,m);
        attempts=0;
        maxAttempts=max(1000,50*remaining);
        while size(block,1)<remaining && attempts<maxAttempts
            attempts=attempts+1;
            perm=randperm(rs,m);
            supp=sort(perm(1:s));
            row=zeros(1,m); row(supp)=1;
            if ~ismember(row,C,'rows') && ~ismember(row,block,'rows')
                block=[block;row]; %#ok<AGROW>
            end
        end
        C=[C;block]; %#ok<AGROW>
        if size(C,1)>=maxCandidates, break; end
    end
end
if size(C,1)>maxCandidates
    C=C(1:maxCandidates,:);
end
end

% =========================================================================
function g=vectorGCD(v)
v=round(abs(v(:)));
v=v(v>0);
if isempty(v), g=0; return; end
g=v(1);
for i=2:numel(v)
    g=gcd(g,v(i));
end
end

% =========================================================================
function F=projectionFisherInformation(b,mu,zDataMax,tailSigma,probFloor)
% Fisher information evaluated from a log-domain compound-Poisson PMF.
% The score ratios for all z and all parameters are formed in one matrix;
% this is algebraically identical to the scalar loops but much faster in MATLAB.
b=round(b(:));
mu=mu(:);
meanZ=b.'*mu;
varZ=(b.^2).'*mu;
maxJump=max(b);
zMax=max([double(zDataMax),maxJump,ceil(meanZ+tailSigma*sqrt(max(varZ,0))+maxJump)]);
logP=compoundPoissonLogPmf(b,mu,zMax);
pnum=numel(mu);

valid=isfinite(logP) & logP>log(probFloor);
if ~any(valid)
    F=zeros(pnum);
    return
end

zVals=find(valid)-1;
lp=logP(valid);
idx=bsxfun(@minus,zVals,b.');
ok=idx>=0;
idxSafe=idx+1;
idxSafe(~ok)=1;
lpShift=logP(idxSafe);
lr=bsxfun(@minus,lpShift,lp);
ratio=zeros(size(idx));
if any(ok(:))
    lrGood=lr(ok);
    if any(lrGood>log(realmax))
        F=zeros(pnum); return
    end
    ratio(ok)=exp(lrGood);
end
score=ratio-1;
pv=exp(lp);
W=bsxfun(@times,score,sqrt(pv));
F=W.'*W;
F=0.5*(F+F.');
if any(~isfinite(F(:)))
    F=zeros(pnum);
end
end

% =========================================================================
function [selected,gains,scores,Fsum]=greedyDOptimal(Fcand,Hwanted,ridge,Bcand,pTarget,stopWhenIdentified)
p=size(Fcand,1);
Q=size(Fcand,3);
if nargin<4, Bcand=[]; end
if nargin<5 || isempty(pTarget), pTarget=p; end
if nargin<6, stopWhenIdentified=false; end
selected=zeros(Hwanted,1);
gains=zeros(Hwanted,1);
scores=zeros(Hwanted,1);
used=false(Q,1);
Fsum=zeros(p);
current=logdetSPD(ridge*eye(p)+Fsum);
for h=1:Hwanted
    best=-Inf;
    bestIdx=0;
    base=ridge*eye(p)+Fsum;          % same additions, in the same order, as before
    for q=1:Q
        if used(q), continue; end
        sc=logdetSPD(base+Fcand(:,:,q));
        if sc>best
            best=sc;
            bestIdx=q;
        end
    end
    if bestIdx==0
        selected=selected(1:h-1);
        gains=gains(1:h-1);
        scores=scores(1:h-1);
        return
    end
    selected(h)=bestIdx;
    gains(h)=best-current;
    scores(h)=best;
    used(bestIdx)=true;
    Fsum=Fsum+Fcand(:,:,bestIdx);
    current=best;

    % In fast mode, stop the expensive greedy search as soon as the current
    % prefix is already locally identifying.  No suffix is selected and then
    % discarded, so both design time and later likelihood time are reduced.
    if stopWhenIdentified && ~isempty(Bcand)
        Rprefix=jumpClassDesign(Bcand(selected(1:h),:));
        rJump=rank(Rprefix);
        rFish=rank(Fsum,1e-10*max(1,norm(Fsum,2)));
        if rJump>=pTarget && rFish>=pTarget
            selected=selected(1:h);
            gains=gains(1:h);
            scores=scores(1:h);
            return
        end
    end
end
end

% =========================================================================
function val=logdetSPD(M)
M=0.5*(M+M.');
[R,flag]=chol(M);
if flag==0
    val=2*sum(log(diag(R)));
else
    e=eig(M);
    e=max(real(e),realmin);
    val=sum(log(e));
end
end

% =========================================================================
function R=jumpClassDesign(B)
% Rows encode the identifiable aggregate intensities lambda_{h,k}.
[H,p]=size(B);
R=zeros(0,p);
for h=1:H
    bh=B(h,:);
    for k=1:max(bh)
        row=double(bh==k);
        if any(row)
            R=[R;row]; %#ok<AGROW>
        end
    end
end
R=unique(R,'rows','stable');
end

% =========================================================================
function c=matrixCondition(M,p)
if isempty(M)
    c=Inf; return
end
s=svd(M,'econ');
if numel(s)<p || s(p)<=max(size(M))*eps(max(s))
    c=Inf;
else
    c=s(1)/s(p);
end
end

% =========================================================================
function logP=compoundPoissonLogPmf(b,mu,zMax)
% Exact Panjer-type coefficient recursion in log probability.
%
%   p(0)=exp(-sum_k lambda_k)
%   z p(z)=sum_k k lambda_k p(z-k)
%
% Taking log-sum-exp in the second relation prevents the p(0) underflow
% that occurs for high-volume traffic data. No approximation is introduced.
b=round(b(:));
mu=mu(:);
K=max(b);
logP=-inf(zMax+1,1);
if K<=0
    logP(1)=0;
    return
end

lambda=zeros(K,1);
for j=1:numel(mu)
    if b(j)>0
        lambda(b(j))=lambda(b(j))+mu(j);
    end
end
logP(1)=-sum(lambda);

% log(k)+log(lambda_k) and log(z) do not depend on the step: computing them
% once gives the same values, added in the same order, as in the loop.
ksAll=(1:K).';
logC=log(double(ksAll))+log(lambda);
posK=lambda>0;
logZ=log(double((1:zMax).'));
if K==1
    % Single jump size: terms has one element t, and m+log(sum(exp(t-m)))
    % equals t exactly, so the update reduces to t-log(z).
    c1=logC(1);
    if posK(1)
        for z=1:zMax
            logP(z+1)=(c1+logP(z))-logZ(z);
        end
    end
    return
end
persistent useScalar
if isempty(useScalar)
    % MATLAB compiles scalar loops (JIT), which makes them faster than small
    % vector operations; Octave does not, so it keeps the vector version.
    useScalar=exist('OCTAVE_VERSION','builtin')==0;
end
if K<=4 && useScalar
    % Scalar version for few jump sizes (the usual case): the same terms,
    % maximum and sum in the same order, without small-vector overheads.
    for z=1:zMax
        kmax=min(K,z);
        m=-Inf; any_=false;
        for k=1:kmax
            if posK(k)
                t=logC(k)+logP(z-k+1);
                if t>-Inf
                    if ~any_ || t>m, m=t; end
                    any_=true;
                end
            end
        end
        if ~any_, continue; end
        acc=0;
        for k=1:kmax
            if posK(k)
                t=logC(k)+logP(z-k+1);
                if t>-Inf, acc=acc+exp(t-m); end
            end
        end
        logP(z+1)=m+log(acc)-logZ(z);
    end
    return
end
for z=1:zMax
    kmax=min(K,z);
    ks=ksAll(1:kmax);
    prev=logP(z-ks+1);
    use=posK(1:kmax) & isfinite(prev);
    if ~any(use)
        continue
    end
    terms=logC(ks(use))+prev(use);
    m=max(terms);
    logP(z+1)=m+log(sum(exp(terms-m)))-logZ(z);
end
end

% =========================================================================
function v=evaluateComposite(mu,proj,weights,pnum,needDerivatives)
% Evaluate the projected composite likelihood in the log-PMF domain.
% Exact shift identities are vectorized over all distinct observed projected
% counts.  This removes the previous z-by-parameter-by-parameter MATLAB loops
% without changing the composite likelihood, gradient, or Hessian.
mu=mu(:);
Hproj=numel(proj);
ell=0;
g=zeros(pnum,1);
Hess=zeros(pnum);
valid=true;
minLogProb=Inf;
logRealMax=log(realmax);

for h=1:Hproj
    b=round(proj(h).b(:));
    zVals=double(proj(h).zValues(:));
    counts=double(proj(h).counts(:));
    logP=compoundPoissonLogPmf(b,mu,proj(h).maxZ);
    wh=weights(h);

    lpz=logP(zVals+1);
    if isfield(proj,'eps') && ~isempty(proj(h).eps) && proj(h).eps>0
        % Contaminated likelihood p~=(1-eps)p+eps*g (see option Contamination).
        [lh,gh,Hh,ok]=contaminatedTerms(logP,lpz,b,zVals,wh*counts,proj(h).eps, ...
            proj(h).logG,pnum,needDerivatives,logRealMax);
        if ~ok, valid=false; continue; end
        minLogProb=min(minLogProb,min(lh.lmix));
        ell=ell+lh.ell;
        if needDerivatives, g=g+gh; Hess=Hess+Hh; end
        continue
    end
    minLogProb=min(minLogProb,min(lpz));
    if any(~isfinite(lpz))
        valid=false;
        continue
    end

    a=wh*counts;
    ell=ell+sum(a.*lpz);
    if ~needDerivatives
        continue
    end

    % ratio1(r,j)=p(z_r-b_j)/p(z_r).  Matrix indexing handles every
    % observed z and parameter in one operation.
    idx1=bsxfun(@minus,zVals,b.');
    ok1=idx1>=0;
    idxSafe=idx1+1;
    idxSafe(~ok1)=1;
    lpShift=logP(idxSafe);
    lr1=bsxfun(@minus,lpShift,lpz);
    ratio1=zeros(size(idx1));
    if any(ok1(:))
        lrGood=lr1(ok1);
        if any(lrGood>logRealMax)
            valid=false;
            continue
        end
        ratio1(ok1)=exp(lrGood);
    end
    score=ratio1-1;
    if any(~isfinite(score(:)))
        valid=false;
        continue
    end
    g=g+score.'*a;

    % First Hessian term depends only on the integer shift b_j+b_k.
    % There are at most 2*max(b)+1 distinct shifts, usually only a handful,
    % so accumulate each shift once rather than once for every (j,k,z).
    maxShift=2*max(b);
    shiftWeighted=zeros(maxShift+1,1);
    for sh=0:maxShift
        idx=zVals-sh;
        ok=idx>=0;
        if ~any(ok), continue; end
        lr=logP(idx(ok)+1)-lpz(ok);
        if any(lr>logRealMax)
            valid=false;
            break
        end
        shiftWeighted(sh+1)=sum(a(ok).*exp(lr));
    end
    if ~valid, continue; end
    sumShift=bsxfun(@plus,b,b.')+1;
    firstTerm=shiftWeighted(sumShift);
    secondTerm=ratio1.'*bsxfun(@times,ratio1,a);
    Hterm=firstTerm-secondTerm;
    if any(~isfinite(Hterm(:)))
        valid=false;
        continue
    end
    Hess=Hess+Hterm;
end

if ~valid || ~isfinite(ell)
    ell=-Inf;
    if needDerivatives
        g(:)=NaN;
        Hess(:)=NaN;
    end
elseif needDerivatives
    Hess=0.5*(Hess+Hess.');
end

if isfinite(minLogProb)
    minProb=exp(max(minLogProb,log(realmin)));
else
    minProb=0;
end
v=struct('ell',ell,'g',g,'H',Hess,'valid',valid, ...
    'minProbability',minProb,'minLogProbability',minLogProb);
end

% =========================================================================
function [muHat,info]=multiStartEstimate(A,Y,args,o)
% Run the estimator from several initial points with the same projections
% and report how the solution depends on the start.
self=@estimateMuOptimizedProjections;
p=size(A,2); N=size(Y,1);
[~,d]=self(A,Y,args{:},'NumStarts',1,'DesignOnly',true);
C=d.C; pilot=d.muPilot(:); nnls=d.muNNLS(:);
base=max(pilot,o.PilotFloor);
rs=RandStream('mt19937ar','Seed',o.MultiStartSeed);
K=o.NumStarts;
S=zeros(p,K); labels=cell(1,K);
if isempty(o.Mu0), S(:,1)=nnls; labels{1}='NNLS'; else, S(:,1)=o.Mu0(:); labels{1}='Mu0'; end
S(:,2)=base; labels{2}='floored pilot';
for k=3:K
    S(:,k)=base.*exp(o.MultiStartSpread*randn(rs,p,1));
    labels{k}=sprintf('random %d',k-2);
end
M=nan(p,K); crit=nan(1,K); iters=nan(1,K); flags=nan(1,K); msgs=cell(1,K); secs=zeros(1,K);
parWorkers=0; if o.UseParallel, parWorkers=Inf; end
parfor (k=1:K, parWorkers)
    try
        [m,ik]=self(A,Y,args{:},'NumStarts',1,'ProjectionMatrix',C,'PilotMu',pilot, ...
            'Mu0',S(:,k),'ComputeSandwich',false,'UseParallel',false,'DesignOnly',false);
        M(:,k)=m(:); crit(k)=ik.criterion; iters(k)=ik.nIterations; flags(k)=ik.exitflag;
        secs(k)=ik.runtimeSeconds;
        msgs{k}=ik.message;
    catch ME
        msgs{k}=['failed: ' ME.message];
    end
end
valid=isfinite(crit);
if ~any(valid)
    error('estimateMuOptimizedProjections:MultiStartFailed', ...
        'All %d starts failed; first message: %s',K,msgs{1});
end
tmp=crit; tmp(~valid)=-Inf; [~,kb]=max(tmp);
% Final solve from the best solution: recomputes the diagnostics and the
% sandwich covariance at the reported estimate.
[muHat,info]=self(A,Y,args{:},'NumStarts',1,'ProjectionMatrix',C,'PilotMu',pilot, ...
    'Mu0',M(:,kb),'ComputeSandwich',o.ComputeSandwich,'DesignOnly',false);
nb=max(norm(muHat),eps);
rd=nan(1,K);
for k=find(valid), rd(k)=norm(M(:,k)-muHat(:))/nb; end
% Distinct local maxima: greedy clustering of the solutions by relative distance.
order=find(valid); [~,ix]=sort(crit(order),'descend'); order=order(ix);
reps=[];
for k=order
    if isempty(reps) || all(arrayfun(@(r)norm(M(:,k)-M(:,r))/nb,reps)>o.MultiStartTol)
        reps(end+1)=k; %#ok<AGROW>
    end
end
R=struct();
R.starts=S; R.labels=labels; R.solutions=M; R.criterion=crit;
R.relativeDistance=rd; R.iterations=iters; R.exitflag=flags; R.messages=msgs;
R.nStarts=K; R.nValid=nnz(valid);
R.fractionAtBest=nnz(rd(valid)<=o.MultiStartTol)/nnz(valid);
R.nDistinctMaxima=numel(reps);
R.maxRelativeDistance=max(rd(valid));
R.criterionGap=max(crit(valid))-min(crit(valid));
R.criterionGapPerObservation=R.criterionGap/max(N,1);
R.bestStart=kb;
info.multistart=R;
info.designRuntimeSeconds=d.designRuntimeSeconds;
% Fit time of all starts plus the final solve (the final solve alone starts at
% the best solution and is almost instantaneous). With parfor this is the sum
% of the workers' times, i.e. the computing effort, not the elapsed time.
info.finalSolveSeconds=info.runtimeSeconds;
info.runtimeSeconds=sum(secs)+info.runtimeSeconds;
info.nIterations=sum(iters(isfinite(iters)))+info.nIterations;
if isfield(d,'nSignatureCandidates'), info.nSignatureCandidates=d.nSignatureCandidates; end
if isfield(d,'nFallbackProjections'), info.nFallbackProjections=d.nFallbackProjections; end
end

% =========================================================================
function f=downweightedShare(mu,proj)
% Share of (observation, projection) pairs with weight (1-eps)p/p~ < 0.5.
num=0; den=0;
for h=1:numel(proj)
    cnt=double(proj(h).counts(:)); den=den+sum(cnt);
    if ~isfield(proj,'eps') || isempty(proj(h).eps) || proj(h).eps<=0, continue; end
    logP=compoundPoissonLogPmf(round(proj(h).b(:)),mu(:),proj(h).maxZ);
    lpz=logP(double(proj(h).zValues(:))+1);
    l1=log1p(-proj(h).eps)+lpz; l2=log(proj(h).eps)+proj(h).logG(:);
    mx=max(l1,l2); lmix=mx+log(exp(l1-mx)+exp(l2-mx));
    w=zeros(size(lpz)); fin=isfinite(lpz); w(fin)=exp(l1(fin)-lmix(fin));
    num=num+sum(cnt(w<0.5));
end
f=num/max(den,1);
end

% =========================================================================
function [L,g,H,ok]=contaminatedTerms(logP,lpz,b,zVals,a,eps,logG,pnum,needDer,logRealMax)
% Log-likelihood, gradient and Hessian of sum_r a_r log p~(z_r), with
% p~=(1-eps)p+eps*g. With w=(1-eps)p/p~ and r_j=p(z-b_j)/p(z):
%   d log p~ / d mu_j      = w (r_j-1)
%   d2 log p~ / d mu_j mu_k = w (r_jk - r_j - r_k + 1) - w^2 (r_j-1)(r_k-1).
% Observations with p(z)=0 have w=0: they enter only through eps*g.
ok=true; g=zeros(pnum,1); H=zeros(pnum);
l1=log1p(-eps)+lpz; l2=log(eps)+logG(:);
mx=max(l1,l2); lmix=mx+log(exp(l1-mx)+exp(l2-mx));
w=zeros(size(lpz)); fin=isfinite(lpz); w(fin)=exp(l1(fin)-lmix(fin));
L=struct('ell',sum(a.*lmix),'lmix',lmix,'w',w);
if ~needDer, return; end
ratio1=zeros(numel(zVals),pnum);
if any(fin)
    zf=zVals(fin); lf=lpz(fin);
    idx1=bsxfun(@minus,zf,b.'); ok1=idx1>=0; idxSafe=idx1+1; idxSafe(~ok1)=1;
    lr1=bsxfun(@minus,logP(idxSafe),lf);
    if any(lr1(ok1)>logRealMax), ok=false; return; end
    r=zeros(size(idx1)); r(ok1)=exp(lr1(ok1));
    ratio1(fin,:)=r;
end
score0=ratio1-1;
aw=a.*w;
g=(score0.'*aw);
maxShift=2*max(b);
shiftWeighted=zeros(maxShift+1,1);
for sh=0:maxShift
    idx=zVals-sh; use=idx>=0 & fin;
    if ~any(use), continue; end
    lr=logP(idx(use)+1)-lpz(use);
    if any(lr>logRealMax), ok=false; return; end
    shiftWeighted(sh+1)=sum(aw(use).*exp(lr));
end
firstTerm=shiftWeighted(bsxfun(@plus,b,b.')+1);
rsum=ratio1.'*aw;
H=firstTerm-rsum-rsum.'+sum(aw)-score0.'*bsxfun(@times,score0,a.*w.^2);
if any(~isfinite(H(:))) || any(~isfinite(g)), ok=false; end
end

% =========================================================================
function [c,nb]=sensitivityConditioning(mu,Hess)
% Conditioning of the composite criterion at the solution, restricted to
% the interior components (those at the bound are fixed by the constraint).
mu=mu(:); tol=1e-8*max(1,max(mu));
ii=find(mu>tol); nb=numel(mu)-numel(ii);
c=NaN;
if isempty(ii), return; end
if any(~isfinite(Hess(:))), c=Inf; return; end
S=-Hess(ii,ii); S=0.5*(S+S.');
% Scale invariance: intensities differ by orders of magnitude, which alone
% would make the raw condition number huge. Normalize to unit diagonal
% (a correlation-type matrix), so that c measures genuine near-dependence.
dg=diag(S);
if any(~(dg>0)), c=Inf; return; end
S=S./sqrt(dg*dg.');
e=eig(0.5*(S+S.'));
if min(e)<=0, c=Inf; else, c=max(e)/min(e); end
end

% =========================================================================
function [Sigma,se,interior]=compositeSandwich(mu,proj,weights,Z,Hess)
% Godambe (sandwich) covariance of the composite estimator,
%   Var(muHat) ~ Ihat^{-1} Vhat Ihat^{-1},  Ihat=-Hess,  Vhat=sum_r s_r s_r',
% where s_r is the composite score of observation r (all projections of the
% same observation are summed, so their dependence is accounted for).  The
% formula is valid only for interior components; components at the bound
% mu_j=0 are excluded and reported as NaN (nonstandard asymptotics).
mu=mu(:);
p=numel(mu);
N=size(Z,1);
tol=1e-8*max(1,max(mu));
interior=mu>tol;
S=zeros(N,p);
for h=1:numel(proj)
    b=round(proj(h).b(:));
    zVals=double(proj(h).zValues(:));
    logP=compoundPoissonLogPmf(b,mu,proj(h).maxZ);
    lpz=logP(zVals+1);
    fin=isfinite(lpz);
    idx=bsxfun(@minus,zVals,b.');
    ok=idx>=0 & repmat(fin,1,p);
    idx(~ok)=0;
    lr=logP(idx+1)-repmat(lpz,1,p);
    ratio=zeros(size(idx));
    ratio(ok)=exp(lr(ok));
    w=double(fin);
    if isfield(proj,'eps') && ~isempty(proj(h).eps) && proj(h).eps>0
        l1=log1p(-proj(h).eps)+lpz; l2=log(proj(h).eps)+proj(h).logG(:);
        mx=max(l1,l2); lmix=mx+log(exp(l1-mx)+exp(l2-mx));
        w=zeros(size(lpz)); w(fin)=exp(l1(fin)-lmix(fin));
    end
    scoreZ=bsxfun(@times,w,ratio-1);         % distinct z values by p
    [~,loc]=ismember(double(Z(:,h)),zVals);
    S=S+weights(h)*scoreZ(loc,:);
end
Sigma=nan(p);
se=nan(p,1);
ii=find(interior);
if isempty(ii) || any(~isfinite(Hess(:))), return; end
Ihat=-Hess(ii,ii);
Vhat=S(:,ii).'*S(:,ii);
if rcond(Ihat)<1e-12, return; end
Sig=Ihat\Vhat/Ihat;
Sig=0.5*(Sig+Sig.');
Sigma(ii,ii)=Sig;
se(ii)=sqrt(max(diag(Sig),0));
end

% =========================================================================
function r=projectedGradientResidual(mu,g)
if any(~isfinite(g))
    r=Inf;
else
    r=norm(mu-max(0,mu+g),Inf);
end
end

% =========================================================================
function [s,flag,nIter]=solveQuadraticModelFast(mu,g,H,Delta)
% Trust-region subproblem: maximize g's+0.5*s'Hs subject to ||s||<=Delta and
% lb<=s<=ub, lb=max(-mu,-Delta), ub=Delta. An active-set method fixes the
% variables at their lower bound and solves the ball-constrained problem on
% the free ones exactly (eigendecomposition + secular equation), adding
% violated bounds and releasing wrong-sign multipliers. It is started from
% three initial active sets; the best result, and the projected Cauchy step,
% are then polished by feasible projected-gradient steps that are accepted
% only if they improve the model. flag=1 (the outer ratio test decides).
p=numel(mu); mu=mu(:); g=g(:);
M=-0.5*(H+H.'); b=-g;                         % minimize f(s)=b's+0.5*s'Ms
lb=max(-mu,-Delta); ub=Delta*ones(p,1);
f=@(x) b.'*x+0.5*x.'*M*x;
% Projected Cauchy step
d=-b; d(mu<=1e-12 & d<0)=0; sc=zeros(p,1);
if norm(d)>0
    curv=d.'*M*d;
    if curv>0, t=min(Delta/norm(d),(d.'*d)/curv); else, t=Delta/norm(d); end
    sc=projectBoxBall(t*d,lb,ub,Delta);
end
inits={(mu<=1e-12)&(b>=0), false(p,1), sc<=lb+1e-12*max(1,Delta)};
best=sc; fbest=f(sc); nIter=0;
for k=1:numel(inits)
    [x,it]=activeSetBallBox(M,b,lb,ub,Delta,inits{k});
    nIter=nIter+it;
    if f(x)<fbest, best=x; fbest=f(x); end
end
% Feasible projected-gradient polish
alpha=1/max(norm(M,2),eps);
for k=1:50
    gr=b+M*best;
    x=projectBoxBall(best-alpha*gr,lb,ub,Delta);
    fx=f(x);
    if fx<fbest-1e-15*max(1,abs(fbest)), best=x; fbest=fx; else, alpha=alpha/2; end
    if alpha<1e-12/max(norm(M,2),eps), break; end
end
s=best; flag=1;
end

function [s,nIter]=activeSetBallBox(M,b,lb,ub,Delta,active)
p=numel(b); s=zeros(p,1); tol=1e-12*max(1,Delta); nIter=0;
for it=1:3*p+3
    nIter=it;
    s(active)=lb(active);
    F=~active;
    if ~any(F), break; end
    r2=Delta^2-sum(s(active).^2);
    if r2<=0, s(F)=0; break; end
    bF=b(F)+M(F,active)*s(active);
    s(F)=ballQuadMin(M(F,F),bF,sqrt(r2));
    viol=F & (s<lb-tol | s>ub+tol);
    if any(viol)
        s(viol)=min(max(s(viol),lb(viol)),ub(viol));
        active=active | viol;
        continue
    end
    gr=b+M*s;
    release=active & (s<=lb+tol) & (gr<-tol);
    if any(release), active(release)=false; continue; end
    break
end
s=projectBoxBall(s,lb,ub,Delta);
end

function x=projectBoxBall(x,lb,ub,Delta)
% Feasible point close to x: clip to the box, then scale into the ball
% (the box contains 0, so scaling keeps the bounds).
x=min(max(x,lb),ub);
n=norm(x); if n>Delta, x=x*(Delta/n); end
end

function x=ballQuadMin(M,b,r)
% Minimize b'x+0.5*x'Mx subject to ||x||<=r (exact, More-Sorensen type).
M=0.5*(M+M.');
[V,L]=eig(M); lam=diag(L); [lam,ord]=sort(lam); V=V(:,ord);
c=V.'*b;
if lam(1)>0
    x=-V*(c./lam);
    if norm(x)<=r, return; end
end
% Boundary solution: find sigma>=max(0,-lam1) with ||x(sigma)||=r
lo=max(0,-lam(1));
nrm=@(sig) norm(c./(lam+sig));
small=abs(c)<=1e-14*max(1,norm(c));
hard=small(1) && nrm(lo+1e-14*max(1,abs(lo)))<r;
if hard
    % Hard case: move along the leftmost eigenvector to reach the boundary
    sig=lo; w=zeros(size(c)); idx=lam>lam(1)+1e-12; w(idx)=-c(idx)./(lam(idx)+sig);
    x=V*w; extra=sqrt(max(r^2-norm(x)^2,0));
    x=x+extra*V(:,1);
    return
end
a=lo; hi=max(lo,1)+norm(b)/r+abs(lam(1));
while nrm(hi)>r, hi=2*hi+1; end
for k=1:100                                   % bisection on the monotone secular equation
    mid=0.5*(a+hi);
    if mid<=lo, mid=lo+eps(max(1,lo)); end
    if nrm(mid)>r, a=mid; else, hi=mid; end
    if hi-a<=1e-12*max(1,hi), break; end
end
x=-V*(c./(lam+hi));
end

function [s,flag,nCalls]=solveQuadraticModel(mu,g,H,Delta,options)
p=numel(mu);
H=0.5*(H+H.');
lb=-mu;
ub=Delta*ones(p,1);
lb=max(lb,-Delta);
% Projected gradient seed.
d=g;
d(mu<=1e-12 & d<0)=0;
if norm(d)>0
    s0=Delta*d/norm(d);
    s0=min(max(s0,lb),ub);
    if norm(s0)>Delta
        s0=s0*(Delta/norm(s0));
    end
else
    s0=zeros(p,1);
end
obj=@(s)quadObjective(s,g,H);
nonlcon=@(s)ballConstraint(s,Delta);
nCalls=0;
try
    [s,~,flag,out]=fmincon(obj,s0,[],[],[],[],lb,ub,nonlcon,options);
    if isstruct(out) && isfield(out,'funcCount'), nCalls=out.funcCount; end
catch
    s=zeros(p,1);
    flag=-3;
end
end

function [f,grad]=quadObjective(s,g,H)
f=-(g.'*s+0.5*s.'*H*s);
grad=-(g+H*s);
end

function [c,ceq,dc,dceq]=ballConstraint(s,Delta)
c=s.'*s-Delta^2;
ceq=[];
dc=2*s;
dceq=[];
end

% =========================================================================
function out=checkCompositeDerivatives(mu,proj,weights,pnum)
v=evaluateComposite(mu,proj,weights,pnum,true);
h=1e-5;
gFD=zeros(pnum,1);
HFD=zeros(pnum);
for j=1:pnum
    e=zeros(pnum,1); e(j)=h;
    vp=evaluateComposite(mu+e,proj,weights,pnum,true);
    mm=mu-e;
    if mm(j)<0
        vm=evaluateComposite(mu,proj,weights,pnum,true);
        gFD(j)=(vp.ell-vm.ell)/h;
        HFD(:,j)=(vp.g-vm.g)/h;
    else
        vm=evaluateComposite(mm,proj,weights,pnum,true);
        gFD(j)=(vp.ell-vm.ell)/(2*h);
        HFD(:,j)=(vp.g-vm.g)/(2*h);
    end
end
out=struct();
out.gradientRelativeError=norm(v.g-gFD)/max(1,norm(gFD));
out.hessianRelativeError=norm(v.H-HFD,'fro')/max(1,norm(HFD,'fro'));
out.gradientFD=gFD;
out.hessianFD=HFD;
end
