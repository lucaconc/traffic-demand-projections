function tests = testEstimator
%TESTESTIMATOR Unit tests of estimateMuOptimizedProjections on the running
%example of the paper (two sensors, three paths, A=[1 0 1;0 1 1]).
%   Run with runPaperLibTests or runtests('testEstimator').
tests = functiontests(localfunctions);
end

function setupOnce(tc)
addpath(fileparts(fileparts(mfilename('fullpath'))));
setupPaperLib();
tc.TestData.A = [1 0 1;0 1 1];
tc.TestData.mu = [4;3;5];
end

function Y = sample(tc,N,seed)
rng(seed,'twister');
A = tc.TestData.A; mu = tc.TestData.mu;
X = zeros(N,3);
for j = 1:3, X(:,j) = poissonInverse(mu(j),N); end
Y = X*A.';
end

function testRecoversIntensities(tc)
% Asymptotic standard deviation is about 0.415*sqrt(2000/N) per component.
Y = sample(tc,4000,11);
muHat = estimateMuOptimizedProjections(tc.TestData.A,Y,'Display','off');
verifyLessThan(tc,max(abs(muHat(:)-tc.TestData.mu)),1.25);
end

function testRankConditionHolds(tc)
Y = sample(tc,1000,12);
[~,info] = estimateMuOptimizedProjections(tc.TestData.A,Y,'Display','off');
verifyTrue(tc,info.rankConditionHolds);
verifyGreaterThanOrEqual(tc,info.nSignatureCandidates,1);
end

function testSignatureFallback(tc)
% One requested projection cannot identify three intensities; the fallback
% must append a signature projection and restore the rank condition.
Y = sample(tc,1000,13);
[~,info] = estimateMuOptimizedProjections(tc.TestData.A,Y,'Display','off', ...
    'MaxCombinationSize',1,'NumProjections',1);
verifyTrue(tc,info.rankConditionHolds);
verifyGreaterThanOrEqual(tc,info.nFallbackProjections,1);
end

function testRankFailsWithoutSignatures(tc)
Y = sample(tc,1000,14);
[~,info] = verifyWarning(tc,@()estimateMuOptimizedProjections(tc.TestData.A,Y, ...
    'Display','off','MaxCombinationSize',1,'IncludeSignatureProjections',false, ...
    'SignatureRankFallback',false),'estimateMuOptimizedProjections:ProjectionRank');
verifyFalse(tc,info.rankConditionHolds);
end

function testSandwichStandardErrors(tc)
% Theory (Proposition 4.3) gives about 0.415 at N=2000 for each component.
Y = sample(tc,2000,15);
[~,info] = estimateMuOptimizedProjections(tc.TestData.A,Y,'Display','off');
verifyTrue(tc,all(info.interiorMask));
verifyGreaterThan(tc,info.standardErrors,0.30*ones(3,1));
verifyLessThan(tc,info.standardErrors,0.55*ones(3,1));
end

function testDuplicateColumnsRejected(tc)
verifyError(tc,@()estimateMuOptimizedProjections([1 1 0;0 0 1],ones(10,2), ...
    'Display','off'),'estimateMuOptimizedProjections:DuplicateColumns');
end

function x = poissonInverse(lam,n)
% Inversion sampler, independent of the Statistics Toolbox.
x = zeros(n,1); u = rand(n,1); p0 = exp(-lam);
for i = 1:n
    k = 0; p = p0; F = p0;
    while u(i) > F, k = k+1; p = p*lam/k; F = F+p; end
    x(i) = k;
end
end

function testContaminationHandlesImpossibleCounts(tc)
% Two sensors on the same path plus a third one; extra counts on sensor 1
% make Y1+Y2 odd, which the exact model cannot produce. eps=0 must fail,
% eps>0 must recover the intensities.
rng(21,'twister');
A = [1 0;1 0;0 1]; mu = [6;4]; N = 1500;
X = [poissonInverse(mu(1),N) poissonInverse(mu(2),N)];
Y = X*A.'; Y(:,1) = Y(:,1)+(rand(N,1)<0.15);
C = [1 1 0;0 0 1];
opts = {'Display','off','ProjectionMatrix',C,'PilotMu',[1;1],'Mu0',[1;1]};
verifyError(tc,@()estimateMuOptimizedProjections(A,Y,opts{:}), ...
    'estimateMuOptimizedProjections:InitialPoint');
[muHat,info] = estimateMuOptimizedProjections(A,Y,opts{:},'Contamination',0.01,'CheckDerivatives',true);
verifyLessThan(tc,abs(muHat(:)-mu),[0.4;0.4]);
verifyLessThan(tc,info.derivativeCheck.gradientRelativeError,1e-5);
verifyLessThan(tc,info.derivativeCheck.hessianRelativeError,1e-5);
verifyGreaterThan(tc,info.downweightedFraction,0);
end

function testMultiStartReport(tc)
% On the running example all starts should reach the same maximum.
Y = sample(tc,1000,31);
mu1 = estimateMuOptimizedProjections(tc.TestData.A,Y,'Display','off');
[mu4,info] = estimateMuOptimizedProjections(tc.TestData.A,Y,'Display','off','NumStarts',4);
R = info.multistart;
verifyEqual(tc,R.nStarts,4);
verifyEqual(tc,R.nValid,4);
verifyEqual(tc,R.nDistinctMaxima,1);
verifyEqual(tc,R.fractionAtBest,1);
verifyLessThan(tc,R.maxRelativeDistance,1e-2);
verifyLessThan(tc,max(abs(mu4(:)-mu1(:))),1e-3*max(1,max(mu1)));
verifyTrue(tc,all(isfinite(info.standardErrors)));
end

function testDesignOnlyIgnoresNumStarts(tc)
% The Parma pipeline selects shared projections with DesignOnly=true while
% passing NumStarts>1; the design call must not attempt a multi-start fit.
Y = sample(tc,600,41);
opts = {'Display','off','NumStarts',5,'Contamination',0.01};
[~,d] = estimateMuOptimizedProjections(tc.TestData.A,Y,opts{:},'DesignOnly',true,'PilotMu',[1;1;1]);
verifyTrue(tc,d.designOnly);
[~,info] = estimateMuOptimizedProjections(tc.TestData.A,Y,opts{:},'ProjectionMatrix',d.C, ...
    'PilotMu',[1;1;1],'Mu0',[1;1;1]);
verifyEqual(tc,info.multistart.nValid,5);
end

function testEfficientGMMRecoversIntensities(tc)
% Efficiently weighted projected scores: consistent on the running example.
rs = RandStream('mt19937ar','Seed',41);
A = [1 0 1; 0 1 1]; mu = [5; 3.5; 6.5]; N = 4000;
X = zeros(N,3); for j = 1:3, X(:,j) = poissinv(rand(rs,N,1),mu(j)); end
Y = X*A.';
[muHat,info] = estimateMuEfficientGMM(A,Y);
verifyLessThan(tc,norm(muHat-mu)/norm(mu),0.08);
verifyEqual(tc,info.exitflag,1);
verifyGreaterThan(tc,info.nMoments,size(info.C,1));
end

function testRobustGMMResistsSpikes(tc)
% Robust GMM: close to the truth with 5% of the intervals carrying a large spike,
% where the plain GMM breaks down.
rs = RandStream('mt19937ar','Seed',43);
A = [1 0 1; 0 1 1]; mu = [5; 3.5; 6.5]; N = 3000;
X = zeros(N,3); for j = 1:3, X(:,j) = poissinv(rand(rs,N,1),mu(j)); end
Y = X*A.';
bad = randperm(rs,N,round(0.05*N)); Y(bad,1) = Y(bad,1) + 120;      % detector glitches
muR = estimateMuEfficientGMM(A,Y,'Robust',true);
muP = estimateMuEfficientGMM(A,Y);
verifyLessThan(tc,norm(muR-mu)/norm(mu),0.10);
verifyGreaterThan(tc,norm(muP-mu)/norm(mu),norm(muR-mu)/norm(mu));
end
