function T = benchmarkEstimator(varargin)
%BENCHMARKESTIMATOR Time the estimator settings on a synthetic problem.
%
%   T = benchmarkEstimator()            % 40 sensors, 30 classes, 300 intervals
%   T = benchmarkEstimator('Sensors',60,'Classes',40,'Intervals',300)
%
%   Fits the same data with SubproblemSolver 'fmincon' and 'fast' and reports
%   the fit time, the number of trust-region iterations, the criterion reached
%   and the largest difference between the two estimates. Use it to see the
%   speed-up on your machine; both settings should reach the same criterion.
p = inputParser;
addParameter(p,'Sensors',40); addParameter(p,'Classes',30); addParameter(p,'Intervals',300);
addParameter(p,'Seed',1);
parse(p,varargin{:}); o = p.Results;
rs = RandStream('mt19937ar','Seed',o.Seed);
m = o.Sensors; n = o.Classes;
A = zeros(m,0);
while size(A,2) < n
    c = zeros(m,1); idx = randperm(rs,m); c(idx(1:1+randi(rs,4))) = 1;
    if ~any(all(A==c,1)), A = [A c]; end %#ok<AGROW>
end
mu = 2+18*rand(rs,n,1);
X = poissrnd(repmat(mu.',o.Intervals,1));
Y = X*A.';
solvers = {'fmincon','fast'}; rows = cell(2,1); est = cell(2,1);
for k = 1:2
    t = tic;
    [est{k},info] = estimateMuOptimizedProjections(A,Y,'Display','off','SubproblemSolver',solvers{k});
    rows{k} = table(string(solvers{k}),toc(t),info.runtimeSeconds,info.designRuntimeSeconds, ...
        info.nIterations,info.criterion,'VariableNames', ...
        {'SubproblemSolver','TotalSeconds','FitSeconds','DesignSeconds','Iterations','Criterion'});
end
T = vertcat(rows{:});
T.MaxRelDifference = [0; max(abs(est{2}-est{1}))/max(est{1})];
disp(T);
end
