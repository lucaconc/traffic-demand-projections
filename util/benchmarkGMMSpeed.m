function T = benchmarkGMMSpeed(varargin)
%BENCHMARKGMMSPEED Solution time and error of the weighted projected scores.
%
%   T = benchmarkGMMSpeed()
%   T = benchmarkGMMSpeed('N',10000,'Replicates',3)
%
%   Times estimateMuEfficientGMM on the difficult benchmark matrix (all 15 nonzero
%   columns of {0,1}^4, rate scale 10) and its robust version on a 6-sensor
%   corridor (all 21 intervals of consecutive sensors, rate scale 2, 2% of the
%   intervals with a spike of 10 times the sensor mean), with the pilot options:
%     'pcl'       projected CL on its own small design (default)
%     'pclFull'   projected CL on all GMM projections (previous behaviour)
%     'cumulants' NNLS on the projections' cumulant equations only
%   and reports median seconds and relative error over the replicates, with the
%   projected CL and GLS for reference. Run it in MATLAB: Octave timings of the
%   projected CL are not representative.
p = inputParser;
addParameter(p,'N',10000); addParameter(p,'Replicates',3); addParameter(p,'Seed',1);
parse(p,varargin{:}); o = p.Results;
cases = benchmarkRoutingCases(); A1 = cases(3).A; mu1 = 10*cases(3).relRates;
m = 6; A2 = zeros(m,0);
for i = 1:m, for j = i:m, c = zeros(m,1); c(i:j) = 1; A2 = [A2 c]; end, end %#ok<AGROW>
n2 = size(A2,2); r2 = 0.5 + mod((1:n2)'*0.6180339887,1); mu2 = 2*r2/mean(r2);
rows = {};
for r = 1:o.Replicates
    rs = RandStream('mt19937ar','Seed',o.Seed+r);
    Y1 = poissinv(rand(rs,o.N,numel(mu1)),repmat(mu1.',o.N,1))*A1.';
    Y2 = poissinv(rand(rs,2000,n2),repmat(mu2.',2000,1))*A2.';
    Y2 = contaminateCounts(Y2,0.02,A2*mu2,10,rs);
    runs = {'Difficult: Projected CL',      @() estimateMuOptimizedProjections(A1,Y1,'Display','off'), mu1; ...
            'Difficult: Bootstrap GLS',     @() fitBootstrapGLS(A1,Y1,4,200,rs), mu1; ...
            'Difficult: GMM, pilot pcl',    @() estimateMuEfficientGMM(A1,Y1,'PlainPilot','pcl'), mu1; ...
            'Difficult: GMM, pilot pclFull',@() estimateMuEfficientGMM(A1,Y1,'PlainPilot','pclFull'), mu1; ...
            'Difficult: GMM, pilot cumulants',@() estimateMuEfficientGMM(A1,Y1,'PlainPilot','cumulants'), mu1; ...
            'Corridor 2%: Projected CL, robust',@() estimateMuOptimizedProjections(A2,Y2,'Display','off','Contamination',0.01), mu2; ...
            'Corridor 2%: robust GMM, full pilot',  @() estimateMuEfficientGMM(A2,Y2,'Robust',true), mu2; ...
            'Corridor 2%: robust GMM, small pilot', @() estimateMuEfficientGMM(A2,Y2,'Robust',true,'RobustPilotDesign','small'), mu2};
    for k = 1:size(runs,1)
        t0 = tic; est = runs{k,2}(); sec = toc(t0);
        rows(end+1,:) = {runs{k,1}, sec, norm(est(:)-runs{k,3})/norm(runs{k,3})}; %#ok<AGROW>
    end
end
R = cell2table(rows,'VariableNames',{'Setting','Seconds','RelError'});
[G,names] = findgroups(R.Setting);
T = table(names,splitapply(@median,R.Seconds,G),splitapply(@median,R.RelError,G), ...
    'VariableNames',{'Setting','MedianSeconds','MedianRelError'});
disp(T);
end
