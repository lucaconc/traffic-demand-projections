function [Y,rows] = contaminateCounts(Y,fraction,meanCounts,spikeFactor,rs)
%CONTAMINATECOUNTS Detector glitches: spikes in a fraction of the intervals.
%
%   [Y,rows] = contaminateCounts(Y,fraction,meanCounts,spikeFactor,rs)
%
%   In round(fraction*N) randomly chosen intervals, one randomly chosen sensor i
%   gets spikeFactor*meanCounts(i)+1 extra counts (rounded). The other sensors
%   of the interval are unchanged. rs is a RandStream (or [] for the global one).
N = size(Y,1); m = size(Y,2); k = round(fraction*N);
rows = zeros(k,1);
if k == 0, return; end
if isempty(rs), pick = randperm(N,k); sens = randi(m,k,1);
else, pick = randperm(rs,N,k); sens = randi(rs,m,k,1); end
rows = pick(:);
for t = 1:k
    i = sens(t);
    Y(rows(t),i) = Y(rows(t),i) + round(spikeFactor*meanCounts(i) + 1);
end
end
