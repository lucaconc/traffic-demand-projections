function cases = benchmarkRoutingCases()
%BENCHMARKROUTINGCASES Routing matrices with controlled identification difficulty.
%
%   cases(k).name, .label, .A (sensors x paths), .relRates (mean 1),
%   .identifyingOrder (smallest cumulant order whose system has full rank).
%
%   easy          the running example: 2 sensors, 3 paths; covariance identifies
%   intermediate  all 7 nonzero binary columns of {0,1}^3; needs third order
%   difficult     all 15 nonzero binary columns of {0,1}^4; needs fourth order
%                 (cumulant ranks through orders 1..4: 4, 10, 14, 15; Lev-Ari et al.)
%   complex       all 31 nonzero binary columns of {0,1}^5; needs fifth order
%                 (cumulant ranks through orders 1..5: 5, 15, 25, 30, 31)
%   Relative rates are fixed, strictly positive and distinct (an irregular
%   sequence between 0.5 and 1.5, normalized to mean one), and are scaled by
%   the benchmark's overall rate scale.
names  = {'easy','intermediate','difficult','complex'};
labels = {'Easy (2 sensors, 3 paths)','Intermediate (3 sensors, 7 paths)', ...
          'Difficult (4 sensors, 15 paths)','Complex (5 sensors, 31 paths)'};
A2 = [1 0 1; 0 1 1];
cases = struct('name',{},'label',{},'A',{},'relRates',{},'identifyingOrder',{});
for k = 1:4
    if k == 1, A = A2; else, A = allBinaryColumns(k+1); end
    n = size(A,2);
    % Irregular fixed sequence (fractional parts of multiples of the golden
    % ratio): strictly positive, distinct, and without the special structure of
    % an arithmetic progression, which can place mu in the row space of A and
    % make the minimum-norm solution of the mean equations exact by accident.
    r = 0.5 + mod((1:n)*0.6180339887,1);
    r = r/mean(r);
    ord = find(arrayfun(@(o) rank(cumulantSystem(A,o)) == n, 1:size(A,1)),1);
    cases(k) = struct('name',names{k},'label',labels{k},'A',A,'relRates',r(:), ...
        'identifyingOrder',ord);
end
end

function A = allBinaryColumns(m)
% All 2^m-1 nonzero binary columns of length m, ordered by number of ones.
c = dec2bin(1:2^m-1,m)-'0';
[~,ord] = sortrows([sum(c,2) -c]);
A = c(ord,:).';
end
