function [mu,info] = fitCumulantNNLS(A,Y,order)
%FITCUMULANTNNLS Nonnegative cumulant matching through a given order.
%   Unweighted nonnegative least squares on the sample joint cumulants of all
%   sensor subsets of size <= order (order 1 = NNLS on the sensor means).
[C,subsets] = cumulantSystem(A,order);
k = sampleJointCumulants(Y,subsets);
mu = lsqnonneg(C,k,optimset('Display','off'));
info = struct('rank',rank(C),'equations',size(C,1));
end
