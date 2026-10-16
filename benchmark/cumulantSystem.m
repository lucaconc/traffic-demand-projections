function [C,subsets] = cumulantSystem(A,order)
%CUMULANTSYSTEM Linear system of Poisson joint cumulants through a given order.
%
%   [C,subsets] = cumulantSystem(A,r)
%
%   For independent Poisson path counts X and Y=A*X with binary A, the joint
%   cumulant of the sensors in a set S (of any order, with or without repeated
%   indices) is sum_j mu_j * prod_{i in S} a_ij. Repeated indices give the same
%   row, so the distinct equations through order r are indexed by the nonempty
%   sensor subsets of size at most r: C(s,j)=1 if path j crosses all sensors of
%   subsets{s}. Rows that are identically zero are dropped.
m = size(A,1); A = double(A~=0);
subsets = {}; C = zeros(0,size(A,2));
for s = 1:min(order,m)
    combos = nchoosek(1:m,s);
    for k = 1:size(combos,1)
        row = all(A(combos(k,:),:),1);
        if any(row)
            subsets{end+1,1} = combos(k,:); %#ok<AGROW>
            C(end+1,:) = row; %#ok<AGROW>
        end
    end
end
end
