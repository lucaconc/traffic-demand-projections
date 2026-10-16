function k = sampleJointCumulants(Y,subsets)
%SAMPLEJOINTCUMULANTS Sample joint cumulants of distinct sensors.
%
%   k = sampleJointCumulants(Y,subsets)   Y: N x m counts; subsets: cell of index sets.
%
%   Order 1: mean. Order 2: unbiased covariance. Order 3: unbiased k-statistic
%   N/((N-1)(N-2)) sum c1 c2 c3. Order 4: plug-in m4 - (m12 m34 + m13 m24 + m14 m23)
%   from central moments (bias O(1/N)). Order 5: plug-in m5 minus the ten products
%   m_ab * m_cde over the splits of the five indices into a pair and a triple.
Y = double(Y); N = size(Y,1); c = Y - mean(Y,1);
k = zeros(numel(subsets),1);
for s = 1:numel(subsets)
    S = subsets{s};
    switch numel(S)
        case 1, k(s) = mean(Y(:,S));
        case 2, k(s) = sum(c(:,S(1)).*c(:,S(2)))/(N-1);
        case 3, k(s) = N/((N-1)*(N-2))*sum(c(:,S(1)).*c(:,S(2)).*c(:,S(3)));
        case 4
            m2 = @(a,b) mean(c(:,a).*c(:,b));
            k(s) = mean(prod(c(:,S),2)) - (m2(S(1),S(2))*m2(S(3),S(4)) + ...
                m2(S(1),S(3))*m2(S(2),S(4)) + m2(S(1),S(4))*m2(S(2),S(3)));
        case 5
            m2 = @(a,b) mean(c(:,a).*c(:,b));
            m3 = @(a,b,d) mean(c(:,a).*c(:,b).*c(:,d));
            pr = nchoosek(1:5,2); v = mean(prod(c(:,S),2));
            for t = 1:size(pr,1)
                rest = setdiff(1:5,pr(t,:));
                v = v - m2(S(pr(t,1)),S(pr(t,2)))*m3(S(rest(1)),S(rest(2)),S(rest(3)));
            end
            k(s) = v;
        otherwise
            error('sampleJointCumulants:Order','Orders above 5 are not implemented.');
    end
end
end
