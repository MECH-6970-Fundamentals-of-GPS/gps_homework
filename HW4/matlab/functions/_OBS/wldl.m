function [L,D] = wldl(Q)
%WLDL Summary of this function goes here
%   Detailed explanation goes here
    
    n = length(Q);
    L = zeros(size(Q));
    D = diag(diag(Q));
    for i = n:-1:1
        L(i,1:i) = Q(i,1:i)/sqrt(Q(i,i));
        for j = 1:i-1
            Q(j,1:j) = Q(j,1:j) - L(i,1:j)*L(i,j);
        end
        L(i,1:i) = L(i,1:i)/L(i,i);
    end
end

