function [Z,L] = zTransform(L)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    n = length(L);
    Z = eye(n);
    for i = n:-1:1
        for j = i + 1:n
            mu = round(L(j,i));
            if mu ~= 0
                L(j:n,i) = L(j:n,i) - mu*L(j:n,j);
                Z(1:n,i) = Z(1:n,i) - mu*Z(1:n,j);
            end
        end
    end
end

