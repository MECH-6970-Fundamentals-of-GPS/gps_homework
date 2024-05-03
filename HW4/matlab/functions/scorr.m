function [Rxy] = scorr(a, b)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    if ~exist('b','var')
        b = a;
    end
    N = length(a);
    Rxy = zeros(2*N,1);
    for k = 1:(2*N+1)
        shift = circshift(b,k - N - 1);
        Rxy(k) = sum(a.*shift)/N;
    end
end

