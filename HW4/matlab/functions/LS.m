function [s] = LS(y, H, s0)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    s = s0;
    ds = 1e3*ones(length(s0),1);
    while norm(ds) > 1e-4
        y_hat = H*s;
        ds = pinv(H)*(y - y_hat);
        s = s + ds;
    end
end

