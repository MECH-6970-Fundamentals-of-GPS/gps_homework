function [out] = upsample(vec,n)
%UPSAMPLE Summary of this function goes here
%   Detailed explanation goes here
    out_idx = ceil((1:n) .* ((length(vec))/(n)));
    out = vec(out_idx);
end

