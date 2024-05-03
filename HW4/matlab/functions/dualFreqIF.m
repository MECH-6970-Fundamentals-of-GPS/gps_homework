function [psrIF, psrVar] = dualFreqIF(psr1, psr2, f1, f2, var1, var2)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    psrIF = ((f1^2).*psr1 - (f2^2).*psr2)./(f1^2 - f2^2);

    if exist('var1', 'var') && exist('var2', 'var')
        psrVar = (f1^2/(f1^2 - f2^2)).*var1 + (f2^2/(f1^2 - f2^2)).*var2;
    else
        psrVar = NaN;
    end
end

