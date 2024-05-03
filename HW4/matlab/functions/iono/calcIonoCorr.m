function T_iono = calcIonoCorr(F,AMP,x)
if x<1.57
    T_iono = F*(5e-9 + AMP*(1 - (x^2)/4 + (x^4)/24));
else
    T_iono = F*(5e-9);
end
end