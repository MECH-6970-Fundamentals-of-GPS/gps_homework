function AMP = ampVertDelay(alpha, phi_m)
AMP = 0;
for n = 1:4
    AMP_i = alpha(n)*(phi_m^(n-1));
    AMP = AMP + AMP_i;
    if AMP<0
        AMP = 0;
    end
end
end