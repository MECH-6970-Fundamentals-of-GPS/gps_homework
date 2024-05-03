function PER = modelPeriod(beta, phi_m)
PER = 0;
for n = 1:4
PER_i = beta(n)*(phi_m^(n-1));
PER = PER + PER_i;
end
if PER < 72000
    PER = 72000;
end
end