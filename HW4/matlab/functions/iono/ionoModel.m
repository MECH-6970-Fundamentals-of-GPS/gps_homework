function [psrL1_IM, weightsIM] = ionoModel(XsL1, truth_pos, psrL1, time, alpha, beta)
for i = 1:length(time)
    %DEFINE LAT AND LON OF USER AND SV AZ EL (in radians)
    phi_u = deg2rad(truth_pos(1))/pi;
    lambda_u = deg2rad(truth_pos(2))/pi;
    pos_sv = squeeze(XsL1(i,:,1:3));
    pos_sv = ecef2lla(pos_sv);
    for j = 1:32
        if isnan(pos_sv(j,1))
            A(j,1) = NaN; E(j,1) = NaN;
        else
            [A(j,1),E(j,1),~] = lookangles(truth_pos,pos_sv(j,:));
            A(j,1) = deg2rad(A(j,1))/pi;
            E(j,1) = deg2rad(E(j,1))/pi;
        end
    end

    % calculate amplitude of vertical delay 
    psi = angleUserAndEProjOfIonoIntPoint(E);
    phi_i = geodetLatEProjOfIonoIntPoint(phi_u, psi, A);
    lambda_i = geodetLonEProjOfIonoIntPoint(lambda_u, phi_i, psi, A);
    phi_m = geomagLatEProjOfIonoIntPoint(phi_i, lambda_i);
    AMP = ampVertDelay(alpha, phi_m);

    % calculate phase
    t = calcLocalTime(lambda_i, time(i));
    PER = modelPeriod(beta, phi_m);
    x = calcPhase(t, PER);

    %calculate obliquity factor
    F = obFact(E);
    weightsIM(i,:) = 1./F;

    %calculate corrections
    T_iono(i,:) = calcIonoCorr(F,AMP,x)';
    IM(i,:) = T_iono(i,:)*physconst('LightSpeed');
end

%correct and solve
psrL1_IM = psrL1 - IM;
end