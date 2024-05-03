function [Xu, H, LLA] = pvtSolve(Xs, svPrns, psr, dop, weights)
%PVTSOLVE Summary of this function goes here
%   Detailed explanation goes here
    flag = false;
    if exist('weights', 'var')
        flag = true;
    end
    c = physconst('LightSpeed');
    fL1 = 1575.42e6;

    X0 = zeros(8,1);
    Xu = zeros(length(X0), length(psr));
    for i = 1:size(psr,1)
        Xs_i = squeeze(Xs(i,svPrns(i,:),:));
        rho = psr(i, svPrns(i,:)) + c*Xs_i(:,4)';
        rho_dot = dop(i, svPrns(i,:)) + c*Xs_i(:,8)';
        if flag
            W = diag(weights(i, svPrns(i,:)));
        else
            W = eye(length(psr(i, svPrns(i,:))));
        end

        if i == 1
            [Xu(:,i), H] = GPS_LS(Xs_i, Xu(:,i), rho', rho_dot', W);
        else
            [Xu(:,i), H] = GPS_LS(Xs_i, Xu(:,i-1), rho', rho_dot', W);
        end
    end
    LLA = ecef2lla(Xu(1:3,:)');
end

