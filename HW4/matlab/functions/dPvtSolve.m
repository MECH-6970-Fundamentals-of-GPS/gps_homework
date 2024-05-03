function [Xu, Xru, LLA] = dPvtSolve(Xs_r, Xr, psr_r, psr_u, dop_r, dop_u, weights)
%PVTSOLVE Summary of this function goes here
%   Detailed explanation goes here
    flag = false;
    if exist('weights', 'var')
        flag = true;
    end
    
    N = length(psr_r);
    Xru = zeros(8,N);
    Xu = Xru;
    for i = 1:N
        rho_ru = psr_r(i,:) - psr_u(i,:);
        rho_dot_ru = dop_r(i,:) - dop_u(i,:);
        prns = ~isnan(rho_ru) & ~isnan(rho_dot_ru);
        Xs_i = squeeze(Xs_r(i,prns,:));

        if flag
            W = diag(weights(i,prns));
        else
            W = eye(sum(prns));
        end

        if sum(prns) >= 4
            r = vecnorm(Xs_i(:,1:3) - Xr(1:3,i)', 2, 2);
            U = (Xs_i(:,1:3) - Xr(1:3,i)')./r;
            H = [-U ones(length(U),1)];
            Xru(1:4,i) = (H'*W*H)\H'*W*rho_ru(prns)';
            Xru(5:8,i) = (H'*W*H)\H'*W*rho_dot_ru(prns)';
            Xu(:,i) = Xr(:,i) - Xru(:,i);
        else
            Xru(:,i) = Xru(:,i-1);
            Xu(:,i) = Xu(:,i-1);
        end
    end
    LLA = ecef2lla(Xu(1:3,:)');
end