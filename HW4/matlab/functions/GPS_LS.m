function [Xu, H] = GPS_LS(Xs, X0, rho, rho_dot, W)
    Xu = X0;
    dX = 1e3 * ones(4,1);
    while norm(dX) > 1e-6
        dr = Xs(:,1:3) - Xu(1:3)';
        r = vecnorm(dr, 2, 2);
        rho_hat = r + Xu(4);
        U = dr./r;
        H = [-U ones(length(U),1)];
        dX = (H'*W*H)\H'*W*(rho - rho_hat);
        Xu(1:4) = Xu(1:4) + dX;
    end
    r_dot = rho_dot - dot(Xs(:,5:7), U, 2);
    Xu(5:8) = (H'*W*H)\H'*W*r_dot;
end