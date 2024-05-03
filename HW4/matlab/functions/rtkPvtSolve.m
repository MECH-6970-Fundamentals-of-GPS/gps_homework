function [Xu, Xru, LLA] = rtkPvtSolve(Xs_r, X_r, dPsr, dCarr1, dCarr2, lambda1, lambda2, weights)
%RTKPVTSOLVE Summary of this function goes here
%   Detailed explanation goes here
    flag = false;
    if exist('weights', 'var')
        flag = true;
    end

    N = length(dPsr);
    for i = 1:N
        dPrns = ~isnan(dPsr(i,:)) & ~isnan(dCarr1(i,:)) & ~isnan(dCarr2(i,:));
        
        y = [dPsr(i,dPrns), dCarr1(i,dPrns), dCarr2(i,dPrns)]';
        Xs_i = squeeze(Xs_r(i,dPrns,1:3));
        if length(dPrns) < 6
            Xru(:,i) = Xru(:,i-1);
            Xu(:,i) = Xu(:,i-1);
        else
            r = vecnorm(Xs_i - X_r(1:3,i)',2,2);
            U = (Xs_i - X_r(1:3,i)')./r;
            H = [U, ones(length(U),1)];
            H = repmat(H,3,1);
            LN = null(H')';
            lambda = [zeros(length(dPsr(i,dPrns)), 2*length(dPsr(i,dPrns)));
                      lambda1*eye(length(dCarr1(i,dPrns))), zeros(length(dPsr(i,dPrns)), length(dPsr(i,dPrns)))
                      zeros(length(dPsr(i,dPrns)), length(dPsr(i,dPrns))), lambda2*eye(length(dCarr1(i,dPrns)))];
            LNlam = LN*lambda;
            LNy = LN*y;
            
            Nf = pinv(LNlam)*LNy;
            
            y2 = [dCarr1(i,dPrns), dCarr2(i,dPrns)]';
            H = [U, ones(length(U),1)];
            H = repmat(H,2,1);
            lambda = [lambda1*eye(length(dCarr1(i,dPrns))), zeros(length(dPsr(i,dPrns)), length(dPsr(i,dPrns)))
                      zeros(length(dPsr(i,dPrns)), length(dPsr(i,dPrns))), lambda2*eye(length(dCarr1(i,dPrns)))];
            if flag
                W = diag([weights(i,dPrns) weights(i+N,dPrns)]);
            else
                W = eye(2*sum(dPrns));
            end
            if sum(dPrns) >= 6
                Xru(:,i) = (H'*W*H)\H'*W*(y2-diag(lambda).*Nf);
                Xu(:,i) = X_r(1:3,i) + Xru(1:3,i);
            else
                Xru(:,i) = Xru(:,i-1);
                Xu(:,i) = Xu(:,i-1);
            end
        end
    end

    LLA = ecef2lla(Xu(1:3,:)');
end

