function [DOP, error] = gpsStats(sigma2, H, lla0)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    P = sigma2*inv(H'*H);
    C = eceftoenu(lla0(1), lla0(2));
    P(1:3, 1:3) = C*P(1:3, 1:3)*C';
    DOP.G = sqrt(sum(diag(P)./sigma2));
    DOP.P = sqrt(sum(diag(P(1:3,1:3))./sigma2));
    DOP.H = sqrt(sum(diag(P(1:2,1:2))./sigma2));
    DOP.V = sqrt(sum(diag(P(3,3))./sigma2));
    error.G = sqrt(sigma2)*DOP.G;
    error.P = sqrt(sigma2)*DOP.P;
    error.H = sqrt(sigma2)*DOP.H;
    error.V = sqrt(sigma2)*DOP.V;
    if size(H) > 3 
        DOP.T = sqrt(sum(diag(P(4,4))./sigma2));
        error.T = sqrt(sigma2)*DOP.T;
    end
end

