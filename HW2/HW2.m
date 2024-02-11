%% QUESTION 1
clear; close all; clc;
a = 3;                          % True State
N_max = 50;                     % Number of Samples
M = 1e4;                        % Number of Monte Carlos
y_var = 1;                      % Measurement Variance

P = zeros(N_max,1);             % True Covariance
P_hat = zeros(N_max,1);         % Monte Carlo Covariance
for N = 1:N_max
    H = ones(N,1);              % Geometry Matrix
    P(N) = y_var/(H'*H);        % True Covariance

    a_hat = zeros(M,1);         % State Estimates
    for i = 1:M
        y = H*a + sqrt(y_var)*randn(N,1);       % Noisy Measurements
        a_hat(i) = LS(y, H, a_hat(i));          % Least Squares Estimate
    end

    P_hat(N) = var(a_hat);
end

figure();
hold("on");
title("Covariance v. Number of Samples");
plot(P, '-o');
plot(P_hat, '-*');
xlabel("Number of Samples");
ylabel("Variance");
legend("True Covariance", "Estimated Covariance");

%% QUESTION 2
clear;
x = 0:4;
y = [0.181 2.680 3.467 3.101 3.437]';
y_var = 0.4^2;

Hab = [ones(length(x),1) x'];
ab = LS(y, Hab, [0 0]');
Pab = y_var*inv(Hab'*Hab);
fprintf('Linear Fit Coefficients: [%0.4g %0.4g]\n', ab);

Habc = [Hab (x.^2)'];
abc = LS(y, Habc, [0 0 0]');
Pabc = y_var*inv(Habc'*Habc);
fprintf('Quadratic Fit Coefficients: [%0.4g %0.4g %0.4g]\n', abc);

Habcd = [Habc (x.^3)'];
abcd = LS(y, Habcd, [0 0 0 0]');
Pabcd = y_var*inv(Habcd'*Habcd);
fprintf('Cubic Fit Coefficients: [%0.4g %0.4g %0.4g %0.4g]\n', abcd);

fprintf('Estimation error for a:\n');
fprintf('\tLinear: %0.4g\n', sqrt(Pab(1,1)));
fprintf('\tQuadratic: %0.4g\n', sqrt(Pabc(1,1)));
fprintf('\tCubic: %0.4g\n\n', sqrt(Pabcd(1,1)));

%% QUESTION 3
clear;
a = [0 10 0 10]';
b = [0 0 10 10]';
r2 = [25 65 45 85]';
r_var = 0.5^2;

H_func = @(s) [2*(s(1)-a) 2*(s(2)-b)];
s = [0 0]';
ds = 1e3;
while ds > 1e-4
    r2_hat = (s(1) - a).^2 + (s(2) - b).^2;
    H = H_func(s);
    ds = pinv(H)*(r2 - r2_hat);
    s = s + ds;
end

P = r_var*inv(H'*H);
err = sqrt(diag(P));

M = 1e3;
s_hat = zeros(2,M);
for i = 1:M
    r2 = [25 65 45 85]' + sqrt(r_var)*randn(4,1);
    s = [0 0]';
    ds = 1e3;
    while ds > 1e-4
        r2_hat = (s(1) - a).^2 + (s(2) - b).^2;
        H = H_func(s);
        ds = pinv(H)*(r2 - r2_hat);
        s = s + ds;
    end
    s_hat(:,i) = s;
end
err_hat = std(s_hat')';

fprintf('The actual error on the estimate: [%0.3g %0.3g]\n', err);
fprintf('The estimated error on the estimate: [%0.3g %0.3g]\n\n', err_hat);

%% QUESTION 4
clear;
filename = 'HW2_data.txt';
T = readlines(filename);
Xs = zeros(length(T)-4, 3);
rho = zeros(length(T)-4, 1);
r_var = 0.5^2;
for i = 4:length(T)
    data = strsplit(T(i));
    Xs(i-3,1) = str2double(data(2));
    Xs(i-3,2) = str2double(data(3));
    Xs(i-3,3) = str2double(data(4));
    rho(i-3) = str2double(data(5));
end

X0 = [0 0 0 0];
[Xu_4, H_4] = GPS_LS(rho(1:4), Xs(1:4,:), X0);
[Xu_9, H_9] = GPS_LS(rho(1:9), Xs(1:9,:), X0);
[Xu_cl, H_cl] = GPS_LS(rho(1:4), Xs(1:4,:), [X0(1:3),Xu_9(4)]);
% [Xu_SOOP, H_SOOP] = GPS_LS(rho(1:9), Xs(1:9,:), X0);      % Infinite Loop
% Poor Geometry
X_guess = [423000 -5362000 3417000 0];
[Xu_SOOP, H_SOOP] = GPS_LS([rho(1:2);rho(10:11)], [Xs(1:2,:);Xs(10:11,:)], X_guess);

lla0 = ecef2lla(Xu_9(1:3));
C = ECEF_ENU(lla0(1), lla0(2))';

P_4 = r_var*inv(H_4'*H_4);
P_4(1:3,1:3) = C*P_4(1:3,1:3)*C';
H_error4 = norm(diag(P_4(1:2,1:2)));
V_error4 = norm(P_4(3,3));
G_error4 = norm(diag(P_4));
P_9 = r_var*inv(H_9'*H_9);
P_9(1:3,1:3) = C*P_9(1:3,1:3)*C';
H_error9 = norm(diag(P_9(1:2,1:2)));
V_error9 = norm(P_9(3,3));
G_error9 = norm(diag(P_9));
P_cl = r_var*inv(H_cl'*H_cl);
P_cl(1:3,1:3) = C*P_cl(1:3,1:3)*C';
H_errorCl = norm(diag(P_cl(1:2,1:2)));
V_errorCl = norm(P_cl(3,3));
G_errorCl = norm(diag(P_cl));
P_SOOP = r_var*inv(H_SOOP'*H_SOOP);
P_SOOP(1:3,1:3) = C*P_SOOP(1:3,1:3)*C';
H_errorSOOP = norm(diag(P_SOOP(1:2,1:2)));
V_errorSOOP = norm(P_SOOP(3,3));
G_errorSOOP = norm(diag(P_SOOP));

x = ["4 SVs" "9 SVs" "4 SVs + Perfect Clock" "2 SVs + 2 SOOPs"];
y = [G_error4 G_error9 G_errorCl G_errorSOOP;
     H_error4 H_error9 H_errorCl H_errorSOOP;
     V_error4 V_error9 V_errorCl V_errorSOOP];

figure();
hold("on");
title('GDOP for Varying Scenarios')
bar(x,y, "stacked");
xlabel('Scenario');
ylabel('GDOP');
legend('PDOP', 'HDOP', 'GDOP')

%% QUESTION 5
c = physconst("lightspeed");
A1 = 5e-9;
I_func = @(theta) A1 * (1 + 16*(0.53 - theta/180)^3);

rho9 = rho(1:9);
Xs9 = Xs(1:9,:);

[E, N, U] = ecef2enu(Xs9(:,1), Xs9(:,2), Xs9(:,3), lla0(1), lla0(2), lla0(3), wgs84Ellipsoid('kilometer'));
h = vecnorm([E, N]');
ang = atan2d(U, h');

i = 1;
for a = 1:5:50
    Xs9_filt = Xs9(ang > a, :);
    rho_filt = rho(ang > a);
    if length(rho_filt) > 4
        dx = 1e3*ones(4,1);
        state = [0 0 0 0];
        while norm(dx) > 1e-4
            r = vecnorm((Xs9_filt - state(1:3)), 2, 2);
            U = (Xs9_filt - state(1:3))./r;
            H = [-U ones(length(U),1)];
            rho_hat = r + state(4) + c*I_func(i);
            dx = pinv(H)*(rho_filt - rho_hat);
            state = state + dx';
        end
        bleh(i,:) = state;
        P = r_var*inv(H'*H);
        P(1:3,1:3) = P(1:3,1:3)*ECEF_ENU(lla0(1), lla0(2));
        a_range(i) = a;
        GDOP(i) = vecnorm(diag(P));
        PDOP(i) = vecnorm(diag(P(1:3,1:3)));
        HDOP(i) = vecnorm(diag(P(1:2,1:2)));
        VDOP(i) = P(3,3);
        TDOP(i) = P(4,4);
        i = i + 1;
    else
        break;
    end
end

X_err = abs(bleh - Xu_9);

figure();
hold("on");
title("GDOP vs. Mask Angle");
plot(a_range, GDOP, '*');
plot(a_range, PDOP, '*');
plot(a_range, HDOP, '*');
plot(a_range, VDOP, '*');
plot(a_range, TDOP, '*');
xlabel('Mask Angle');
ylabel('DOP');
legend('GDOP', 'PDOP', 'HDOP', 'VDOP', 'TDOP');

%% FUNCTIONS
function [s] = LS(y, H, s0)
s = s0;
ds = 1e3*ones(length(s0),1);
while norm(ds) > 1e-4
    y_hat = H*s;
    ds = pinv(H)*(y - y_hat);
    s = s + ds;
end
end

function [Xu, H] = GPS_LS(rho, Xs, Xu)
    dx = 1e3*ones(4,1);
    flag = false;
    if Xu(4) ~= 0
        b = Xu(4);
        Xu = [0 0 0];
        flag = true;
    end
    while norm(dx) > 1e-4
        r = vecnorm((Xs - Xu(1:3)), 2, 2);
        U = (Xs - Xu(1:3))./r;
        if flag
            H = [-U];
            rho_hat = r + b;
        else
            H = [-U ones(length(U),1)];
            rho_hat = r + Xu(4);
        end
        dx = pinv(H)*(rho - rho_hat);
        Xu = Xu + dx';
    end
end

function [C] = ECEF_ENU(lat, lon)
    C = [-sind(lon), -cosd(lon)*sind(lat), cosd(lon)*cosd(lat);
          cosd(lon), -sind(lon)*sind(lat), sind(lon)*cosd(lat);
                  0,            cosd(lat),           sind(lat)];
end