%% GPS HOMEWORK 1
clear; close all; clc;

% Q1.3 Variables
dt = 0.1;                           % Time Step [s]
x_dot = 100;                        % Velocity [m/s]
c = physconst("LightSpeed");        % Speed of Light [m/s]
x = [1 1]';                         % Initial State [m]
f = 100e6;                          % Transmitter Frequency [Hz]
z = [-33.1679 -33.1711 -33.1743];   % Doppler Shifts [Hz]

% Q2 Variables
N = 100;
a1 = 2*ceil(rand(N,1)-0.5)-1;
a2 = 2*ceil(rand(N,1)-0.5)-1;
N_long = 1000;
a1_long = 2*ceil(rand(N_long,1)-0.5)-1;
a2_long = 2*ceil(rand(N_long,1)-0.5)-1;

% Q3 Variables
A = 3 + 3*randn(1000,1);
B = 5 + 5*randn(1000,1);
C = A + B;
DATA=[A B C];

%% QUESTION 1.3
lambda = c/f;                       % Wavelength [m]
y = (-z*lambda)';                   % Range Rates [m/s]
dx = 10;                            % State Delta [m]
% Least Squares
while norm(dx) > 1e-4
    p1 = x(1) + x_dot*dt;
    p2 = p1 + x_dot*dt;
    r0_hat = sqrt(x(1)^2 + x(2)^2);
    r1_hat = sqrt(p1^2 + x(2)^2);
    r2_hat = sqrt(p2^2 + x(2)^2);
    y_hat = [(x(1)*x_dot)/r0_hat, (p1*x_dot)/r1_hat, (p2*x_dot)/r2_hat]';
    dy = y - y_hat;
    H = [-(x(1)^2 + x(2)^2)^(-3/2)*(x(1)^2)*x_dot + x_dot/r0_hat, -(x(1)^2 + x(2)^2)^(-3/2)*x(2)*x(1)*x_dot;
         -(p1^2 + x(2)^2)^(-3/2)*(p1^2)*x_dot + x_dot/r1_hat, -(p1^2 + x(2)^2)^(-3/2)*x(2)*p1*x_dot;
         -(p2^2 + x(2)^2)^(-3/2)*(p2^2)*x_dot + x_dot/r2_hat, -(p2^2 + x(2)^2)^(-3/2)*x(2)*p2*x_dot;];
    dx = pinv(H)*dy;
    x = x + dx;
end

fprintf('1.3) Initial X & Y Position of the Aircraft: [%0.3g %0.3g] km\n\n', x./1000);

%% QUESTION 2A
figure();
tiledlayout(2,1);
nexttile();
histogram(a1);
xlabel("Value");
ylabel("Occurences");
title("Histogram of Random Sequence 1");
nexttile();
histogram(a2);
xlabel("Value");
ylabel("Occurences");
title("Histogram of Random Sequence 2");

%% QUESTION 2B
figure();
tiledlayout(2,1);
nexttile();
plot(abs(fft(a1)));
title("PSD of Random Sequence 1");
nexttile();
plot(abs(fft(a2)));
title("PSD of Random Sequence 2");

%% QUESTION 2C
acorr1 = supercorr(a1);
acorr2 = supercorr(a2);
figure();
tiledlayout(2,1);
nexttile();
plot(-N:N, acorr1);
title("Autocorrelation of Random Sequence 1");
nexttile();
plot(-N:N, acorr2);
title("Autocorrelation of Random Sequence 2");

%% QUESTION 2D
xcorr12 = supercorr(a1,a2);
xcorr21 = supercorr(a2,a1);
figure();
tiledlayout(2,1);
nexttile();
plot(-N:N, xcorr12);
title("Cross Correlation of Random Sequence 1 w/ 2");
nexttile();
plot(-N:N, xcorr21);
title("Cross Correlation of Random Sequence 2 w/ 1");

%% QUESTION 2 - BONUS
acorr1_long = supercorr(a1_long);
acorr2_long = supercorr(a2_long);
figure();
tiledlayout(2,1);
nexttile();
plot(-N_long:N_long, acorr1_long);
title("Autocorrelation of Random Sequence 1");
nexttile();
plot(-N_long:N_long, acorr2_long);
title("Autocorrelation of Random Sequence 2");

xcorr12_long = supercorr(a1_long,a2_long);
xcorr21_long = supercorr(a2_long,a1_long);
figure();
tiledlayout(2,1);
nexttile();
plot(-N_long:N_long, xcorr12_long);
title("Cross Correlation of Random Sequence 1 w/ 2");
nexttile();
plot(-N_long:N_long, xcorr21_long);
title("Cross Correlation of Random Sequence 2 w/ 1");

%% QUESTION 3A
meanA = mean(A);
meanB = mean(B);
meanC = mean(C);
fprintf("3a) Mean of A: %0.3g\n", meanA);
fprintf("3a) Mean of B: %0.3g\n", meanB);
fprintf("3a) Mean of C: %0.3g\n\n", meanC);

varA = var(A);
varB = var(B);
varC = var(C);
fprintf("3a) Variance of A: %0.3g\n", varA);
fprintf("3a) Variance of B: %0.3g\n", varB);
fprintf("3a) Variance of C: %0.3g\n\n", varC);

%% QUESTION 3B
meanDATA = mean(DATA, 'all');
fprintf("3b) Mean of Data: %0.3g\n\n", meanDATA);

%% QUESTION 3C
P_DATA = cov(DATA);
fprintf("3c) Covariance of Data:\n");
fprintf("\t[%0.3g\t %0.3g\t %0.3g]\n", P_DATA)

%% FUNCTIONS
function [Rxy] = supercorr(a,b)
    if ~exist('b','var')
        b = a;
    end
    N = length(a);
    Rxy = zeros(2*N,1);
    for k = 1:(2*N+1)
        shift = circshift(b,k - N - 1);
        Rxy(k) = sum(a.*shift)/N;
    end
end