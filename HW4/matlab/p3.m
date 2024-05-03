clear; close all; clc;

currentFile = mfilename('fullpath');
currentFolder = fileparts(currentFile);
addpath(genpath(currentFolder + "/.."));

signal = generate_signal(2);
Ts = 1e-6;
Tint = 0.01;
t = 0:Ts:(Tint-Ts);
time = 0:Ts:Ts*(length(signal)-2);

[sumI, sumQ] = costas(signal, 1/Ts, Tint, 2);

figure();
hold on;
plot(sumI, 'LineWidth', 2);
plot(sumQ, '--', 'LineWidth', 2);
title('In-Phase & Quadrature');
xlabel('Sample');
legend('In-Phase', 'Quadrature');

exportgraphics(gcf, currentFolder + "/../figures/p3_IQ.png", 'Resolution', 300);

D = (1 + sign(sumI))./2;
D1 = D(1:1/Tint:end);
D2 = reshape(D1, 8, [])';
decode = bin2dec(num2str(D2));
fprintf('The message reads: %s\n', decode);

function [sumI, sumQ] = costas(signal, Fs, Tint, bandwidth, f0, phi0)
    if ~exist('Tint', 'var'); Tint = 10e-3; end
    if ~exist('bandwidth', 'var'); bandwidth = 1; end
    if ~exist('f0', 'var'); f0 = 100; end
    if ~exist('phi0', 'var'); phi0 = 0; end

    wn = 2*pi*bandwidth;
    zeta = 0.9;
    Ki = wn^2;
    Kp = 2*zeta*wn;
    
    Ts = 1/Fs;
    t = 0:Ts:(Tint-Ts);
    time = 0:Ts:Ts*(length(signal)-2);
    N = length(time);
    M = Tint/Ts;
    
    freq = f0.*ones(1,round(N/M));
    phase = phi0.*ones(1,round(N/M));
    err = zeros(1,round(N/M));
    ierr = zeros(1,round(N/M));
    sumI = zeros(1,round(N/M));
    sumQ = zeros(1,round(N/M));
    for k = 1:N/M-1
        X = signal(1  + (k-1)*M : k*M);
        replicaI = sin(2.*pi.*freq(k).*t + phase(k));
        replicaQ = cos(2.*pi.*freq(k).*t + phase(k));
        I = X.*replicaI;
        Q = X.*replicaQ;
        sumI(k) = sum(I);
        sumQ(k) = sum(Q);
        err(k+1) = atan(sumQ(k)/sumI(k));
        ierr(k+1) = ierr(k) + err(k+1)*Tint;
        freq(k+1) = freq(k) + (Kp*err(k+1) + Ki*ierr(k+1));
        phase(k+1) = rem(2.*pi*freq(k+1)*t(end) + phase(k), 2*pi);
    end
end