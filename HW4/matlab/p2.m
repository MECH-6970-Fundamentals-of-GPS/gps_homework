clear; close all; clc;

currentFile = mfilename('fullpath');
currentFolder = fileparts(currentFile);
addpath(genpath(currentFolder + "/.."));

signal = generate_signal(1);
Ts = 1e-6;
Tint = 0.01;
time = 0:Ts:Ts*(length(signal)-2);
sampledTime = 0:Tint/Ts:length(time)-1;
intTime = 0:Ts:Tint;

% (A)
[freq1, phase1, err1] = costas(signal, 1/Ts, Tint);

figure();
hold('on');
plot(sampledTime, rad2deg(phase1), 'LineWidth', 2);
plot(sampledTime, rad2deg(err1), '--', 'LineWidth', 2);
title('Phase & Error vs. Time');
subtitle('Bandwidth: 1Hz');
xlabel('Time (s)');
ylabel('Phase/Error (rad)')
legend('Phase', 'Phase Error');
ax = gca;
ax.FontSize = 18;

exportgraphics(gcf, currentFolder + "/../figures/p2a_phase.png", 'Resolution', 300);

figure();
hold('on');
plot(length(time)*Ts - Tint*Ts + intTime, signal(length(signal)-Tint/Ts:end), 'LineWidth', 2);
plot(length(time)*Ts - Tint*Ts + intTime, sin(2.*pi.*freq1(end).*intTime + phase1(end)), '--', 'LineWidth', 2);
title('Signal & Replica vs. Time');
subtitle('Bandwidth: 1Hz');
xlabel('Time (s)');
ylabel('Signal (rad)');
legend('Signal', 'Replica');
ax = gca;
ax.FontSize = 18;

exportgraphics(gcf, currentFolder + "/../figures/p2a_signal.png", 'Resolution', 300);

% (B)
[freq2, phase2, err2] = costas(signal, 1/Ts, Tint, 2);

figure();
hold('on');
plot(sampledTime, rad2deg(phase2), 'LineWidth', 2);
plot(sampledTime, rad2deg(err2), '--', 'LineWidth', 2);
title('Phase & Error vs. Time');
subtitle('Bandwidth: 2Hz');
xlabel('Time (s)');
ylabel('Phase/Error (rad)')
legend('Phase', 'Phase Error');
ax = gca;
ax.FontSize = 18;

exportgraphics(gcf, currentFolder + "/../figures/p2b_phase.png", 'Resolution', 300);

figure();
hold('on');
plot(length(time)*Ts - Tint*Ts + intTime, signal(length(signal)-Tint/Ts:end), 'LineWidth', 2);
plot(length(time)*Ts - Tint*Ts + intTime, sin(2.*pi.*freq2(end).*intTime + phase2(end)), '--', 'LineWidth', 2);
title('Signal & Replica vs. Time');
subtitle('Bandwidth: 2Hz');
xlabel('Time (s)');
ylabel('Signal (rad)');
legend('Signal', 'Replica');
ax = gca;
ax.FontSize = 18;

exportgraphics(gcf, currentFolder + "/../figures/p2b_signal.png", 'Resolution', 300);

% (D)
[freq10, phase10, err10] = costas(signal, 1/Ts, Tint, 2, 10);

figure();
hold('on');
plot(sampledTime, rad2deg(phase10), 'LineWidth', 2);
plot(sampledTime, rad2deg(err10), '--', 'LineWidth', 2);
title('Phase & Error vs. Time');
subtitle('Bandwidth: 2Hz');
xlabel('Time (s)');
ylabel('Phase/Error (rad)')
legend('Phase', 'Phase Error');
ax = gca;
ax.FontSize = 18;

exportgraphics(gcf, currentFolder + "/../figures/p2d_phase.png", 'Resolution', 300);

figure();
hold('on');
plot(length(time)*Ts - Tint*Ts + intTime, signal(length(signal)-Tint/Ts:end), 'LineWidth', 2);
plot(length(time)*Ts - Tint*Ts + intTime, sin(2.*pi.*freq10(end).*intTime + phase10(end)), '--', 'LineWidth', 2);
title('Signal & Replica vs. Time');
subtitle('Initial Frequency 10Hz');
xlabel('Time (s)');
ylabel('Signal (rad)');
legend('Signal', 'Replica');
ax = gca;
ax.FontSize = 18;

exportgraphics(gcf, currentFolder + "/../figures/p2d_signal.png", 'Resolution', 300);

%% FUNCTIONS
function [freq, phase, err] = costas(signal, Fs, Tint, bandwidth, f0, phi0)
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
    
    freq = f0.*ones(1,N/M);
    phase = phi0.*ones(1,N/M);
    err = zeros(1,N/M);
    ierr = zeros(1,N/M);
    for k = 1:N/M-1
        X = signal(1  + (k-1)*M : k*M);
        replicaI = sin(2.*pi.*freq(k).*t + phase(k));
        replicaQ = cos(2.*pi.*freq(k).*t + phase(k));
        I = X.*replicaI;
        Q = X.*replicaQ;
        err(k+1) = atan2(sum(Q), sum(I));
        ierr(k+1) = ierr(k) + err(k+1)*Tint;
        freq(k+1) = freq(k) + (Kp*err(k+1) + Ki*ierr(k+1));
        phase(k+1) = rem(2.*pi*freq(k+1)*t(end) + phase(k), 2*pi);
    end
end