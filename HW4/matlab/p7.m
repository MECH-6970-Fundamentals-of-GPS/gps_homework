clear; close all; clc;

currentFile = mfilename('fullpath');
currentFolder = fileparts(currentFile);
addpath(genpath(currentFolder + "/.."));

[signalData, samplesRead] = ifen_parser('gpsBase_IFEN_IF.bin');
fint = 5000445.88565834;
prn = genCA(7); 
pnr7(prn==0) = -1;
prn = -prn;

Tint = 0.001;
Tchip = (1e-3)/1023;
Fs = 20e6;
Ts = 1/Fs;
t = 0:Ts:(Tint-Ts);
spc = Ts/Tchip;

M = round(Tint/Ts);

dopRange = -10e3:500:10e3;
tauRange = 1:1023;

prnUp = upsample(prn, M);
batch = signalData(1:M)';
idx = 1;
for dop = dopRange
    replicaI = sin(2*pi*(fint + dop)*t);
    replicaQ = cos(2*pi*(fint + dop)*t);
    I = batch.*replicaI;
    Q = batch.*replicaQ;
    result = I + Q;
    fft_result = fft(result);
    fft_prn = fft(prnUp);
    conj_prn = conj(fft_prn);
    combined = fft_result.*conj_prn;
    corr(idx,:) = abs(ifft(combined)).^2;
    idx = idx+1;
end

figure();
surf(corr, 'EdgeColor', 'none');
title('PRN #7');
subtitle('1ms Integration Period');
xlabel('Doppler (Hz)');
ylabel('Code Phase (Chips)');

exportgraphics(gcf, currentFolder + "/../figures/p7_prn7.png", 'Resolution', 300);