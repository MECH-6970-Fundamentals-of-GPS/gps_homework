clear; close all; clc;

currentFile = mfilename('fullpath');
currentFolder = fileparts(currentFile);
addpath(genpath(currentFolder + "/.."));

% (A)
spc = 16;

prn4 = genCA(4);
prn4(prn4 == 0) = -1;
prn4 = -prn4;
prn4up = upsample(prn4, spc*1023);

tauRange = -5:(1/16):5;
idx = 1;
for shift = tauRange
    R(idx) = sum(prn4up.*circshift(prn4up, shift*spc));
    idx = idx + 1;
end

prn4_noisy = prn4up + 0.2*randn(1,spc*1023);

idx = 1;
for shift = tauRange
    R_noisy(idx) = sum(prn4_noisy.*circshift(prn4_noisy, shift*spc));
    idx = idx + 1;
end

figure();
hold('on');
plot(R, 'LineWidth', 2);
plot(R_noisy, '--', 'LineWidth', 2);
title('Autocorrelation');
subtitle('No Noise');
xlabel('Shift (chips)');
ylabel('Correlation');
legend('No Noise', '\sigma=0.2');

exportgraphics(gcf, currentFolder + "/../figures/p6_corr.png", 'Resolution', 300);