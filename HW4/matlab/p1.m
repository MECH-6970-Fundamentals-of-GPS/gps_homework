clear; close all; clc;

currentFile = mfilename('fullpath');
currentFolder = fileparts(currentFile);
addpath(genpath(currentFolder + "/.."));

time = 0:0.01:10;

s = tf('s');
plant = 1/s;
plant_eigs = pole(plant);
[plant_GM, plant_PM] = margin(plant);
[Y,~] = step(plant, time);
uc_ess = 1 - Y(end);

figure('Renderer', 'painters', 'Position', [10 10 900 600])
tiledlayout(2,2);
nexttile();
hold('on');
plot(time, Y, 'LineWidth', 2);
title('Uncontrolled Step Response');
xlabel('Time (s)');
ylabel('State');
ax = gca;
ax.FontSize = 18;

% (A) Proportional Controller
Kp = 2;
Pol = Kp*plant;
Pcl = Pol / (1 + Pol);
Pcl_eigs = pole(Pcl);
[P_GM, P_PM] = margin(Pol);
[YP,~] = step(Pcl, time);
P_ess = 1 - YP(end);

fprintf('A) Gain Margin: %0.4g\n', P_GM);
fprintf('A) Phase Margin: %0.4g\n', P_PM);
fprintf('A) EigenValues: %0.4g & %0.4g\n', Pcl_eigs);
fprintf('A) SS Error: %0.4g\n', P_ess);

nexttile();
hold('on');
plot(time, YP, 'LineWidth', 2);
title('Proportional Control Step Response');
xlabel('Time (s)');
ylabel('State');
ax = gca;
ax.FontSize = 18;

% (B)
[YP_ramp,~] = step(Pcl / s, time);
[Y_ramp,~] = step(1/s,time);

YP_ramp_ess = Y_ramp(end) - YP_ramp(end);
fprintf('B) SS Error: %0.4g\n', YP_ramp_ess);

nexttile();
hold('on');
plot(time, YP_ramp, 'LineWidth', 2);
plot(time, Y_ramp, '--', 'LineWidth', 2);
title('Proportional Control Ramp Response');
xlabel('Time (s)');
ylabel('State');
legend('System', 'R(t)');
ax = gca;
ax.FontSize = 18;

% (C)
Ki = 1;
PIol = ((Kp*s + Ki)/s)*plant;
PIcl = PIol / (1 + PIol);
PIcl_eigs = pole(PIcl);
[PI_GM, PI_PM] = margin(Pol);
[YPI_ramp,~] = step(PIcl / s, time);
PI_ramp_ess = Y_ramp(end) - YPI_ramp(end);


fprintf('C) Gain Margin: %0.4g\n', PI_GM);
fprintf('C) Phase Margin: %0.4g\n', PI_PM);
fprintf('C) EigenValues: %0.4g & %0.4g\n', PIcl_eigs);
fprintf('C) SS Error: %0.4g\n', PI_ramp_ess);

nexttile();
hold('on');
plot(time, YPI_ramp, 'LineWidth', 2);
plot(time, Y_ramp, '--', 'LineWidth', 2);
title('PI Control Ramp Response');
xlabel('Time (s)');
ylabel('State');
legend('System', 'R(t)');
ax = gca;
ax.FontSize = 18;

exportgraphics(gcf, currentFolder + "/../figures/p1.png", 'Resolution', 300);