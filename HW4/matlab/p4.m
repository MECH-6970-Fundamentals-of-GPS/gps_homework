clear; close all; clc;

currentFile = mfilename('fullpath');
currentFolder = fileparts(currentFile);
addpath(genpath(currentFolder + "/.."));

signal = generate_signal(3);
prn = [1 -1 -1 -1 1 -1 1 1];
spc = 10;
Tchip = 1e-3/1023;
Ts = Tchip*spc;
time = 0:Ts:Ts*(length(signal)-2);

wn = 0.5*2*pi;
zeta = 0.9;
Ki = wn^2;
Kp = 2*zeta*wn;

N = length(signal);
prnUp = upsample(prn, spc*length(prn));

M = length(prnUp);
shift = -M:M;
s = 0.5*spc;
Tint = M/spc*Tchip;

tau = ones(1,N/M);
err = zeros(1,N/M);
ierr = zeros(1,N/M);

figure();
for k = 1:N/M-1
    clf;
    hold('on');

    X = signal(1 + (k-1)*M : k*M);
    R = scorr(X, prnUp);
    idx = find(tau(k)==shift);
    RE = R(idx+s);
    RL = R(idx-s);
    err(k+1) = (RE - RL)/(RE + RL);
    ierr(k+1) = ierr(k) + err(k+1)*Tint;
    tau(k+1) = tau(k) + round(Kp*err(k+1) + Ki*ierr(k+1));
    prnUp = circshift(prnUp, tau(k+1));

    plot(X);
    plot(prnUp);
    pause(0.1);
end

figure();
plot(time(1:M:end), tau, 'LineWidth', 2);
title('Code Phase vs. Time');
xlabel('Time (s)');
ylabel('Code Phase');

exportgraphics(gcf, currentFolder + "/../figures/p4_tau.png", 'Resolution', 300);

%% FUNCTIONS

function [CE, CP, CL] = shift_code(code, shift, spc, spacing)
%shift_code.m Shift C/A Code
%   Shifts the C/A code by a specified amount
% Inputs:
%    code       : PRN Code (already upsampled)
%    shift      : shift to the code (in terms of chips)
%    spc        : samples per chip
%    spacing    : I don't remember
%
% Outputs:
%    CE         : C/A Early
%    CP         : C/A Prompt
%    CL         : C/A Late
%
% Author: Walter Livingston

    if ~exist('spc', 'var')
        spc = 1;
    end
    if ~exist('spacing', 'var')
        spacing = 1/2;
    end
    mult = round(spacing*spc);
    tau = mult*shift;
    CE = circshift(code, tau);
    CP = code;
    CL = circshift(code, -tau);
end