%% PART I
clear; close all; clc;

Nmc = 1000;             % Number of Monte Carlo Runs

sigma_a = 0.5;          % Standard Deviation of a
sigma_b = 0.2;          % Standard Deviation of b

y_ = 0;                 % Expected Mean
sigma_y = sqrt(9*sigma_a^2 ...
    + 16*sigma_b^2);    % Expected Standard Devation

for i = 1:Nmc
    a = sigma_a*randn(1);   % State
    b = sigma_b*randn(1);   % State
    y(i) = 3*a - 4*b;       % Combined State
end
y_s = mean(y);           % Sample Mean
sigma_ys = std(y);       % Sample Standard Deviation

fprintf('Sample Mean from Monte Carlo: %0.3g\n', y_s);
fprintf('Sample Standard Deviation from Monte Carlo: %0.3g\n\n', sigma_ys);
fprintf('Expected Mean: %0.3g\n', y_);
fprintf('Expected Standard Deviation: %0.3g\n\n', sigma_y);

figure();
hold('on')
title('Histogram of Y');
histogram(y);
xlabel('Y-Value');
ylabel('Occurrences');

%% PART II
clear; close all; clc;

Nmc = 1000;

dt = 0.1;               % Time Step [s]
time = 0:dt:10;         % Time [s]
sigma1 = 0.1;           % Standard Deviation 1
sigma2 = 0.01;          % Standard Deviation 2

for i = 1:Nmc
    x1 = zeros(length(time),1);     % Integrated White Noise Vector
    x2 = zeros(length(time),1);     % Integrated White Noise Vector
    for j = 1:(length(time)-1)
        x1(j+1) = x1(j) + (sigma1*randn(1))*dt;
        x2(j+1) = x2(j) + (sigma2*randn(1))*dt;
    end
    x1_(i) = mean(x1);
    sigma_x1(i) = std(x1);
    x2_(i) = mean(x2);
    sigma_x2(i) = std(x2);

    if mod(i,250) == 0
        figure();
        hold('on');
        title(['Histogram of MC Run: ',num2str(i)]);
        histogram(x1);
        histogram(x2);
        xlabel('X-Value');
        ylabel('Occurrences');
    end
end

figure();
hold('on');
title('Integrated White Noise Mean Comparison');
plot(1:Nmc, x1_);
plot(1:Nmc, x2_);
xlabel('Monte Carlo Run');
ylabel('Mean (mu)');
legend('Sigma = 0.1', 'Sigma = 0.01');

figure();
hold('on');
title('Integrated White Noise Sigma Comparison');
plot(1:Nmc, sigma_x1);
plot(1:Nmc, sigma_x2);
xlabel('Monte Carlo Run');
ylabel('Standard Deviation (sigma)');
legend('Sigma = 0.1', 'Sigma = 0.01');

figure();
hold('on');
title('Final Monte Carlo Run');
plot(time,x1);
plot(time,x2);
xlabel('Time (s)');
ylabel('Integrated Noise');
legend('Sigma = 0.1', 'Sigma = 0.01');

tau1 = 1;               % Time Constant [s]
tau2 = 100;             % Time Constant [s]

for i = 1:Nmc
    x1t1 = zeros(length(time),1);     % 1st Order GM Vector; sig1; tau1
    x1t2 = zeros(length(time),1);     % 1st Order GM Vector; sig1; tau2
    x2t1 = zeros(length(time),1);     % 1st Order GM Vector; sig2; tau1
    x2t2 = zeros(length(time),1);     % 1st Order GM Vector; sig2; tau2
    for j = 1:length(time)-1
        x1t1_dot = -(1/tau1)*x1t1(j) + sigma1*randn(1);
        x1t2_dot = -(1/tau2)*x1t2(j) + sigma1*randn(1);
        x2t1_dot = -(1/tau1)*x2t1(j) + sigma2*randn(1);
        x2t2_dot = -(1/tau2)*x2t2(j) + sigma2*randn(1);
        
        x1t1(j+1) = x1t1(j) + x1t1_dot*dt;
        x1t2(j+1) = x1t2(j) + x1t2_dot*dt;
        x2t1(j+1) = x2t1(j) + x2t1_dot*dt;
        x2t2(j+1) = x2t2(j) + x2t2_dot*dt;
    end
    x1t1_(i) = mean(x1t1);
    sigma_x1t1(i) = std(x1t1);
    x1t2_(i) = mean(x1t2);
    sigma_x1t2(i) = std(x1t2);
    x2t1_(i) = mean(x2t1);
    sigma_x2t1(i) = std(x2t1);
    x2t2_(i) = mean(x2t2);
    sigma_x2t2(i) = std(x2t2);
end

figure();
hold('on');
title('1st Order Gauss Markov Mean Comparison');
plot(1:Nmc, x1t1_);
plot(1:Nmc, x1t2_);
plot(1:Nmc, x2t1_);
plot(1:Nmc, x2t2_);
xlabel('Monte Carlo Run');
ylabel('Mean (mu)');
legend('Sigma = 0.1; Tau = 1', 'Sigma = 0.1; Tau = 100', ...
    'Sigma = 0.01; Tau = 1', 'Sigma = 0.01; Tau = 100');

figure();
hold('on');
title('1st Order Gauss Markov Sigma Comparison');
plot(1:Nmc, sigma_x1t1);
plot(1:Nmc, sigma_x1t2);
plot(1:Nmc, sigma_x2t1);
plot(1:Nmc, sigma_x2t2);
xlabel('Monte Carlo Run');
ylabel('Standard Deviation (sigma)');
legend('Sigma = 0.1; Tau = 1', 'Sigma = 0.1; Tau = 100', ...
    'Sigma = 0.01; Tau = 1', 'Sigma = 0.01; Tau = 100');

%% PART V
clear; close all; clc;

Xs = [  0 300;
      100 400;
      700 400;
      800 300];
BS = [400 0];
Xu = [401 0];

rBu = norm(BS - Xu);