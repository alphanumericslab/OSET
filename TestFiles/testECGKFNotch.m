%
% Test program for ECG Kalman notch filtering
%
% Dependencies: The baseline wander toolbox of the Open Source ECG Toolbox
%
% Sample code for testing naive maternal template subtraction technique
%
% Reza Sameni (C)
% Email: rsameni@shirazu.ac.ir
% Web: www.sameni.info
%
% Crated 2007
% Modified June 2018
%
% Reference:
% Please cite this code as follows:
% 
% 1- R. Sameni, “A Linear Kalman Notch Filter for Power-Line Interference
% Cancellation,” in Proceedings of the 16th CSI International Symposium on 
% Artificial Intelligence and Signal Processing (AISP), Shiraz, Iran, 2-3 
% May 2012 2012, pp. 604-610. [Online]. 
% Available: https://doi.org/10.1109/AISP.2012.6313817
% 
% 2- R. Sameni, OSET: The open-source electrophysiological toolbox,
% URL http://www. oset. ir, v. 2012
% 

clc
clear all
close all;

load('SampleECG1kHz1.mat'); fs = 1000; % sample ECG 1
% load('SampleECG1kHz2.mat'); data = data(1:15000,6)'; fs = 1000; % sample ECG 2
% load('SamplePowerlineNoise1kHz'); data = sig(1:15000,1)'; fs = 1000; % sample powerline artifact
% load('SampleFECGDaISy'); data = cell2mat(data); data = data(1, :); fs = cell2mat(fs);

f0 = 50;

% baseline wander removal
data = data - LPFilter(data, 0.5/fs);

n = (0:length(data)-1);
x = data + .01*sin(2*pi*n*0.1/fs).*sin(2*pi*n*f0/fs)/std(data); % nonstationary powerline noise scenario 0
% x = data + (.055 + .02*sin(2*pi*n*0.1/fs)).*sin(2*pi*n*f0/fs)/std(data); % nonstationary powerline noise scenario 1
% x = data + (.055 + .02*sin(2*pi*n*1.0/fs)).*sin(2*pi*n*f0/fs)/std(data); % nonstationary powerline noise scenario 2 (notice the modulation around the notch frequency)
% x = data + 0.1*sin(2*pi*n*f0/fs + pi/9)/std(data); % stationary powerline noise
% x = data + .002*randn(size(data))/std(data); % no powerline

Q = 1e-6;
R = .1*var(x);
gamma = 0.99;

% standard implementation:
[KF, KS, Pbar, Phat, PSmoothed, Kgain] = KFNotch(x, f0, fs, Q, R, gamma);

% simplified implementation (works with the Kalman filter Q/R ratio):
% wlen = 10;
% [KF, KS, Pbar, Phat, PSmoothed, Kgain] = KFNotch2(x, f0, fs, 1000*Q/R, wlen); % w


% Display results
t = n/fs;

figure;
hold on;
plot(t, x, 'b');
plot(t, KF, 'm');
plot(t, KS, 'g');
plot(t, data, 'r');
plot(t, Kgain(:,:)', 'k');
grid;
xlabel('time(s)');
ylabel('Amplitude(mV)');
legend('noisy ECG', 'Kalman filter', 'Kalman smoother', 'Original ECG', 'Kalman Gain');
title('Powerline Cancellation using the Kalman Filter');

figure;
psd(x, length(x)/2,fs);
title('noisy spectrum');

figure;
psd(KF, length(KF)/2,fs);
title('Kalman filter output spectrum');

figure;
psd(KS, length(KS)/2,fs);
title('Kalman smoother output spectrum');

figure;
plot(t, [Kgain(1,:) ; Kgain(2,:)]);
grid
xlabel('time(s)');
ylabel('K');
title('The Kalman Gain');

figure;
hold on;
plot(t, sqrt(squeeze(Pbar(1,1,:))));
plot(t, sqrt(squeeze(Phat(1,1,:))), 'r');
plot(t, sqrt(squeeze(PSmoothed(1,1,:))), 'g');
legend('Pbar', 'Phat', 'Psmoothed');
grid
xlabel('time(s)');
ylabel('Magnitude(mV)');
title('Estimation Standard Deviation');

figure;
hold on;
plot(squeeze(Pbar(1,1,:)), squeeze(Pbar(2,2,:)), 'bo');
plot(squeeze(Phat(1,1,:)), squeeze(Phat(2,2,:)), 'ro');
plot(squeeze(PSmoothed(1,1,:)), squeeze(PSmoothed(2,2,:)), 'go');
legend('Pbar', 'Phat', 'Psmoothed');
grid
xlabel('P1');
ylabel('P2');
title('Estimation variance scatter plot (the less scattered, the better!)');
