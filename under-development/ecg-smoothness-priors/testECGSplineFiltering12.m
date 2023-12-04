clear
close all
clc

% load ECG1
% s = data;

load ECG2
s = data(14, 1:end);

% s = randn(1, 10000); fs = 1000;

% load ECG3
% s = data(2, 1:end);

% s = resample(s, 0.3*fs, fs);
% fs = 0.3*fs;
s = LPFilter(s - LPFilter(s, 1.0/fs), 100.0/fs);
s = s(:)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% warning('off', 'MATLAB:nearlySingularMatrix'); % THIS IS JUST FOR AVOIDING THE MATRIX SINGULARITY WARNINGS. A MORE NUMERICALLY STABLE VERSION OF THE CODE MIGHT BE NEEDED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
randn('seed', 0);
rand('seed', 0);

wlen = 100e-3; % window length (s)
f0 = 1.0; % Heart Beart (Hz)
DiffOrder = 2; % smoothness constraint order greater than 1
% lambda2 = 500*.5e-3*sqrt(fs^DiffOrder); % smoothness factor
% guardlen_l = 2*DiffOrder; % lower guard window lengh, >= DiffOrder
% guardlen_u = 2*DiffOrder; % upper guard window lengh, >= DiffOrder
snr = 10; % dB
% % % % adapt = 0; % adaptive(1), naive adaptation(2), or fixed(0) lambda
withwindow = 0; % with(1) or without(0) windowing
mode = 3; % 0..5

nvar = var(s)/10^(snr/10);
x = s + sqrt(nvar)*randn(size(s));

% lambda2 = 1e6*nvar;
% lambda2 = sqrt(nvar);
% [x_filtered1, x_filtered2] = ECGSmoothingPriors(x, DiffOrder, wlen*fs, lambda, guardlen_l, guardlen_u, adapt, withwindow);
% % % [x_filtered1, x_filtered2] = ECGSmoothnessPriorsDenoiser(x, DiffOrder, round(wlen*fs), 0, lambda2, adapt);
% [x_filtered1, x_filtered2] = ECGSmoothnessPriorsDenoiserBW(x, lambda2, mode, DiffOrder, round(wlen*fs));
[x_filtered1, x_filtered2, optim_gammas1, optim_gammas2] = ECGSmoothnessPriorsDenoiserBW(x, 10*(nvar), mode, DiffOrder, round(wlen*fs));
% [x_KF,x_KS,Pbar,Phat,PSmoothed,Kgain,a, x_HI, PPhat, KKgain, margin] = KalmanSmoothingPriorsFilter(x, .1*nvar);%, b, a, beta, var(x - BPFilter(x,4/fs,25/fs)), gamma, wlen, 'obsrv');


tt = (0:length(x)-1)/fs;

snr0 = 10*log10(sum(s.^2)/sum((s - x).^2));
snr1 = 10*log10(sum(s.^2)/sum((s - x_filtered1).^2));
snr2 = 10*log10(sum(s.^2)/sum((s - x_filtered2).^2));
% snr3 = 10*log10(sum(s.^2)/sum((s - x_KF).^2));
% snr4 = 10*log10(sum(s.^2)/sum((s - x_KS).^2));
% snr5 = 10*log10(sum(s.^2)/sum((s - x_HI).^2));

disp(['snr_input = ' num2str(snr0) '(dB)']);
disp(['snr_filtered1 = ' num2str(snr1) '(dB)']);
disp(['snr_filtered2 = ' num2str(snr2) '(dB)']);
% disp(['snr_KF = ' num2str(snr3) '(dB)']);
% disp(['snr_KS = ' num2str(snr4) '(dB)']);
% disp(['snr_HI = ' num2str(snr5) '(dB)']);

% peaks_0 = PeakDetection(s, f0/fs);
% I_0 = find(peaks_0);
% peaks_1 = PeakDetection(x_filtered1, f0/fs);
% I_1 = find(peaks_1);
% peaks_2 = PeakDetection(x_filtered2, f0/fs);
% I_2 = find(peaks_2);

% figure
% mesh(inv(eye(size(D,2)) + (lambda)^2*(D'*D)));
% mesh(inv(eye(size(Dd,2)) + (lambda)^2*(Dd'*Dd)));

% figure
% subplot(211);
% plot(tt, s, 'k');
% hold on
% plot(tt(I_1), s(I_1), 'bo');
% plot(tt(I_2), s(I_2), 'ro');
% grid
% subplot(212);
% HR_0 = 60*fs./diff(I_0);
% HR_1 = 60*fs./diff(I_1);
% HR_2 = 60*fs./diff(I_2);
% plot(tt(I_0(2:end)), HR_0, 'ko');
% hold on
% plot(tt(I_1(2:end)), HR_1, 'b');
% plot(tt(I_2(2:end)), HR_2, 'r');
% grid

% figure
% hold on
% plot(tt, x, 'c');
% plot(tt, s, 'k', 'linewidth', 2);
% plot(tt, x - x_filtered1, 'b');
% plot(tt, x - x_filtered2, 'r');
% grid

figure
hold on;
plot(tt, x, 'c');
plot(tt, s, 'k', 'linewidth', 2)
plot(tt, x_filtered1, 'b', 'linewidth', 2);
plot(tt, x_filtered2, 'r', 'linewidth', 2);
grid

figure
hold on
plot(optim_gammas1, 'b');
plot(optim_gammas2, 'r');
grid
% figure
% hold on
% plot(tt, x, 'c');
% plot(tt, s, 'k', 'linewidth', 2)
% plot(tt, x_KF, 'b', 'linewidth', 2);
% plot(tt, x_KS, 'r', 'linewidth', 2);
% plot(tt, x_HI, 'm', 'linewidth', 2);
% grid

% figure
% plot(a');
% grid
% 
% figure
% plot(Kgain');
% grid
% 
% figure
% hold on
% plot(squeeze(Pbar(1,1,:)), 'b');
% plot(squeeze(Phat(2,2,:)), 'r');
% plot(squeeze(PSmoothed(2,2,:)),'c');
% grid
% 
% figure
% hold on
% plot(squeeze(Pbar(1,1,:)), squeeze(Pbar(2,2,:)), 'b.');
% plot(squeeze(Phat(1,1,:)), squeeze(Phat(2,2,:)), 'r.');
% plot(squeeze(PSmoothed(1,1,:)), squeeze(PSmoothed(2,2,:)),'c.');
% grid
% 
% figure
% plot(margin');
% grid
