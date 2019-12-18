clear
close all
clc

% load ECG1
% s = data;

load ECG2
s = data(1, 1250:4250);

% s = randn(1, 10000); fs = 1000;

% load ECG3
% s = data(2, 1:end);

% addpath 'G:\Sameni\Documents\Papers\Journal 40 (Spline ECG Filtering)\CODES'
% pth = 'J:\ECGData\Arrhythmia';
% fs = 360;
% 
% wd = cd;
% cd(pth);
% fname = '233';
% segstart = 11.2 + .6;
% segstop = 13.2 + .375;
% system(['rdsamp -r ' fname ' -f ' num2str(segstart) ' -t ' num2str(segstop) ' > sampledata.txt']);
% data_all_channels = load([pth '\sampledata.txt'])';
% s = data_all_channels(2,:);

% s = resample(s, 0.3*fs, fs);
% fs = 0.3*fs;
s = LPFilter(s - LPFilter(s, 2.0/fs), 100.0/fs);
s = s(:)';

% figure
% plot(s')
% grid 

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
snr = 50; % dB
% % % % adapt = 0; % adaptive(1), naive adaptation(2), or fixed(0) lambda
withwindow = 0; % with(1) or without(0) windowing
% mode = 1; % 0..5
mode = 0; % 0..5
SmoothnessFactor = 1e1;

nvar = var(s)/10^(snr/10);
x = s + sqrt(nvar)*randn(size(s));

% lambda2 = 1e6*nvar;
% lambda2 = sqrt(nvar);
% [x_filtered1, x_filtered2] = ECGSmoothingPriors(x, DiffOrder, wlen*fs, lambda, guardlen_l, guardlen_u, adapt, withwindow);
% % % [x_filtered1, x_filtered2] = ECGSmoothnessPriorsDenoiser(x, DiffOrder, round(wlen*fs), 0, lambda2, adapt);
% [x_filtered1, x_filtered2] = ECGSmoothnessPriorsDenoiserBW(x, lambda2, mode, DiffOrder, round(wlen*fs));
[x_filtered1, x_filtered2, optim_gammas1, optim_gammas2] = ECGSmoothnessPriorsDenoiserBW(x, SmoothnessFactor, mode, DiffOrder, round(wlen*fs), 1e-8, 250, 1);
[xx_filtered1, xx_filtered2, xx_filtered3, ~, ~] = ECGSmoothnessPriorsDenoiserBWFiltFilt(x, SmoothnessFactor, mode, DiffOrder, round(wlen*fs), 1e-8, 250, 1);

% [x_KF,x_KS,Pbar,Phat,PSmoothed,Kgain,a, x_HI, PPhat, KKgain, margin] = ECGSmoothnessPriorsDenoiserKF(x, 2, 1*(nvar), 1);%, b, a, beta, var(x - BPFilter(x,4/fs,25/fs)), gamma, wlen, 'obsrv');
[x_KF,x_KS,Pbar,Phat,PSmoothed,Kgain,a] = ECGSmoothnessPriorsDenoiserKF(x, DiffOrder, 1*(nvar), .01*var(s), 1, 1, round(wlen*fs));%, b, a, beta, var(x - BPFilter(x,4/fs,25/fs)), gamma, wlen, 'obsrv');
x_LTI = ECGSmoothnessPriorsDenoiserLTI(x, .1*(nvar), [], DiffOrder);

tt = (0:length(x)-1)/fs;

snr0 = 10*log10(sum(s.^2)/sum((s - x).^2));
snr1 = 10*log10(sum(s.^2)/sum((s - x_filtered1).^2));
snr2 = 10*log10(sum(s.^2)/sum((s - x_filtered2).^2));
snr3 = 10*log10(sum(s.^2)/sum((s - x_KF).^2));
snr4 = 10*log10(sum(s.^2)/sum((s - x_KS).^2));
snr5 = 10*log10(sum(s.^2)/sum((s - x_LTI).^2));

disp(['snr_input = ' num2str(snr0) '(dB)']);
disp(['snr_filtered1 = ' num2str(snr1) '(dB)']);
disp(['snr_filtered2 = ' num2str(snr2) '(dB)']);
disp(['snr_KF = ' num2str(snr3) '(dB)']);
disp(['snr_KS = ' num2str(snr4) '(dB)']);
disp(['snr_LTI = ' num2str(snr5) '(dB)']);

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
hold on
plot(optim_gammas1, 'b');
plot(optim_gammas2, 'r');
grid

figure;%('position', [500, 200, 500, 500])  % create new figure with specified size  
hold on;
plot(tt, x, 'color', .5*ones(1,3), 'linewidth', 2);
plot(tt, s, 'k', 'linewidth', 3)
plot(tt, x_filtered1, 'b', 'linewidth', 3);
plot(tt, x_filtered2, 'r', 'linewidth', 2);
plot(tt, xx_filtered3, 'c', 'linewidth', 2);
plot(tt, xx_filtered3, 'b.');
set(gca, 'fontsize', 16);
xlabel('time(s)', 'fontsize', 16);
ylabel('Amplitude(mV)', 'fontsize', 16);
axis tight
set(gca, 'box', 'on')
grid

figure
hold on
plot(tt, abs(x_filtered2 - xx_filtered3), 'b', 'linewidth', 2);

% axis([0.4739    0.5996   -0.0940    0.0634]);

% set(gcf,'PaperPositionMode','auto');
% set(gcf, [426 201 550 581]);
% set(gca, 'Position', [0.13 0.011 0.775 0.815]);
% set(gca, 'OuterPosition', [0 1.38778e-017 1 1]);
% aa = axis;
% aa(1) = 0;
% aa(2) = tt(end);
% axis(aa);

% set(gca, 'YGrid', 'off')
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
