clear
close all
clc

% load ECG1
% s = data;

load ECG2
s = data(5, 1:2*fs);

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
% s = LPFilter(s - LPFilter(s, 2.0/fs), 100.0/fs);
s = s - LPFilter(s, 2.0/fs);%, 100.0/fs);
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
snr = 6; % dB
% % % % adapt = 0; % adaptive(1), naive adaptation(2), or fixed(0) lambda
withwindow = 0; % with(1) or without(0) windowing
mode = 0; % 0..5
gamma = 10;

nvar = var(s)/10^(snr/10);
noise =  NoiseGeneratorSimple('WHITE', nvar, fs, length(s));
% noise = noise.*cos(2*pi*(1:length(noise))*0.3/fs).*cos(2*pi*(1:length(noise))*0.21/fs);
% x = s + sqrt(nvar)*randn(size(s));
x = s + noise;

[x_filtered1, x_filtered2, optim_gammas1, optim_gammas2] = ECGSmoothnessPriorsDenoiserBW(x, gamma, mode, DiffOrder, [1 length(x)], 1e-8, 250, 1);
x_LTI = ECGSmoothnessPriorsDenoiserLTI(x, gamma, mode, DiffOrder);

tt = (0:length(x)-1)/fs;

snr0 = 10*log10(sum(s.^2)/sum((s - x).^2));
snr1 = 10*log10(sum(s.^2)/sum((s - x_filtered1).^2));
snr2 = 10*log10(sum(s.^2)/sum((s - x_filtered2).^2));
snr5 = 10*log10(sum(s.^2)/sum((s - x_LTI).^2));
mse1 = 100*mean((x_LTI - x_filtered1).^2)/mean(x_filtered1.^2)
mse2 = 100*mean((x_LTI - x_filtered2).^2)/mean(x_filtered2.^2)

disp(['snr_input = ' num2str(snr0) '(dB)']);
disp(['snr_filtered1 = ' num2str(snr1) '(dB)']);
disp(['snr_filtered2 = ' num2str(snr2) '(dB)']);
disp(['snr_LTI = ' num2str(snr5) '(dB)']);

bias = 0;
figure('position', [500, 200, 500, 500])  % create new figure with specified size  
hold on;
plot(tt, x, 'color', .5*ones(1,3), 'linewidth', 2);
plot(tt, s, 'k', 'linewidth', 2)
% plot(tt, x_filtered1 - bias, 'b', 'linewidth', 2);
plot(tt, x_filtered2 - 1*bias, 'r', 'linewidth', 2);
plot(tt, x_LTI - 2*bias, 'b', 'linewidth', 2);

bb = -bias*(0:5)'*ones(1, length(x));
plot(tt, bb, '--', 'color', .5*ones(1,3))
set(gca, 'YTickLabel', [])
set(gca, 'YTick', []);
set(gca, 'fontsize', 16);
xlabel('time(s)', 'fontsize', 16);
axis tight

figure
hold on
plot(tt, 10*log10(abs(x_filtered1 - x_LTI)));
plot(tt, 10*log10(abs(x_filtered2 - x_LTI)), 'r');
grid
