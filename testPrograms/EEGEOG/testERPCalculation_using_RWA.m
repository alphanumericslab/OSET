% Event related potential extraction using Robust Weighted Averaging
%
% Reza Sameni (reza.sameni@gmail.com)
% Copyright June, 2019
%

clear
close all
clc

load BPFilter_1_20Hz.mat h
template = load('SampleEEGwithERPtemplatesRefMethod.txt')'; % Reference ERP templates provided by Dr. Marco Congedo (marco.congedo@gmail.com)
raw = load('SampleEEGwithERP.mat'); % Sample ERP data provided by Dr. Marco Congedo (marco.congedo@gmail.com)
x = raw.data(:, 1:16)';
trigs = raw.data(:, 17)';

target_indexes = find(trigs == 33285); % target indexes
nontarget_indexes = find(trigs == 33286); % non-target indexes

N = 16; % number of channels
T = size(x, 2); % time length
fs = 512; % sampling frequency
erp_wlen = (1.0*fs); % ERP window length
baseline_wlen1 = round(0.300*fs);
baseline_wlen2 = round(0.450*fs);

% % wavelet denoising parameters used for smoothing the average template beat
% TPTR = 'rigrsure';
% SORH = 's';
% SCAL = 'sln';
% NDEN = 7; % decrease to follow the average as it is, increase to make smoother
% WNAME = 'db4';%'coif5';

X_targets = zeros(N, erp_wlen, length(target_indexes));
X_nontargets = zeros(N, erp_wlen, length(nontarget_indexes));

% x_den = zeros(N, T);
% for ch = 1 : N
%     x_den(ch, :) = x(ch, :) - wden(x(ch, :), TPTR, SORH, SCAL, NDEN, WNAME);
% end

% x_den = LPFilter(x_den - LPFilter(x_den, 1.0/fs), 20.0/fs);

bl = BaseLine1(x, baseline_wlen1, 'md');
baseline = BaseLine1(bl, baseline_wlen2, 'mn');
xx = x - baseline;

% twlen = 0.3; % in seconds
% wwlen = round(twlen*fs);
% th = 1;
% M = 1;
% L = 3;
% x_den = EOGRemoval(xx, xx(1, :), wwlen, th, M, L, 0, 0);

% x_den = LPFilter(xx - LPFilter(xx, 1.0/fs), 20.0/fs);
x_den = zeros(N, T);
for i = 1 : N
    x_den(i, :) = filtfilt(h, 1, xx(i, :));
end

for k = 1 : length(target_indexes)
    start = target_indexes(k);
    segment = x_den(:, start : start + erp_wlen - 1);
    X_targets(:, :, k) = segment - repmat(mean(segment, 2), 1, erp_wlen);
end

for k = 1 : length(nontarget_indexes)
    start = nontarget_indexes(k);
    segment = x_den(:, start : start + erp_wlen - 1);
    X_nontargets(:, :, k) = segment - repmat(mean(segment, 2), 1, erp_wlen);
end

X_targets_mn = zeros(N, erp_wlen);
X_nontargets_mn = zeros(N, erp_wlen);
for i = 1 : N
    X_targets_mn(i, :) = RWAverage(squeeze(X_targets(i, :, :))');
    X_nontargets_mn(i, :) = RWAverage(squeeze(X_nontargets(i, :, :))');
    %     X_targets_mn(i, :) = mean(X_targets(i, :, :), 3);
    %     X_nontargets_mn(i, :) = mean(X_nontargets(i, :, :), 3);
end

% PlotECG(x, 4, 'b', fs, 'raw data');
% PlotECG(x_den, 4, 'r', fs, 'raw data');

t = (0:T-1)/fs;
plots_per_figure = 4;
for i = 1:N
    if(mod(i, plots_per_figure)==1 || plots_per_figure==1)
        figure;
    end
    subplot(plots_per_figure, 1, mod(i-1,plots_per_figure) + 1);
    plot(t, x(i,:),'b');
    hold on
    plot(t, xx(i,:),'m');
    plot(t, baseline(i,:),'g');
    plot(t, x_den(i,:),'r');
    ylabel(num2str(i));
    grid;
    a = axis;
    %     axis([a(1) a(2) -50 +50]);
    %     axis([a(1) a(2) -1e-4 1e-4]);
    if(mod(i,plots_per_figure)==1 || plots_per_figure==1)
        title('Noisy vs. Denoised');
    end
    if(mod(i,plots_per_figure)==0 || plots_per_figure==1)
        xlabel('time(s)');
    end
end

tt = (0:erp_wlen-1)/fs;
for i = 1:N
    figure
    %     periodogram(x(i, :),[],'twosided', 1024, fs);
    %     plot(tt, squeeze(X_targets(i, :, :))', 'b');
    hold on
    %     plot(tt, squeeze(X_nontargets(i, :, :))', 'r');
    plot(tt, X_targets_mn(i, :), 'k', 'linewidth', 2);
    plot(tt, X_nontargets_mn(i, :), 'g', 'linewidth', 2);
    plot((0:length(template(i, :))-1)/128.0, std(X_targets_mn(i, :))*template(i, :)/std(template(i, :)), 'c', 'linewidth', 2);
    grid
    a = axis;
    %     axis([a(1) a(2) -50 +50]);
    axis([a(1) a(2) -0.5e-5 0.5e-5]);
    %     axis([a(1) a(2) -0.5 0.5]);
    title(['Channel: ', num2str(i)]);
    legend('Targets (proposed method)', 'Non-targets', 'Targets (reference method)');
end

