% Event related potential extraction using Robust Weighted Averaging
%
% Reza Sameni (reza.sameni@gmail.com)
% Copyright June, 2019
%

clear
close all
clc

if 1
    % dat = edfread('/Users/rsameni/Documents/erp-sample/sub-001_task-P300trainrun1_run-6_eeg_EXPORTED.edf');
    dat = edfread('/Users/rsameni/Documents/erp-sample/sub-001_task-P300trainrun2_run-7_eeg_EXPORTED.edf');
    % dat = edfread('/Users/rsameni/Documents/erp-sample/sub-001_task-RSVPtask_run-3_eeg_EXPORTED.edf');
    fs = 512; % sampling frequency
    dat = timetable2table(dat);
    x = cell2mat(table2array(dat(:, 2:end-1)))';
    colnames = dat.Properties.VariableNames(2:end-1);
    trigs = cell2mat(table2array(dat(:, end)))';
    target_indexes = find(trigs >= 0.9 & trigs < 1.1); % target indexes
    nontarget_indexes = find(trigs >= 1.9 & trigs < 2.1); % non-target indexes
else
    load BPFilter_1_20Hz.mat h
    template = load('SampleEEGwithERPtemplatesRefMethod.txt')'; % Reference ERP templates provided by Dr. Marco Congedo (marco.congedo@gmail.com)
    raw = load('SampleEEGwithERP.mat'); % Sample ERP data provided by Dr. Marco Congedo (marco.congedo@gmail.com)
    fs = 512; % sampling frequency
    colnames = [];
    T = 110;
    x = raw.data(1:round(T*fs)-1, 1:16)';
    trigs = raw.data(1:round(T*fs)-1, 17)';
    % x = raw.data(:, 1:16)';
    % trigs = raw.data(:, 17)';

    target_indexes = find(trigs == 33285); % target indexes
    nontarget_indexes = find(trigs == 33286); % non-target indexes
end

N = size(x, 1); % number of channels
T = size(x, 2); % time length
erp_wlen = round(0.6*fs); % ERP window length

% spikes = zeros(1, T);
% spikes(target_indexes) = 1;
% pulse = filter(ones(1, erp_wlen), 1, spikes);
% pulse = find(pulse > 0);
% [x_decomposed, W, A, B, C, lambda] = nonstationary_component_analysis(x, pulse, 1:T);

params = struct;

params.plot_results = false;

params.baseline_filter = 'SINGLE-ORDER-IIR'; % 'NONE', 'DC', 'MDMN', 'SINGLE-ORDER-IIR'
switch params.baseline_filter
    case 'MDMN'
        params.baseline_wlen1 = 0.3;
        params.baseline_wlen2 = 0.45;
    case 'SINGLE-ORDER-IIR'
        params.hp_cuttoff = 0.1;
end

params.enhancer_filter = true;
if params.enhancer_filter
    % load lp_filter_coefs.mat h_lp
    % load lp_filter_coefs_512_9hz h_lp
    f_pass = 8.0;                % passband frequency
    f_stop = 9.5;              % stopband frequency
    passband_ripple = 0.1; % peak-to-peak passband ripple in dB
    d_pass = 10^((passband_ripple/2)/20) - 1;  % passband ripple
    stopband_attenuation = 30.0; % stopband attenuation in dB
    d_stop = 10^(-stopband_attenuation/20);     % stopband attenuation
    dens  = 20;               % density factor
    % calculate the order from the parameters using FIRPMORD.
    [NN, fo, Ao, WW] = firpmord([f_pass, f_stop]/(fs/2), [1 0], [d_pass, d_stop]);

    % calculate the coefficients using the FIRPM function.
    h_lp  = firpm(NN, fo, Ao, WW, {dens});

    params.enhancer_filt_b = h_lp;
    params.enhancer_filt_a = 1;
end

params.enhancer_nsca = true;
if params.enhancer_nsca
    params.nsca_max_components = 6;
end

params.erp_vertical_offset_alignment = 'ZERO-MEAN-FRACTION'; % 'NONE', 'ZERO-MEAN', 'ZERO-MEAN-FRACTION'
if isequal(params.erp_vertical_offset_alignment , 'ZERO-MEAN-FRACTION')
    params.erp_vertical_offset_fraction = 0.15;
end

[erps_rwa_target, erps_mean_target, erps_rwmd_target, erps_var_target, stacked_targets] = erp_extractor_rwa(x, fs, target_indexes, erp_wlen, params);
[erps_rwa_nontarget, erps_mean_nontarget, erps_rwmd_nontarget, erps_var_nontarget, stacked_nontargets] = erp_extractor_rwa(x, fs, nontarget_indexes, erp_wlen, params);


tt = (0:erp_wlen-1)/fs;
for i = 1:size(x, 1)
    figure
    hold on
    % plot(tt, squeeze(stacked_targets(i, :, :))', 'color', [1, 0.6, 0.6]);
    % plot(tt, squeeze(stacked_nontargets(i, :, :))', 'color', [0.6, 0.6, 1]);
    h1 = plot(tt, erps_mean_target(i, :), 'r:', 'linewidth', 2);
    h2 = plot(tt, erps_rwa_target(i, :), 'r--', 'linewidth', 2);
    h3 = plot(tt, erps_rwmd_target(i, :), 'r', 'linewidth', 2);
    h4 = plot(tt, erps_mean_nontarget(i, :), 'b:', 'linewidth', 2);
    h5 = plot(tt, erps_rwa_nontarget(i, :), 'b--', 'linewidth', 2);
    h6 = plot(tt, erps_rwmd_nontarget(i, :), 'b', 'linewidth', 2);
    grid
    legend([h1, h2, h3, h4, h5, h6], 'Mean Targets', 'RWA Targets', 'RWMD Targets', 'Mean Non-targets', 'RWA Non-targets', 'RWMD Non-targets');
    set(gca, 'fontsize', 16)
    if isempty(colnames)
        title(['Channel ', num2str(i)]);
    else
        title(['Channel ', num2str(i), ': ', colnames{i}]);
    end
end

% % % % % % % if 0
% % % % % % % 
% % % % % % %     % % wavelet denoising parameters used for smoothing the average template beat
% % % % % % %     % TPTR = 'rigrsure';
% % % % % % %     % SORH = 's';
% % % % % % %     % SCAL = 'sln';
% % % % % % %     % NDEN = 7; % decrease to follow the average as it is, increase to make smoother
% % % % % % %     % WNAME = 'db4';%'coif5';
% % % % % % % 
% % % % % % %     X_targets = zeros(N, erp_wlen, length(target_indexes));
% % % % % % %     X_nontargets = zeros(N, erp_wlen, length(nontarget_indexes));
% % % % % % % 
% % % % % % %     % x_den = zeros(N, T);
% % % % % % %     % for ch = 1 : N
% % % % % % %     %     x_den(ch, :) = x(ch, :) - wden(x(ch, :), TPTR, SORH, SCAL, NDEN, WNAME);
% % % % % % %     % end
% % % % % % % 
% % % % % % %     % x_den = LPFilter(x_den - LPFilter(x_den, 1.0/fs), 20.0/fs);
% % % % % % % 
% % % % % % %     bl = baseline_sliding_window(x, baseline_wlen1, 'md');
% % % % % % %     baseline = baseline_sliding_window(bl, baseline_wlen2, 'mn');
% % % % % % %     xx = x - baseline;
% % % % % % % 
% % % % % % %     % twlen = 0.3; % in seconds
% % % % % % %     % wwlen = round(twlen*fs);
% % % % % % %     % th = 1;
% % % % % % %     % M = 1;
% % % % % % %     % L = 3;
% % % % % % %     % x_den = EOGRemoval(xx, xx(1, :), wwlen, th, M, L, 0, 0);
% % % % % % % 
% % % % % % %     % x_den = LPFilter(xx - LPFilter(xx, 1.0/fs), 20.0/fs);
% % % % % % % 
% % % % % % %     x_den = zeros(N, T);
% % % % % % %     for i = 1 : N
% % % % % % %         x_den(i, :) = filtfilt(h, 1, xx(i, :));
% % % % % % %     end
% % % % % % % 
% % % % % % % 
% % % % % % %     % x_den = xx;
% % % % % % % 
% % % % % % % 
% % % % % % %     for k = 1 : length(target_indexes)
% % % % % % %         start = target_indexes(k);
% % % % % % %         segment = x_den(:, start : start + erp_wlen - 1);
% % % % % % %         X_targets(:, :, k) = segment - repmat(mean(segment, 2), 1, erp_wlen);
% % % % % % %     end
% % % % % % % 
% % % % % % %     for k = 1 : length(nontarget_indexes)
% % % % % % %         start = nontarget_indexes(k);
% % % % % % %         segment = x_den(:, start : start + erp_wlen - 1);
% % % % % % %         X_nontargets(:, :, k) = segment - repmat(mean(segment, 2), 1, erp_wlen);
% % % % % % %     end
% % % % % % % 
% % % % % % %     X_targets_mn = zeros(N, erp_wlen);
% % % % % % %     X_nontargets_mn = zeros(N, erp_wlen);
% % % % % % %     for i = 1 : N
% % % % % % %         X_targets_mn(i, :) = robust_weighted_average(squeeze(X_targets(i, :, :))');
% % % % % % %         X_nontargets_mn(i, :) = robust_weighted_average(squeeze(X_nontargets(i, :, :))');
% % % % % % %         %     X_targets_mn(i, :) = mean(X_targets(i, :, :), 3);
% % % % % % %         %     X_nontargets_mn(i, :) = mean(X_nontargets(i, :, :), 3);
% % % % % % %     end
% % % % % % % 
% % % % % % %     % PlotECG(x, 4, 'b', fs, 'raw data');
% % % % % % %     % PlotECG(x_den, 4, 'r', fs, 'raw data');
% % % % % % % 
% % % % % % %     t = (0:T-1)/fs;
% % % % % % %     plots_per_figure = 4;
% % % % % % %     for i = 1:N
% % % % % % %         if(mod(i, plots_per_figure)==1 || plots_per_figure==1)
% % % % % % %             figure;
% % % % % % %         end
% % % % % % %         subplot(plots_per_figure, 1, mod(i-1,plots_per_figure) + 1);
% % % % % % %         plot(t, x(i,:),'b');
% % % % % % %         hold on
% % % % % % %         plot(t, xx(i,:),'m');
% % % % % % %         plot(t, baseline(i,:),'g');
% % % % % % %         plot(t, x_den(i,:),'r');
% % % % % % %         ylabel(num2str(i));
% % % % % % %         grid;
% % % % % % %         a = axis;
% % % % % % %         %     axis([a(1) a(2) -50 +50]);
% % % % % % %         %     axis([a(1) a(2) -1e-4 1e-4]);
% % % % % % %         if(mod(i,plots_per_figure)==1 || plots_per_figure==1)
% % % % % % %             title('Noisy vs. Denoised');
% % % % % % %         end
% % % % % % %         if(mod(i,plots_per_figure)==0 || plots_per_figure==1)
% % % % % % %             xlabel('time(s)');
% % % % % % %         end
% % % % % % %     end
% % % % % % % 
% % % % % % %     tt = (0:erp_wlen-1)/fs;
% % % % % % %     for i = 1:N
% % % % % % %         figure
% % % % % % %         %     periodogram(x(i, :),[],'twosided', 1024, fs);
% % % % % % %         %     plot(tt, squeeze(X_targets(i, :, :))', 'b');
% % % % % % %         hold on
% % % % % % %         %     plot(tt, squeeze(X_nontargets(i, :, :))', 'r');
% % % % % % %         plot(tt, X_targets_mn(i, :), 'k', 'linewidth', 2);
% % % % % % %         plot(tt, X_nontargets_mn(i, :), 'g', 'linewidth', 2);
% % % % % % %         plot((0:length(template(i, :))-1)/128.0, std(X_targets_mn(i, :))*template(i, :)/std(template(i, :)), 'c', 'linewidth', 2);
% % % % % % %         grid
% % % % % % %         a = axis;
% % % % % % %         %     axis([a(1) a(2) -50 +50]);
% % % % % % %         axis([a(1) a(2) -0.5e-5 0.5e-5]);
% % % % % % %         %     axis([a(1) a(2) -0.5 0.5]);
% % % % % % %         title(['Channel: ', num2str(i)]);
% % % % % % %         legend('Targets (proposed method)', 'Non-targets', 'Targets (reference method)');
% % % % % % %     end
% % % % % % % 
% % % % % % % end
