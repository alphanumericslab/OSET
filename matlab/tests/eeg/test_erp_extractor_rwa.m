% A script to test event related potential extraction using Robust Weighted Synchronous Averaging
%
% Reference:
%   Sameni, R. and Gouy-Pailler, C. (2014). An iterative subspace denoising
%   algorithm for removing electroencephalogram ocular artifacts. Journal of
%   neuroscience methods, 225, 97-105.
%
%   Sameni, R., Jutten, C., & Shamsollahi, M. B. (2010). A Deflation
%   Procedure for Subspace Decomposition. In IEEE Transactions on Signal
%   Processing (Vol. 58, Issue 4, pp. 2363–2374). Institute of Electrical
%   and Electronics Engineers (IEEE).
%   https://doi.org/10.1109/tsp.2009.2037353
%
% Reza Sameni, 2011-2024
% The Open-Source Electrophysiological Toolbox (OSET): https://github.com/alphanumericslab/OSET

clear
close all
clc

%% Select dataset or provide new ones as follows:
%   x: an EEG matrix of channels x time samples
%   fs: sampling frequency in Hz
%   target_indexes: target ERP indexes corresponding to the time indexes
%       when the target stimuli were given to the subject
%   nontarget_indexes: non-target ERP indexes corresponding to the time indexes
%       when the non-target stimuli were given to the subject
%   colnames: cell array of EEG column names used for plots (enumerates the channels if empty)

CASE_STUDY = 1; % change or add new ones
switch CASE_STUDY
    case 1
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
    case 2
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

%% Set algorithm parameters

N = size(x, 1); % number of channels
T = size(x, 2); % time length
erp_wlen = round(1.2*fs); % ERP window length
params = struct;
params.plot_results = false; % plot results or skip
% baseline removal method
params.baseline_filter = 'SINGLE-ORDER-IIR'; % 'NONE', 'DC', 'MDMN', 'SINGLE-ORDER-IIR'
switch params.baseline_filter
    case 'MDMN'
        params.baseline_wlen1 = 0.3;
        params.baseline_wlen2 = 0.45;
    case 'SINGLE-ORDER-IIR'
        params.hp_cuttoff = 0.1;
end

% ERP enhancer filter design (change filter specs or replace with custom design)
params.enhancer_filter = false; % apply or skip the filter
if params.enhancer_filter
    f_pass = 8.0;   % passband frequency
    f_stop = 9.5;   % stopband frequency
    passband_ripple = 0.1; % peak-to-peak passband ripple in dB
    d_pass = 10^((passband_ripple/2)/20) - 1;  % passband ripple
    stopband_attenuation = 30.0; % stopband attenuation in dB
    d_stop = 10^(-stopband_attenuation/20);     % stopband attenuation
    dens  = 20;               % density factor
    [NN, fo, Ao, WW] = firpmord([f_pass, f_stop]/(fs/2), [1 0], [d_pass, d_stop]); % calculate the order from the parameters using FIRPMORD.
    h_lp  = firpm(NN, fo, Ao, WW, {dens}); % calculate the coefficients using the FIRPM function.

    % set the numerator and denominator of the filter
    params.enhancer_filt_b = h_lp;
    params.enhancer_filt_a = 1;
end

% spatial filtering via nonstationary component analysis (NSCA)
params.enhancer_nsca = false; % apply or skip NSCA
if params.enhancer_nsca
    params.nsca_max_components = 6; % number of generalized eigenvalues of multichannel signals to preserve
end

% vertical offset adjuestment of each trial using the mean of the trial or a fraction at the beginnning of the trial
params.erp_vertical_offset_alignment = 'ZERO-MEAN-FRACTION'; % 'NONE', 'ZERO-MEAN', 'ZERO-MEAN-FRACTION'
if isequal(params.erp_vertical_offset_alignment , 'ZERO-MEAN-FRACTION')
    params.erp_vertical_offset_fraction = 0.15; % percentage of the beginning of each time window to use for baseline estimation and correction
end

%% Apply robust weighted averaging to the target and non-target trials
[erps_rwa_target, erps_mean_target, erps_rwmd_target, erps_median_target, erps_var_target, stacked_targets] = erp_extractor_rwa(x, fs, target_indexes, erp_wlen, params);
[erps_rwa_nontarget, erps_mean_nontarget, erps_rwmd_nontarget, erps_median_nontarget, erps_var_nontarget, stacked_nontargets] = erp_extractor_rwa(x, fs, nontarget_indexes, erp_wlen, params);

%% Plot the results
tt = (0:erp_wlen-1)/fs;
for i = 1:size(x, 1)
    figure
    hold on
    H = [];
    lgnd = {};
    % plot(tt, squeeze(stacked_targets(i, :, :))', 'color', [1, 0.6, 0.6]); % plot the stacked target events in shades
    % plot(tt, squeeze(stacked_nontargets(i, :, :))', 'color', [0.6, 0.6, 1]); % plot the stacked non-target events in shades
    h = plot(tt, erps_mean_target(i, :), 'r:', 'linewidth', 2); H = cat(2, H, h); lgnd = cat(2, lgnd, 'Mean Targets');
    h = plot(tt, erps_median_target(i, :), 'r-.', 'linewidth', 2); H = cat(2, H, h); lgnd = cat(2, lgnd, 'Median Targets');
    h = plot(tt, erps_rwa_target(i, :), 'r--', 'linewidth', 2); H = cat(2, H, h); lgnd = cat(2, lgnd, 'RWA Targets');
    h = plot(tt, erps_rwmd_target(i, :), 'r', 'linewidth', 2); H = cat(2, H, h); lgnd = cat(2, lgnd, 'RWAMD Targets');
    h = plot(tt, erps_mean_nontarget(i, :), 'b:', 'linewidth', 2); H = cat(2, H, h); lgnd = cat(2, lgnd, 'Mean Non-Targets');
    h = plot(tt, erps_median_nontarget(i, :), 'b-.', 'linewidth', 2); H = cat(2, H, h); lgnd = cat(2, lgnd, 'Median Non-Targets');
    h = plot(tt, erps_rwa_nontarget(i, :), 'b--', 'linewidth', 2); H = cat(2, H, h); lgnd = cat(2, lgnd, 'RWA Non-Targets');
    h = plot(tt, erps_rwmd_nontarget(i, :), 'b', 'linewidth', 2); H = cat(2, H, h); lgnd = cat(2, lgnd, 'RWAMD Non-Targets');
    grid
    legend(H, lgnd);
    set(gca, 'fontsize', 16)
    if isempty(colnames)
        title(['Channel ', num2str(i)]);
    else
        title(['Channel ', num2str(i), ': ', colnames{i}]);
    end
end
