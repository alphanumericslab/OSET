% A script to test EOG cancellation from EEG followed by event related potential extraction using Robust Weighted Synchronous Averaging
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
% Reza Sameni, 2019-2024
% The Open-Source Electrophysiological Toolbox (OSET): https://github.com/alphanumericslab/OSET

clear
close all
clc

CASE_STUDY = 1;
switch CASE_STUDY
    case 1
        % dat = edfread('/Users/rsameni/Documents/erp-sample/sub-001_task-P300trainrun1_run-6_eeg_EXPORTED.edf');
        % dat = edfread('/Users/rsameni/Documents/erp-sample/sub-001_task-P300trainrun2_run-7_eeg_EXPORTED.edf');
        dat = edfread('/Users/rsameni/Documents/erp-sample/sub-001_task-RSVPtask_run-3_eeg_EXPORTED.edf');
        fs = 512; % sampling frequency
        eog_channel = 1; % FP1
        dat = timetable2table(dat);
        x = cell2mat(table2array(dat(:, 2:end-1)))';
        colnames = dat.Properties.VariableNames(2:end-1);
        trigs = cell2mat(table2array(dat(:, end)))';
        target_indexes = find(trigs >= 0.9 & trigs < 1.1); % target indexes
        nontarget_indexes = find(trigs >= 1.9 & trigs < 2.1); % non-target indexes
    case 2
        % load BPFilter_1_20Hz.mat h
        template = load('SampleEEGwithERPtemplatesRefMethod.txt')'; % Reference ERP templates provided by Dr. Marco Congedo (marco.congedo@gmail.com)
        raw = load('SampleEEGwithERP.mat'); % Sample ERP data provided by Dr. Marco Congedo (marco.congedo@gmail.com)
        fs = 512; % sampling frequency
        colnames = [];
        eog_channel = 1; % First EEG channel
        T = 110;
        x = raw.data(1:round(T*fs)-1, 1:16)';
        trigs = raw.data(1:round(T*fs)-1, 17)';
        % x = raw.data(:, 1:16)';
        % trigs = raw.data(:, 17)';
        target_indexes = find(trigs == 33285); % target indexes
        nontarget_indexes = find(trigs == 33286); % non-target indexes
end

%% APPLY PREPROCESSING FILTER (change filter design or replace with cutom filter as needed)
params.apply_preprocessing = true;
if params.apply_preprocessing
    fc1 = 0.1; % lower frequency cutoff
    x = x - lp_filter_zero_phase(x, fc1/fs);
    fc2 = 50; % upper frequency cutoff
    x = lp_filter_zero_phase(x, fc2/fs);
end

%% Select a subset of trials if needed
if 0 % set to false to process all trials or change the trial indexes as needed
    target_indexes = target_indexes(1:50);
    nontarget_indexes = nontarget_indexes(1:50);
end

%% REMOVE EOG ARTIFACTS
% Set eog_canceller_nsca parameters
twlen = 0.3; % Time window length in seconds for EOG power envelope calculation
num_itr = 3; % Number of iterations for the artifact removal process
th_method = 'envelope-prctile'; % Threshold method for EOG segment detection. Options: 'absolute', 'envelope-mean-fraction', 'envelope-median-fraction', 'envelope-prctile'
th_value = 75.0; % Threshold value for the chosen method. Interpretation depends on th_method. In 'envelope-prctile' mode, this is the envelope amplitude percentile above which is EOG
num_comp_est_method = 'opt'; % Method for estimating the number of components to denoise in each iteration. Options: 'fixed', 'opt'
shrinkage_params = struct; % Initialize shrinkage parameters structure
switch num_comp_est_method % Configure parameters based on the component estimation method
    case 'fixed'
        % Use a fixed number of channels based on the largest generalized eigenvalues
        shrinkage_params.num_noise_gevds = 1;
    case 'opt'
        % Dynamically estimate the EOG dimensions embedded in background EEG noise
        % Additional parameters can be set here if required for dynamic estimation
end

shrinkage_method = 'wden'; % Choose the method for EOG component shrinkage. Options: 'remove', 'wden' for wavelet denoising
switch shrinkage_method % Configure shrinkage parameters based on the chosen method
    case 'remove'
        % Simply remove components without filtering, reducing dimensionality
        % This mode expects shrinkage_params.num_noise_gevds to be set
    case 'wden'
        % Wavelet denoising parameters
        shrinkage_params.TPTR = 'heursure'; % Threshold selection rule
        shrinkage_params.SORH = 's';        % Soft or hard thresholding
        shrinkage_params.SCAL = 'mln';      % Level-dependent threshold scaling
        shrinkage_params.NUM = 5;           % Wavelet decomposition level
        shrinkage_params.WNAME = 'sym5';    % Wavelet name
end
adapt_reference = 'no'; % Adaptation method for the reference EOG signal across iterations. Options: 'no', 'top-gevd', 'top-few-gevd'
flagplot = 1; % Plotting flag to control the verbosity of output plots

%% Call the EOG artifact removal function with configured parameters
eog = x(eog_channel, :); % set the EOG channel
x_post_eog_canceller = eog_canceller_nsca(x, eog, fs, twlen, num_itr, th_method, th_value, num_comp_est_method, shrinkage_method, shrinkage_params, adapt_reference, flagplot);

%% Plot the EOG cancellation results
t = (0:length(x)-1)/fs;
L = 4;
ch = 1:size(x,1);
for i = 1 : length(ch)
    if mod(i,L)==1
        figure;
    end
    subplot(L,1,mod(i-1,L)+1);
    plot(t,x(ch(i),:),'k');
    hold on;
    plot(t,x_post_eog_canceller(ch(i),:),'color',.6*ones(1,3));
    % % %     grid;
    set(gca,'Box','On','FontSize',16);
    if i < L
        set(gca,'XTickLabel',[]);
    end
    ylabel(['EEG',num2str(i)]);
    axis tight
end
xlabel('time(s)','FontSize',16);


%% Extract the ERP for target vs non-targets
% set up the ERP detection algorithm parameters
N = size(x, 1); % number of channels
T = size(x, 2); % time length
erp_wlen = round(0.6*fs); % ERP window length
params = struct;
params.plot_results = false;
params.baseline_filter = 'SINGLE-ORDER-IIR'; % 'NONE', 'DC', 'MDMN', 'SINGLE-ORDER-IIR'
switch params.baseline_filter
    case 'MDMN'
        params.baseline_wlen1 = 0.3;
        params.baseline_wlen2 = 0.45;
    case 'SINGLE-ORDER-IIR'
        % params.hp_cuttoff = 0.5;
        params.hp_cuttoff = 0.1; % used for single trial tests (under-dev.); keeps the low-frequency trend of the signal
end

% ERP enhancer filter design (change filter specs or replace with custom design)
% ISSUE: Applying causal filters changes the trigger and data latencies;
%   non-causal filters (by forward-backward filtfilt) do not result in such
%   lags, but make the ERP activity non-causal (appearing before their actual
%   expected times.
% TODO: All temporal filtering should be converted to causal FIR filters
%   and their group-delays with the triggers (half filter length) should be
%   compensated before synchronous averaging.
params.enhancer_filter = false; % apply a preprocessing filter of skip it
if params.enhancer_filter
    f_pass = 8.0;                % passband frequency
    % f_pass = 15.0;                % passband frequency
    f_stop = 9.5;              % stopband frequency
    % f_stop = 18.5;              % stopband frequency
    passband_ripple = 0.1; % peak-to-peak passband ripple in dB
    d_pass = 10^((passband_ripple/2)/20) - 1;  % passband ripple
    stopband_attenuation = 30.0; % stopband attenuation in dB
    d_stop = 10^(-stopband_attenuation/20);     % stopband attenuation
    dens  = 20;               % density factor
    [NN, fo, Ao, WW] = firpmord([f_pass, f_stop]/(fs/2), [1 0], [d_pass, d_stop]); % calculate the order from the parameters using FIRPMORD.
    h_lp  = firpm(NN, fo, Ao, WW, {dens}); % calculate the coefficients using the FIRPM function.
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
    params.erp_vertical_offset_fraction = 0.1; % percentage of the beginning of each time window to use for baseline estimation and correction
end

%% Apply robust weighted averaging to the target and non-target trials on the original EEG and the EOG-removed data
% perform the calculations on a rescaled version of the post-ECG
% cancellation signals to make the ERP amplitudes comparable before and
% after EOG removal
scale = std(x, [], 2) ./ std(x_post_eog_canceller, [], 2);
x_post_eog_canceller_norm = (x_post_eog_canceller - mean(x_post_eog_canceller, 2)).*(scale * ones(1, T));
[erps_rwa_target, erps_mean_target, erps_rwmd_target, erps_median_target, erps_var_target, stacked_targets] = erp_extractor_rwa(x, fs, target_indexes, erp_wlen, params);
[erps_rwa_nontarget, erps_mean_nontarget, erps_rwmd_nontarget, erps_median_nontarget, erps_var_nontarget, stacked_nontargets] = erp_extractor_rwa(x, fs, nontarget_indexes, erp_wlen, params);
[erps_rwa_target_eog_rem, erps_mean_target_eog_rem, erps_rwmd_target_eog_rem, erps_median_target_eog_rem, erps_var_target_eog_rem, stacked_targets_eog_rem] = erp_extractor_rwa(x_post_eog_canceller_norm, fs, target_indexes, erp_wlen, params);
[erps_rwa_nontarget_eog_rem, erps_mean_nontarget_eog_rem, erps_rwmd_nontarget_eog_rem, erps_median_nontarget_eog_rem, erps_var_nontarget_eog_rem, stacked_nontargets_eog_rem] = erp_extractor_rwa(x_post_eog_canceller_norm, fs, nontarget_indexes, erp_wlen, params);

%% Plot the results
tt = (0:erp_wlen-1)/fs;
if 1
    for i = 1:size(x, 1)
        figure
        hold on
        % plot(tt, squeeze(stacked_targets(i, :, :))', 'color', [1, 0.6, 0.6]);
        % plot(tt, squeeze(stacked_nontargets(i, :, :))', 'color', [0.6, 0.6, 1]);
        H = [];
        lgnd = {};
        h = plot(tt, erps_rwa_target(i, :), 'm--', 'linewidth', 2); H = cat(2, H, h); lgnd = cat(2, lgnd, 'RWA Targets');
        h = plot(tt, erps_mean_target(i, :), 'm:', 'linewidth', 2); H = cat(2, H, h); lgnd = cat(2, lgnd, 'Mean Targets');
        h = plot(tt, erps_rwmd_target(i, :), 'm', 'linewidth', 2); H = cat(2, H, h); lgnd = cat(2, lgnd, 'RWMD Targets');
        % h = plot(tt, cumsum(erps_rwmd_target(i, :)), 'm-.', 'linewidth', 4); H = cat(2, H, h); lgnd = cat(2, lgnd, 'RWMD Targets Cumsum');
        h = plot(tt, erps_rwa_nontarget(i, :), 'c--', 'linewidth', 2); H = cat(2, H, h); lgnd = cat(2, lgnd, 'RWA Non-targets');
        h = plot(tt, erps_mean_nontarget(i, :), 'c:', 'linewidth', 2); H = cat(2, H, h); lgnd = cat(2, lgnd, 'Mean Non-targets');
        h = plot(tt, erps_rwmd_nontarget(i, :), 'c', 'linewidth', 2); H = cat(2, H, h); lgnd = cat(2, lgnd, 'RWMD Non-targets');
        % h = plot(tt, cumsum(erps_rwmd_nontarget(i, :)), 'c-.', 'linewidth', 4); H = cat(2, H, h); lgnd = cat(2, lgnd, 'RWMD Non-targets Cumsum');
        h = plot(tt, erps_rwa_target_eog_rem(i, :), 'r--', 'linewidth', 2); H = cat(2, H, h); lgnd = cat(2, lgnd, 'RWA Targets (post EOG removal)');
        h = plot(tt, erps_mean_target_eog_rem(i, :), 'r:', 'linewidth', 2); H = cat(2, H, h); lgnd = cat(2, lgnd, 'Mean Targets (post EOG removal)');
        h = plot(tt, erps_rwmd_target_eog_rem(i, :), 'r', 'linewidth', 2); H = cat(2, H, h); lgnd = cat(2, lgnd, 'RWMD Targets (post EOG removal)');
        % h = plot(tt, cumsum(erps_rwmd_target_eog_rem(i, :)), 'r-.', 'linewidth', 4); H = cat(2, H, h); lgnd = cat(2, lgnd, 'RWMD Targets Cumsum (post EOG removal)');
        h = plot(tt, erps_rwa_nontarget_eog_rem(i, :), 'b--', 'linewidth', 2); H = cat(2, H, h); lgnd = cat(2, lgnd, 'RWA Non-targets (post EOG removal)');
        h = plot(tt, erps_mean_nontarget_eog_rem(i, :), 'b:', 'linewidth', 2); H = cat(2, H, h); lgnd = cat(2, lgnd, 'Mean Non-targets (post EOG removal)');
        h = plot(tt, erps_rwmd_nontarget_eog_rem(i, :), 'b', 'linewidth', 2); H = cat(2, H, h); lgnd = cat(2, lgnd, 'RWMD Non-targets (post EOG removal)');
        % h = plot(tt, cumsum(erps_rwmd_nontarget_eog_rem(i, :)), 'b-.', 'linewidth', 4); H = cat(2, H, h); lgnd = cat(2, lgnd, 'RWMD Non-targets Cumsum (post EOG removal)');
        grid
        legend(H, lgnd);
        set(gca, 'fontsize', 16)
        if isempty(colnames)
            title(['Channel ', num2str(i)]);
        else
            title(['Channel ', num2str(i), ': ', colnames{i}]);
        end
    end
end

if 0 % for additional plots
    for i = 1:size(x, 1)
        figure
        hold on
        plot(tt, cumsum(squeeze(stacked_targets(i, :, :)), 2)', 'color', [1, 0.6, 0.6]);
        plot(tt, cumsum(squeeze(stacked_nontargets(i, :, :)), 2)', 'color', [0.6, 0.6, 1]);
        H = [];
        lgnd = {};
        % h = plot(tt, erps_mean_target(i, :), 'm:', 'linewidth', 2); H = cat(2, H, h); lgnd = cat(2, lgnd, 'Mean Targets');
        % h = plot(tt, erps_rwa_target(i, :), 'm--', 'linewidth', 2); H = cat(2, H, h); lgnd = cat(2, lgnd, 'RWA Targets');
        % h = plot(tt, erps_rwmd_target(i, :), 'm', 'linewidth', 2); H = cat(2, H, h); lgnd = cat(2, lgnd, 'RWMD Targets');
        h = plot(tt, cumsum(erps_rwmd_target(i, :)), 'm.-', 'linewidth', 4); H = cat(2, H, h); lgnd = cat(2, lgnd, 'RWMD Targets Cumsum');
        % h = plot(tt, erps_mean_nontarget(i, :), 'c:', 'linewidth', 2); H = cat(2, H, h); lgnd = cat(2, lgnd, 'Mean Non-targets');
        % h = plot(tt, erps_rwa_nontarget(i, :), 'c--', 'linewidth', 2); H = cat(2, H, h); lgnd = cat(2, lgnd, 'RWA Non-targets');
        % h = plot(tt, erps_rwmd_nontarget(i, :), 'c', 'linewidth', 2); H = cat(2, H, h); lgnd = cat(2, lgnd, 'RWMD Non-targets');
        h = plot(tt, cumsum(erps_rwmd_nontarget(i, :)), 'c.-', 'linewidth', 4); H = cat(2, H, h); lgnd = cat(2, lgnd, 'RWMD Non-targets Cumsum');
        % h = plot(tt, erps_mean_target_eog_rem(i, :), 'r:', 'linewidth', 2); H = cat(2, H, h); lgnd = cat(2, lgnd, 'Mean Targets (post EOG removal)');
        % h = plot(tt, erps_rwa_target_eog_rem(i, :), 'r--', 'linewidth', 2); H = cat(2, H, h); lgnd = cat(2, lgnd, 'RWA Targets (post EOG removal)');
        % h = plot(tt, erps_rwmd_target_eog_rem(i, :), 'r', 'linewidth', 2); H = cat(2, H, h); lgnd = cat(2, lgnd, 'RWMD Targets (post EOG removal)');
        h = plot(tt, cumsum(erps_rwmd_target_eog_rem(i, :)), 'r.-', 'linewidth', 4); H = cat(2, H, h); lgnd = cat(2, lgnd, 'RWMD Targets Cumsum (post EOG removal)');
        % h = plot(tt, erps_mean_nontarget_eog_rem(i, :), 'b:', 'linewidth', 2); H = cat(2, H, h); lgnd = cat(2, lgnd, 'Mean Non-targets (post EOG removal)');
        % h = plot(tt, erps_rwa_nontarget_eog_rem(i, :), 'b--', 'linewidth', 2); H = cat(2, H, h); lgnd = cat(2, lgnd, 'RWA Non-targets (post EOG removal)');
        % h = plot(tt, erps_rwmd_nontarget_eog_rem(i, :), 'b', 'linewidth', 2); H = cat(2, H, h); lgnd = cat(2, lgnd, 'RWMD Non-targets (post EOG removal)');
        h = plot(tt, cumsum(erps_rwmd_nontarget_eog_rem(i, :)), 'b.-', 'linewidth', 4); H = cat(2, H, h); lgnd = cat(2, lgnd, 'RWMD Non-targets Cumsum (post EOG removal)');
        grid
        legend(H, lgnd);
        set(gca, 'fontsize', 16)
        if isempty(colnames)
            title(['Channel ', num2str(i)]);
        else
            title(['Channel ', num2str(i), ': ', colnames{i}]);
        end
    end
end

figure
hold on
plot(tt, squeeze(mean(mean(cumsum(stacked_targets, 3), 2), 1)), 'color', 'r');
plot(tt, squeeze(mean(mean(cumsum(stacked_nontargets, 3), 2), 1)), 'color', 'b');
legend('Targets', 'Non-Targets');
grid
xlabel('time(s)');
set(gca, 'fontsize', 16)
title('Cumsum of the averege stackec-events ')