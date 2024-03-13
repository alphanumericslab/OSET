function [peak_indexes, hr, hr_smoothed, samples_sat, samples_sat_bp, samples_sat_bp_env, amp_inst, f_inst, f_inst_smoothed, f_hilbert, amp_hilbert, bw_inst] = peak_det_acoustic(samples, fs, varargin)
% PEAK_DET_ACOUSTIC Detects peaks in acoustic cardiac signals for heart rate analysis.
%
%   This function is designed for the analysis of acoustic cardiac signals, such as
%   seismocardiograms (SCG), phonocardiograms (PCG), and one-dimensional Doppler signals.
%   It detects signal peaks, calculates heart rate (HR), applies signal smoothing, and
%   performs various signal processing steps including saturation, filtering, and
%   envelope calculation. The function is highly customizable through parameters, allowing
%   adjustments for specific applications.
%
%   Syntax:
%   [peak_indexes, hr, hr_smoothed, samples_sat, samples_sat_bp, samples_sat_bp_env, amp_inst, f_inst, f_inst_smoothed, f_hilbert, amp_hilbert, bw_inst] = peak_det_acoustic(samples, fs)
%   Uses internal default parameters for processing the acoustic signal samples with sampling frequency fs.
%
%   [peak_indexes, hr, hr_smoothed, samples_sat, samples_sat_bp, samples_sat_bp_env] = peak_det_acoustic(samples, fs, params)
%   Allows custom parameters to be specified via the 'params' structure for tailored signal
%   analysis, accommodating the unique characteristics of different acoustic cardiac signals.
%
%   Inputs:
%   samples - a vector of acoustic time-series to be analyzed.
%   fs - sampling frequency of the signal in Hz.
%   params - (Optional) a structure with custom parameters. If not provided, default values
%            are used. The structure may contain the following fields:
%       verbose (bool) - enable verbose mode displaying all default parameters. default: true.
%       plot_results (bool) - plot signal processing stages and results. default: true.
%       plot_inst_freq (bool) - plot instantaneous frequency analysis. default: true.
%       inst_params_wlen (float) - window length for instantaneous frequency, in seconds. default: 0.5.
%       inst_params_smoothing_wlen (float) - smoothing window length for instantaneous frequency, in seconds. default: 3.0.
%       inst_params_smoothing_lambda (float) - smoothing parameter lambda for instantaneous frequency (using piece-wise Tikhonov regularization). default: 10.0.
%       tanh_saturation_ksigma (float) - saturation level for signal amplitude, in k-times standard deviation of filtered signal envelope. default: 5.0.
%       plot_spectrogram (bool) - enable plotting of the signal spectrogram. Computationally demanding; advised only for short segments and design purposes. default: false.
%       spectrogram_wlen (float) - window length for spectrogram, in seconds. default: 0.15.
%       spectrogram_win_overlap (float) - overlap between windows in spectrogram, in seconds. default: 0.14.
%       spectrogram_nfft (int) - number of DFT points for spectrogram. default: 64.
%       filter_type (string) - Type of filter applied to the signal
%       ('WAVELET', 'MATCHED', 'CUSTOM', 'BYPASS'). default: 'WAVELET'. Note: Conisder using custom filters for improved performance
%       matched_filter_template (vector) - template for matched filtering (relevant if filter_type is 'MATCHED'). default: [].
%       h_bp_b (vector), h_bp_a (vector) - numerator and denominator coefficients for custom filter (relevant if filter_type is 'CUSTOM'). default: [].
%       env_avg_wlen1, env_avg_wlen2, env_avg_wlen3 (float) - window lengths for power envelope averaging, in seconds. defaults: 0.13, 0.23, 0.37.
%       peaks_min_dist (float) - minimum distance between peaks, in seconds. default: 0.6, equivalent to max 100 beats per minute (adjust per application)
%       peaks_prom_fraction (float) - fraction of median of power envelope's local peaks to set as the minimum peak prominence. default: 0.05.
%       hr_max_jumps (float) - maximum allowed jump in HR between consecutive beats for hr_smoothed, in bpm. default: 20.0.
%       hr_smoothing_window (int) - window length for HR smoothing. default: 25.
%       plot_time_unit (string) - unit of time for plotting ('seconds', 'minutes'). default: 'seconds'.
%       Note: Additional parameters for wavelet and other filters are also available.
%
%   Outputs:
%   peak_indexes - indices of detected peaks within the signal.
%   hr - calculated heart rate time-series based on detected peaks, in beats per minute
%   hr_smoothed - smoothed heart rate time-series, in beats per minute
%   samples_sat - signal after amplitude saturation.
%   samples_sat_bp - signal after bandpass or specified filtering.
%   samples_sat_bp_env - power envelope of the filtered signal.
%   amp_inst - instantaneous signal amplitude obtained from signal power over sliding windows.
%   f_inst - instantaneous frequency obtained from Fourier domain center frequency over sliding windows
%   f_inst_smoothed - smoothed version of f_inst using a piece-wise Tikhonov regularization filter
%   f_hilbert - instantaneous frequency obtained from the Hilbert transform
%   amp_hilbert - instantaneous amplitude obtained from the Hilbert transform
%   bw_inst - instantaneous bandwidth obtained from Fourier domain center frequency deviation over sliding windows
%
%   Example Usage:
%   % default parameters for PCG analysis:
%   [pidx, hr, hrs, sats, satbp, satbp_env] = peak_det_acoustic(pcgSignal, 1000);
%
%   % Custom parameters for SCG analysis:
%   % Design filter using filterDesigner or other tools
%   Fs = 100;  % Sampling Frequency
%   Fstop1 = 4.0;              % First Stopband Frequency
%   Fpass1 = 5.0;              % First Passband Frequency
%   Fpass2 = 7.0;              % Second Passband Frequency
%   Fstop2 = 8.0;              % Second Stopband Frequency
%   Dstop1 = 0.01;             % First Stopband Attenuation
%   Dpass  = 0.0057563991496;  % Passband Ripple
%   Dstop2 = 0.01;             % Second Stopband Attenuation
%   dens   = 20;               % Density Factor
% 
%   % Calculate the order from the parameters using FIRPMORD.
%   [N, Fo, Ao, W] = firpmord([Fstop1 Fpass1 Fpass2 Fstop2]/(Fs/2), [0 1 0], [Dstop1 Dpass Dstop2]);
% 
%   % Calculate the coefficients using the FIRPM function.
%   params.h_bp_b = firpm(N, Fo, Ao, W, {dens});
%   params.h_bp_a = 1;
%   params.filter_type = 'CUSTOM';
%   params.verbose = false;
%   [pidx, hr, hrs, sats, satbp, satbp_env] = peak_det_acoustic(scgSignal, 100, params);
% 
%   Notes:
%   - The function uses some functions from the OSET tool.
%   - Make sure to adjust parameters to match the signal characteristics for optimal detection and analysis.
%
%   Revision History:
%       2024: First release
%
%   Reza Sameni, 2024
%   The Open-Source Electrophysiological Toolbox
%   https://github.com/alphanumericslab/OSET

% check if custom params have been provided
if nargin > 2 && ~isempty(varargin{1})
    params = varargin{1};
else
    params = struct;
end

if ~isfield(params, 'verbose') || isempty(params.verbose)
    params.verbose = true;
end
if params.verbose, disp('displaying all default parameters in verbose mode; set params.verbose=false for silent mode'), end

if ~isfield(params, 'plot_results') || isempty(params.plot_results)
    params.plot_results = true;
    if params.verbose, disp(['params.plot_results = ', num2str(params.plot_results)]), end
end

if ~isfield(params, 'plot_inst_freq') || isempty(params.plot_inst_freq)
    params.plot_inst_freq = true;
    if params.verbose, disp(['params.plot_inst_freq = ', num2str(params.plot_inst_freq)]), end
end

if ~isfield(params, 'inst_params_wlen') || isempty(params.inst_params_wlen)
    params.inst_params_wlen = 0.03;
    if params.verbose, disp(['params.inst_params_wlen = ', num2str(params.inst_params_wlen)]), end
end
params.inst_params_wlen = max(1, round(params.inst_params_wlen*fs)); % convert into samples

if ~isfield(params, 'inst_params_smoothing_wlen') || isempty(params.inst_params_smoothing_wlen)
    params.inst_params_smoothing_wlen = 3.0;
    if params.verbose, disp(['params.inst_params_smoothing_wlen = ', num2str(params.inst_params_smoothing_wlen)]), end
end
params.inst_params_smoothing_wlen = round(params.inst_params_smoothing_wlen*fs); % convert into samples

if ~isfield(params, 'inst_params_smoothing_lambda') || isempty(params.inst_params_smoothing_lambda)
    params.inst_params_smoothing_lambda = 10.0;
    if params.verbose, disp(['params.inst_params_smoothing_lambda = ', num2str(params.inst_params_smoothing_lambda)]), end
end

if ~isfield(params, 'tanh_saturation_ksigma') || isempty(params.tanh_saturation_ksigma)
    params.tanh_saturation_ksigma = 5.0;
    if params.verbose, disp(['params.tanh_saturation_ksigma = ', num2str(params.tanh_saturation_ksigma)]), end
end

if ~isfield(params, 'plot_spectrogram') || isempty(params.plot_spectrogram)
    params.plot_spectrogram = false;
    if params.verbose, disp(['params.plot_spectrogram = ', num2str(params.plot_spectrogram)]), end
end

if ~isfield(params, 'spectrogram_wlen') || isempty(params.spectrogram_wlen)
    params.spectrogram_wlen = 0.15;
    if params.verbose, disp(['params.spectrogram_wlen = ', num2str(params.spectrogram_wlen)]), end
end
params.spectrogram_wlen = round(params.spectrogram_wlen*fs); % convert into samples

if ~isfield(params, 'spectrogram_win_overlap') || isempty(params.spectrogram_win_overlap)
    params.spectrogram_win_overlap = 0.14;
    if params.verbose, disp(['params.spectrogram_win_overlap = ', num2str(params.spectrogram_win_overlap)]), end
    params.spectrogram_win_overlap = round(params.spectrogram_win_overlap*fs); % convert into samples
end

if ~isfield(params, 'spectrogram_nfft') || isempty(params.spectrogram_nfft)
    params.spectrogram_nfft = 64;
    if params.verbose, disp(['params.spectrogram_nfft = ', num2str(params.spectrogram_nfft)]), end
end

if ~isfield(params, 'filter_type') || isempty(params.filter_type)
    params.filter_type = 'WAVELET';
    if params.verbose, disp(['params.filter_type = ', params.filter_type]), end
end

if ~isfield(params, 'matched_filter_template') || isempty(params.matched_filter_template)
    params.matched_filter_template = [];
    if params.verbose, disp(['params.matched_filter_template = ', num2str(params.matched_filter_template)]), end
end

if ~isfield(params, 'h_bp_b') || isempty(params.h_bp_b)
    params.h_bp_b = [];
    if params.verbose, disp(['params.h_bp_b = ', num2str(params.h_bp_b)]), end
end

if ~isfield(params, 'h_bp_a') || isempty(params.h_bp_a)
    params.h_bp_a = [];
    if params.verbose, disp(['params.h_bp_a = ', num2str(params.h_bp_a)]), end
end

if ~isfield(params, 'env_avg_wlen1') || isempty(params.env_avg_wlen1)
    params.env_avg_wlen1 = 0.13;
    if params.verbose, disp(['params.env_avg_wlen1 = ', num2str(params.env_avg_wlen1)]), end
end
params.env_avg_wlen1 = round(fs * params.env_avg_wlen1); % convert into samples

if ~isfield(params, 'env_avg_wlen2') || isempty(params.env_avg_wlen2)
    params.env_avg_wlen2 = 0.23;
    if params.verbose, disp(['params.env_avg_wlen2 = ', num2str(params.env_avg_wlen2)]), end
end
params.env_avg_wlen2 = round(fs * params.env_avg_wlen2); % convert into samples

if ~isfield(params, 'env_avg_wlen3') || isempty(params.env_avg_wlen3)
    params.env_avg_wlen3 = 0.37;
    if params.verbose, disp(['params.env_avg_wlen3 = ', num2str(params.env_avg_wlen3)]), end
end
params.env_avg_wlen3 = round(fs * params.env_avg_wlen3); % convert into samples

if ~isfield(params, 'peaks_min_dist') || isempty(params.peaks_min_dist)
    params.peaks_min_dist = 0.6;
    if params.verbose, disp(['params.peaks_min_dist = ', num2str(params.peaks_min_dist)]), end
end
params.peaks_min_dist = round(params.peaks_min_dist * fs);

if ~isfield(params, 'peaks_prom_fraction') || isempty(params.peaks_prom_fraction)
    params.peaks_prom_fraction = 0.05;
    if params.verbose, disp(['params.peaks_prom_fraction = ', num2str(params.peaks_prom_fraction)]), end
end

if ~isfield(params, 'hr_max_jumps') || isempty(params.hr_max_jumps)
    params.hr_max_jumps = 20.0;
    if params.verbose, disp(['params.hr_max_jumps = ', num2str(params.hr_max_jumps)]), end
end

if ~isfield(params, 'hr_smoothing_window') || isempty(params.hr_smoothing_window)
    params.hr_smoothing_window = 25;
    if params.verbose, disp(['params.hr_smoothing_window = ', num2str(params.hr_smoothing_window)]), end
end

if ~isfield(params, 'wden_type') || isempty(params.wden_type)
    params.wden_type = 'sym4'; % mother wavelet used for wavelet denoising
    if params.verbose, disp(['params.wden_type = ', params.wden_type]), end
end
if ~isfield(params, 'wden_upper_level') || isempty(params.wden_upper_level)
    f_high = min(7.0, fs/2); % upper cutoff frequency
    params.wden_upper_level = round(log2(fs/f_high)); % Upper frequency range for the wavelet denoiser
    if params.verbose, disp(['params.wden_upper_level = ', num2str(params.wden_upper_level)]), end
end
if ~isfield(params, 'wden_lower_level') || isempty(params.wden_lower_level)
    f_low = 3.0; % lower cutoff frequency
    params.wden_lower_level = round(log2(fs/f_low)); % Lower frequency range for the wavelet denoiser
    if params.verbose, disp(['params.wden_lower_level = ', num2str(params.wden_lower_level)]), end
end
if params.wden_upper_level > params.wden_lower_level
    error('requested wavelet band too narrow, or sampling frequency is too low.')
end

if ~isfield(params, 'plot_time_unit') || isempty(params.plot_time_unit)
    params.plot_time_unit = 'seconds';
    if params.verbose, disp(['params.plot_time_unit = ', params.plot_time_unit]), end
end

switch params.plot_time_unit
    case 'seconds'
        time = (0:size(samples, 2)-1)/fs;
        time_label = 'time(s)';
    case 'minutes'
        time = (0:size(samples, 2)-1)/fs/60;
        time_label = 'time(min)';
    otherwise
        error('undefined plot time unit')
end

% make sure the signal is zero-mean and a row vector
samples = samples(:)' - mean(samples);

% plot signal's instantaneous frequency and amplitude (useful for designing custom filters)
if params.plot_inst_freq
    noncausal = 1;
    [amp_inst, f_inst, bw_inst, ~, ~, ~, f_hilbert, amp_hilbert] = instantaneous_signal_params(samples, fs, params.inst_params_wlen, ones(1, params.inst_params_wlen), 0, noncausal);

    % smooth the instantaneous frequency
    f_inst_smoothed = ecg_den_seg_wise_smoother(f_inst, 2, params.inst_params_smoothing_wlen, params.inst_params_smoothing_lambda);

    figure
    subplot(211)
    plot(time, samples);
    hold on
    plot(time, amp_hilbert)
    plot(time, amp_inst)
    xlabel(time_label);
    ylabel('Amplitude');
    legend('Signal', 'Hilbert envelope', 'Power envelope');
    grid
    set(gca, 'fontsize', 16)
    subplot(212)
    plot(time, f_hilbert);
    hold on
    plot(time, f_inst);
    plot(time, f_inst_smoothed)
    plot(time, bw_inst)
    xlabel(time_label);
    ylabel('Frequency(Hz)');
    legend('Inst freq Hilbert', 'Inst freq', 'Inst freq smoothed', 'Inst bandwidth');
    grid
    set(gca, 'fontsize', 16)
    sgtitle('Signal vs Hilbert amplitude and instantaneous frequencies using density-based and Hilbert method');
else
    amp_inst = [];
    f_inst = [];
    f_hilbert = [];
    bw_inst = [];
    f_inst_smoothed = [];
end

if params.plot_spectrogram
    figure
    spectrogram(samples, hamming(params.spectrogram_wlen), params.spectrogram_win_overlap, params.spectrogram_nfft, fs, 'yaxis');
    title('Spectrogram');
    set(gca, 'fontsize', 16)
end

% saturate the signal amplitude outliers at k_sigma
samples_sat = tanh_saturation(samples, params.tanh_saturation_ksigma, 'ksigma');

% narrow band filtering of the signal to enhance the local peaks
switch params.filter_type
    case 'BYPASS'
        samples_sat_bp = samples_sat;
    case 'WAVELET' % wavelet denoiser
        samples_sat_bp = zeros(size(samples_sat));
        for kk = 1 : size(samples_sat, 1)
            wt = modwt(samples_sat(kk, :), params.wden_lower_level);
            wtrec = zeros(size(wt));
            wtrec(params.wden_upper_level : params.wden_lower_level, :) = wt(params.wden_upper_level : params.wden_lower_level,:);
            samples_sat_bp(kk, :) = imodwt(wtrec, params.wden_type);
        end
    case 'MATCHED' % matched filter
        samples_sat_match_filt = filter(params.matched_filter_template, 1, samples_sat);
        samples_sat_match_filt = std(samples_sat)*samples_sat_match_filt./std(samples_sat_match_filt);
        samples_sat_match_filt = samples_sat_match_filt(round(length(filter_template)/2) : end);
        samples_sat_match_filt = cat(2, samples_sat_match_filt, zeros(1, length(samples_sat) - length(samples_sat_match_filt)));
        samples_sat_bp = samples_sat_match_filt;
    case 'CUSTOM' % any custom filter
        samples_sat_bp = filtfilt(params.h_bp_b, params.h_bp_a, samples_sat);
    otherwise
        error('Unknown preprocessing filter');
end

% calculate power envelope of the filtered signal over multiple sliding windows
samples_sat_bp_env1 = sqrt(filtfilt(ones(1, params.env_avg_wlen1), params.env_avg_wlen1, samples_sat_bp.^2));
samples_sat_bp_env2 = sqrt(filtfilt(ones(1, params.env_avg_wlen2), params.env_avg_wlen2, samples_sat_bp.^2));
samples_sat_bp_env3 = sqrt(filtfilt(ones(1, params.env_avg_wlen3), params.env_avg_wlen3, samples_sat_bp.^2));

% multiply the envelopes to enhance common local peaks
samples_sat_bp_env = (samples_sat_bp_env1 .* samples_sat_bp_env2 .* samples_sat_bp_env3) .^(1/3);

% run peak detection regardless of local peak prominence
[~, peak_indexes_all_prominence] = findpeaks(samples_sat_bp_env, 'MinPeakDistance', params.peaks_min_dist);

% set peak prominnece relative to the median of all detected beats
peaks_priminence = params.peaks_prom_fraction * median(samples_sat_bp_env(peak_indexes_all_prominence));

% rerun peak detection using the selected prominence
[~, peak_indexes] = findpeaks(samples_sat_bp_env, 'MinPeakDistance', params.peaks_min_dist, 'MinPeakProminence', peaks_priminence);

% calculate the heart rate time-series
hr = 60*fs./diff(peak_indexes); hr = cat(2, hr(1), hr);

% smooth the heart rate time-series using conditional median filtering
hr_smoothed = HRSmoother2(hr, params.hr_max_jumps, params.hr_smoothing_window, 'normal');

% plot results
if params.plot_results
    figure
    lgnd = {};
    hold on
    plot(time, samples); lgnd = cat(1, lgnd, 'samples');
    plot(time, samples_sat); lgnd = cat(1, lgnd, 'samples_sat');
    plot(time, samples_sat_bp); lgnd = cat(1, lgnd, 'samples_sat_bp');
    plot(time, samples_sat_bp_env); lgnd = cat(1, lgnd, 'samples_sat_bp_env');
    plot(time(peak_indexes), samples_sat_bp_env(peak_indexes), 'o', 'markersize', 18); lgnd = cat(1, lgnd, 'peak_indexes');
    legend(lgnd, 'interpreter', 'none');
    grid
    xlabel(time_label);
    ylabel('Amplitude');
    title('Input vs processed signals');
    set(gca, 'fontsize', 16);

    figure
    hold on
    plot(time(peak_indexes), hr);
    plot(time(peak_indexes), hr_smoothed);
    grid
    xlabel(time_label);
    ylabel('Heart rate (bpm)');
    legend('Raw HR', 'Smoothed HR');
    title('Raw and smoothed heart rate time-series');
    set(gca, 'fontsize', 16);
end