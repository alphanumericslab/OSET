function [peaks, peak_indexes, peak_indexes_consensus, qrs_likelihood] = peak_det_probabilistic(data, fs, varargin)
% function [peaks, peak_indexes, qrs_likelihood] = peak_det_probabilistic(data, fs)
% A probabilistic R-peak detector based on local peaks sorting
%
% Note: Under development
%
% Usage:
%   [peaks, peak_indexes, peak_indexes_consensus, qrs_likelihood] = peak_det_probabilistic(data, fs)
%   [peaks, peak_indexes, peak_indexes_consensus, qrs_likelihood] = peak_det_probabilistic(data, fs, params)
%
% Inputs:
%   data: single or multichannel ECG signal with row-wise channels
%   fs: sampling frequency
%   params: a structure containing the algorithm parameters as follows (all
%       parameters have default values if not provided as input):
%       params.verbose: true (all default settings will be displayed as the
%           function is called); or false (silent mode)
%       params.PLOT_DIAGNOSTIC: true (all intermediate peak refinement
%           figures will be plotted); or false (no intermediate figures are generated)
%       params.saturate: true (saturate the channels before R-peak detection) or false (no saturation). Default is true
%       params.filter_type: preprocessing filter type. Options: 'MDMN',
%           'BANDPASS_FILTER', 'MD3MN', 'GAUSSIAN_MATCHED_FILTER',
%           'MULT_MATCHED_FILTER_ENV', or 'WAVELET'. Default is 'WAVELET'
%       Note: each mode has predefined parameters that 

%       params.bp_lower_cutoff: lower cutoff frequency of R-peak detector bandpass filter (in Hz), when params.filter_type = 'BANDPASS_FILTER'. Default = 1
%       params.bp_upper_cutoff: upper cutoff frequency of R-peak detector bandpass filter (in Hz), when params.filter_type = 'BANDPASS_FILTER'. Default = 40
%       params.gauss_match_filt_span: the preprocessing gaussian shape matched filter span, when params.filter_type = 'MATCHED_FILTER'. Default = 0.2
%       params.gaus_match_filt_sigma: the preprocessing gaussian shape matched filter STD, when params.filter_type = 'MATCHED_FILTER'. Default = 0.01
%       params.power_env_wlen: power envelope moving average window length (in seconds). Default = 0.03
%       params.n_bins: number of local peaks histogram bins. Default = max(250, min(500, 10% of signal length))
%       params.power_env_hist_peak_th: the top percentile of the signal's power envelope, considered for R-peak detection. Default = 0.9 (top 10%)
%       params.min_peak_distance: the R-peak detector search window length (in seconds). Default = 0.2
%       params.likelihood_sigma: R-peak likelihood estimator STD in seconds (assuming a Gaussian that drops from the peak maximum). Default = 0.01
%       params.max_likelihood_span: R-peak likelihood estimator max window span in seconds (assuming a Gaussian that drops from the peak maximum). Default = 0.1
%       params.RemoveUnsimilarBeats: remove the R-peaks that do not have similar morphology (true/false). Default = true
%
% Outputs:
%   peaks: a vector with the signal length with 1's at the estimated R-peaks
%   peak_indexes: the estimated R-peak indexes
%   qrs_likelihood: the R-peak likelihood vector (with maximums at the
%       estimated R-peaks useful for classification and scoring purposes)
%
% Revision History:
%   2021: First release
%   2023: Replaced deprecated function PeakDetectionProbabilistic
%   2023: Expanded excess peak rejection functions
%
% Reza Sameni, 2009-2023
% The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET


% use default values when no parameters are set
if nargin > 2
    params = varargin{1};
else
    params = [];
end

if ~isfield(params, 'verbose')
    params.verbose = false;
end
if params.verbose
    disp('Operating in verbose mode. Default settings will be displayed.')
end

if ~isfield(params, 'PLOT_DIAGNOSTIC')
    params.PLOT_DIAGNOSTIC = false;
end
if params.PLOT_DIAGNOSTIC
    disp('Operating in diagnostic plot mode. All peak refinement figures will be plotted.')
end

sig_len = size(data, 2); % signal length

%% pass the channels through a narrow bandpass filter or a matched filter to remove the baseline and the T-waves
if ~isfield(params, 'filter_type')
    params.filter_type = 'WAVELET';
    if params.verbose, disp(['params.filter_type = ', params.filter_type]), end
end
switch params.filter_type
    case 'MDMN' % Median + MA filter
        if ~isfield(params, 'wlen_md')
            params.wlen_md = 0.075;
            if params.verbose, disp(['params.wlen_md = ', num2str(params.wlen_md)]), end
        end
        if ~isfield(params, 'wlen_mn')
            params.wlen_mn = 0.05;
            if params.verbose, disp(['params.wlen_mn = ', num2str(params.wlen_mn)]), end
        end
        data_filtered = zeros(size(data));
        for kk = 1 : size(data, 1)
            bl1 = baseline_sliding_window(data(kk, :), round(params.wlen_md * fs), 'md');
            bl2 = baseline_sliding_window(bl1, round(params.wlen_mn * fs), 'mn');
            data_filtered(kk, :) = data(kk, :) - bl2;
            % figure
            % plot(data(kk, :));
            % hold on
            % plot(bl1);
            % plot(bl2);
            % plot(data_filtered(kk, :));
            % grid
            % legend({'data', 'bl1', 'bl2', 'data_filtered'}, 'interpreter', 'none');
        end
    case 'BANDPASS_FILTER' % Bandpass filter
        if ~isfield(params, 'bp_lower_cutoff')
            params.bp_lower_cutoff = 8.0;
            if params.verbose, disp(['params.bp_lower_cutoff = ', num2str(params.bp_lower_cutoff)]), end
        end
        if ~isfield(params, 'bp_upper_cutoff')
            params.bp_upper_cutoff = 40.0;
            if params.verbose, disp(['params.bp_upper_cutoff = ', num2str(params.bp_upper_cutoff)]), end
        end
        data_lp = lp_filter_zero_phase(data, params.bp_upper_cutoff/fs);
        data_filtered = data_lp - lp_filter_zero_phase(data_lp, params.bp_lower_cutoff/fs);
    case 'MD3MN' % Three median filters followed by a moving average
        if ~isfield(params, 'wlen_md1')
            params.wlen_md1 = 0.07;
            if params.verbose, disp(['params.wlen_md1 = ', num2str(params.wlen_md1)]), end
        end
        if ~isfield(params, 'wlen_md2')
            params.wlen_md2 = 0.08;
            if params.verbose, disp(['params.wlen_md2 = ', num2str(params.wlen_md2)]), end
        end
        if ~isfield(params, 'wlen_md3')
            params.wlen_md3 = 0.09;
            if params.verbose, disp(['params.wlen_md3 = ', num2str(params.wlen_md3)]), end
        end
        if ~isfield(params, 'wlen_mn')
            params.wlen_mn = 0.01;
            if params.verbose, disp(['params.wlen_mn = ', num2str(params.wlen_mn)]), end
        end
        data_filtered = zeros(size(data));
        for kk = 1 : size(data, 1)
            bl1 = baseline_sliding_window(data(kk, :), round(params.wlen_md1 * fs), 'md');
            bl2 = baseline_sliding_window(data(kk, :), round(params.wlen_md2 * fs), 'md');
            bl3 = baseline_sliding_window(data(kk, :), round(params.wlen_md3 * fs), 'md');
            baseline = baseline_sliding_window((bl1 + bl2 + bl3) / 3, round(params.wlen_mn * fs), 'mn');
            data_filtered(kk, :) = data(kk, :) - baseline;
        end

    case 'GAUSSIAN_MATCHED_FILTER' % Matched filter
        if ~isfield(params, 'gauss_match_filt_span')
            params.gauss_match_filt_span = 0.1;
            if params.verbose, disp(['params.gauss_match_filt_span = ', num2str(params.gauss_match_filt_span)]), end
        end
        if ~isfield(params, 'gaus_match_filt_sigma')
            params.gaus_match_filt_sigma = 0.01;
            if params.verbose, disp(['params.gaus_match_filt_sigma = ', num2str(params.gaus_match_filt_sigma)]), end
        end
        t_matched = -params.gauss_match_filt_span/2 : 1 / fs : params.gauss_match_filt_span/2;
        match_filt_template = exp(-t_matched.^2/(2 * params.gaus_match_filt_sigma ^ 2));
        lag_match = round(length(match_filt_template)/2);
        data_filtered = zeros(size(data));
        for kk = 1 : size(data, 1)
            band_passsed = conv(match_filt_template, data(kk, :))/sum(match_filt_template.^2);
            data_filtered(kk, :) = band_passsed(lag_match : sig_len + lag_match - 1);
        end
    case 'MULT_MATCHED_FILTER_ENV' % multi-template matched filter
        if ~isfield(params, 'gauss_match_filt_span')
            params.gauss_match_filt_span = 0.1;
            if params.verbose, disp(['params.gauss_match_filt_span = ', num2str(params.gauss_match_filt_span)]), end
        end
        if ~isfield(params, 'gaus_match_filt_sigma')
            params.gaus_match_filt_sigma = 0.005;
            if params.verbose, disp(['params.gaus_match_filt_sigma = ', num2str(params.gaus_match_filt_sigma)]), end
        end
        t_matched = -params.gauss_match_filt_span/2 : 1 / fs : params.gauss_match_filt_span/2;
        match_filt_template1 = exp(-t_matched.^2/(2 * params.gaus_match_filt_sigma ^ 2));
        match_filt_template2 = diff(match_filt_template1);
        lag_match1 = round(length(match_filt_template1)/2);
        lag_match2 = round(length(match_filt_template2)/2);
        data_filtered = zeros(size(data));
        for kk = 1 : size(data, 1)
            band_passsed1 = conv(match_filt_template1, data(kk, :))/sum(match_filt_template1.^2);
            band_passsed1 = band_passsed1(lag_match1 : sig_len + lag_match1 - 1);

            band_passsed2 = conv(match_filt_template2, data(kk, :))/sum(match_filt_template2.^2);
            band_passsed2 = band_passsed2(lag_match2 : sig_len + lag_match2 - 1);

            data_filtered(kk, :) = sqrt(band_passsed1.^2 + band_passsed2.^2);
        end
    case 'WAVELET' % Wavelet denoiser
        if ~isfield(params, 'wden_type')
            params.wden_type = 'sym4'; % mother wavelet used for wavelet denoising
            if params.verbose, disp(['params.wden_type = ', params.wden_type]), end
        end
        if ~isfield(params, 'wden_upper_level')
            params.wden_upper_level = floor(log2(fs/80.0)); % Upper frequency range for the wavelet denoiser
            if params.verbose, disp(['params.wden_upper_level = ', num2str(params.wden_upper_level)]), end
        end
        if ~isfield(params, 'wden_lower_level')
            params.wden_lower_level = ceil(log2(fs/25.0)); % Lower frequency range for the wavelet denoiser
            if params.verbose, disp(['params.wden_lower_level = ', num2str(params.wden_lower_level)]), end
        end

        data_filtered = zeros(size(data));
        for kk = 1 : size(data, 1)
            wt = modwt(data(kk, :), params.wden_lower_level);
            wtrec = zeros(size(wt));
            wtrec(params.wden_upper_level : params.wden_lower_level,:) = wt(params.wden_upper_level : params.wden_lower_level,:);
            data_filtered(kk, :) = imodwt(wtrec, params.wden_type);
            % figure
            % plot(data(kk, :));
            % hold on
            % plot(data(kk, :) - data_filtered(kk, :));
            % plot(data_filtered(kk, :));
            % grid
            % legend({'data', 'baseline', 'data_filtered'}, 'interpreter', 'none');
        end
    otherwise
        error('Unknown preprocessing filter');
end

%% saturate the channels at k_sigma times the channel-wise STD
if ~isfield(params, 'saturate')
    params.saturate = 1;
    if params.verbose, disp(['params.saturate = ', num2str(params.saturate)]), end
end
if isequal(params.saturate, 1)
    if ~isfield(params, 'sat_k_sigma')
        params.sat_k_sigma = 8.0;
        if params.verbose, disp(['params.sat_k_sigma = ', num2str(params.sat_k_sigma)]), end
    end
    data_filtered = tanh_saturation(data_filtered, params.sat_k_sigma, 'ksigma');
end

data_filtered_mn_all_channels = mean(data_filtered, 1);

%% calculate the power envelope of one or all channels (stage 1)
if ~isfield(params, 'power_env_wlen')
    params.power_env_wlen = 0.025;
    if params.verbose, disp(['params.power_env_wlen = ', num2str(params.power_env_wlen)]), end
end
power_env_wlen = ceil(params.power_env_wlen * fs);
% data_filtered_env1 = filtfilt(ones(power_env_wlen, 1), power_env_wlen, sqrt(mean(data_filtered.^2, 1)));
padding_len = round(1.0*fs);
if params.verbose, disp(['padding_len = ', num2str(padding_len)]), end
% data_filtered_padded = [repmat(data_filtered(:, 1), 1, padding_len), data_filtered, repmat(data_filtered(:, end), 1, padding_len)];
data_filtered_padded = [zeros(size(data_filtered, 1), padding_len), data_filtered, zeros(size(data_filtered, 1), padding_len)];
data_filtered_env1 = filtfilt(ones(power_env_wlen, 1), power_env_wlen, sqrt(mean(data_filtered_padded.^2, 1)));
data_filtered_env1 = data_filtered_env1(padding_len + 1 : padding_len + size(data_filtered, 2));

%% calculate the power envelope of one or all channels (stateg 2)
if ~isfield(params, 'two_stage_env')
    params.two_stage_env = true;
    if params.verbose, disp(['params.two_stage_env = ', num2str(params.two_stage_env)]), end
end
if isequal(params.two_stage_env, 1)
    if ~isfield(params, 'power_env_wlen2')
        params.power_env_wlen2 = 0.075;
        if params.verbose, disp(['   params.power_env_wlen2 = ', num2str(params.power_env_wlen2)]), end
        power_env_wlen2 = ceil(params.power_env_wlen2 * fs);
        % data_filtered_padded = [repmat(data_filtered(:, 1), 1, padding_len), data_filtered, repmat(data_filtered(:, end), 1, padding_len)];
        data_filtered_padded = [zeros(size(data_filtered, 1), padding_len), data_filtered, zeros(size(data_filtered, 1), padding_len)];
        data_filtered_env2 = filtfilt(ones(power_env_wlen2, 1), power_env_wlen2, sqrt(mean(data_filtered_padded.^2, 1)));
        data_filtered_env2 = data_filtered_env2(padding_len + 1 : padding_len + size(data_filtered, 2));
    end
    data_filtered_env = sqrt(data_filtered_env1 .* data_filtered_env2);
    % figure
    % plot(data_filtered')
    % hold on
    % plot(data_filtered_env1)
    % plot(data_filtered_env2)
    % plot(data_filtered_env)
    % grid
    % legend({'data_filtered', 'data_filtered_env1', 'data_filtered_env2', 'data_filtered_env'}, 'interpreter', 'none')
    % set(gca, 'fontsize', 18)
else
    data_filtered_env = data_filtered_env1;
end

%% pick the top percentage of the signal's power envelope for R-peak search
if ~isfield(params, 'power_env_hist_peak_th')
    params.power_env_hist_peak_th = 85.0;
    if params.verbose, disp(['params.power_env_hist_peak_th = ', num2str(params.power_env_hist_peak_th)]), end
end
bumps_indexes = refine_peaks_low_amp_peaks_prctile(data_filtered_env, 1:sig_len, params.power_env_hist_peak_th, params.PLOT_DIAGNOSTIC);

%% search for all local peaks within a given sliding window length
if ~isfield(params, 'min_peak_distance')
    params.min_peak_distance = 0.18;
    if params.verbose, disp(['params.min_peak_distance = ', num2str(params.min_peak_distance)]), end
end
rpeak_search_half_wlen = floor(fs * params.min_peak_distance);
env_pk_detect_mode = 'POS';
if params.verbose, disp(['env_pk_detect_mode = ', env_pk_detect_mode]), end
peak_indexes = refine_peaks_too_close_low_amp(data_filtered_env, bumps_indexes, rpeak_search_half_wlen, env_pk_detect_mode, params.PLOT_DIAGNOSTIC);

% first or last samples are not selected as peaks
if peak_indexes(1) == 1
    peak_indexes = peak_indexes(2:end);
end

if peak_indexes(end) == sig_len
    peak_indexes = peak_indexes(1:end - 1);
end

%% calculate a likelihood function for the R-peaks (useful for classification and scoring purposes)
if ~isfield(params, 'likelihood_sigma')
    %     params.likelihood_sigma = 0.01;
    params.likelihood_sigma = 0.3;
    if params.verbose, disp(['params.likelihood_sigma = ', num2str(params.likelihood_sigma)]), end
end
if ~isfield(params, 'max_likelihood_span')
    params.max_likelihood_span = 0.3;
    if params.verbose, disp(['params.max_likelihood_span = ', num2str(params.max_likelihood_span)]), end
end

% matched filter using average beat shape
if ~isfield(params, 'ENHANCE_MATCHED_FILTER')
    params.ENHANCE_MATCHED_FILTER = false;
    if params.verbose, disp(['params.ENHANCE_MATCHED_FILTER = ', num2str(params.ENHANCE_MATCHED_FILTER)]), end
end
if params.ENHANCE_MATCHED_FILTER
    [data_filtered_enhanced, data_filtered_enhanced_env] = signal_specific_matched_filter(data_filtered, peak_indexes);
    % figure
    % plot(data_filtered);
    % hold on
    % plot(data_filtered_enhanced);
    % plot(data_filtered_enhanced_env);
    % grid
    data_filtered = data_filtered_enhanced;
    data_filtered_env = data_filtered_enhanced_env;
end

%% Refine the extracted R-peaks
peaks = zeros(1, sig_len);
peaks(peak_indexes) = 1;
refinement_methods = {'PRE_REFINEMENT'};

% Refine the extracted beats
if ~isfield(params, 'REFINE_PEAKS')
    params.REFINE_PEAKS = true;
    if params.verbose, disp(['params.REFINE_PEAKS = ', num2str(params.REFINE_PEAKS)]), end
end
if params.REFINE_PEAKS

    if ~isfield(params, 'OMIT_LOW_AMP_PEAKS_PRCTL_FRAC')
        params.OMIT_LOW_AMP_PEAKS_PRCTL_FRAC = true;
        if params.verbose, disp(['params.OMIT_LOW_AMP_PEAKS_PRCTL_FRAC = ', num2str(params.OMIT_LOW_AMP_PEAKS_PRCTL_FRAC)]), end
    end
    % remove excess beats with extreme amplitude deviations from other peaks
    if params.OMIT_LOW_AMP_PEAKS_PRCTL_FRAC
        pparams.percentile = 90.0;
        pparams.percentile_fraction = 0.4;
        if params.verbose, disp(['   pparams.percentile = ', num2str(pparams.percentile)]), end
        if params.verbose, disp(['   pparams.percentile_fraction = ', num2str(pparams.percentile_fraction)]), end
        [~, peaks_OMIT_LOW_AMP_PEAKS_PRCTL_FRAC] = refine_peaks_low_amp_peaks_prctile_fraction(data_filtered_env, peak_indexes, pparams, params.PLOT_DIAGNOSTIC);
        peaks = cat(1, peaks, peaks_OMIT_LOW_AMP_PEAKS_PRCTL_FRAC);
        refinement_methods = cat(2, refinement_methods, 'OMIT_LOW_AMP_PEAKS_PRCTL_FRAC');
    end

    % omit beats with low variance
    if ~isfield(params, 'OMIT_LOW_POWER_BEATS')
        params.OMIT_LOW_POWER_BEATS = true;
        % max_amp_k_sigma = 2.0;
        pparams.beat_std_med_frac_th = 0.5;
        pparams.max_amp_prctile = 90.0;
        if params.verbose, disp(['params.OMIT_LOW_POWER_BEATS = ', num2str(params.OMIT_LOW_POWER_BEATS)]), end
        if params.verbose, disp(['   pparams.beat_std_med_frac_th = ', num2str(pparams.beat_std_med_frac_th)]), end
        if params.verbose, disp(['   pparams.max_amp_prctile = ', num2str(pparams.max_amp_prctile)]), end
    end
    if params.OMIT_LOW_POWER_BEATS
        % peak_indexes_refined = refine_peaks_low_power_beats(data_filtered_env, peak_indexes, max_amp_k_sigma, beat_std_med_frac_th, params.PLOT_DIAGNOSTIC);
        [~, peaks_OMIT_LOW_POWER_BEATS] = refine_peaks_low_power_beats(data_filtered_mn_all_channels, peak_indexes, pparams.max_amp_prctile, pparams.beat_std_med_frac_th, params.PLOT_DIAGNOSTIC);
        peaks = cat(1, peaks, peaks_OMIT_LOW_POWER_BEATS);
        refinement_methods = cat(2, refinement_methods, 'OMIT_LOW_POWER_BEATS');
    end

    % negative correlations
    if ~isfield(params, 'OMIT_NEG_CORRCOEF_BEATS')
        params.OMIT_NEG_CORRCOEF_BEATS = true;
        if params.verbose, disp(['params.OMIT_NEG_CORRCOEF_BEATS = ', num2str(params.OMIT_NEG_CORRCOEF_BEATS)]), end
    end
    if params.OMIT_NEG_CORRCOEF_BEATS
        [~, peaks_OMIT_NEG_CORRCOEF_BEATS] = refine_peaks_waveform_similarity(data_filtered_mn_all_channels, peak_indexes, [], 'NEG-CORR', params.PLOT_DIAGNOSTIC);
        peaks = cat(1, peaks, peaks_OMIT_NEG_CORRCOEF_BEATS);
        refinement_methods = cat(2, refinement_methods, 'OMIT_NEG_CORRCOEF_BEATS');
    end

    % omit beats with low correlations
    if ~isfield(params, 'OMIT_LOW_CORRCOEF_BEATS')
        params.OMIT_LOW_CORRCOEF_BEATS = true;
        pparams.k_sigma = 3.0;
        if params.verbose, disp(['params.OMIT_LOW_CORRCOEF_BEATS = ', num2str(params.OMIT_LOW_CORRCOEF_BEATS)]), end
        if params.verbose, disp(['   pparams.k_sigma = ', num2str(pparams.k_sigma)]), end
    end
    if params.OMIT_LOW_CORRCOEF_BEATS
        [~, peaks_OMIT_LOW_CORRCOEF_BEATS] = refine_peaks_waveform_similarity(data_filtered_env, peak_indexes, pparams, 'BEAT-STD', params.PLOT_DIAGNOSTIC);
        %[~, peaks_OMIT_LOW_CORRCOEF_BEATS] = refine_peaks_waveform_similarity(data_filtered_mn_all_channels, peak_indexes, pparams, 'BEAT-STD', params.PLOT_DIAGNOSTIC);
        peaks = cat(1, peaks, peaks_OMIT_LOW_CORRCOEF_BEATS);
        refinement_methods = cat(2, refinement_methods, 'OMIT_LOW_CORRCOEF_BEATS');
    end

    % remove excess beats based on waveform similarity
    if ~isfield(params, 'OMIT_LOW_CORR_BEATS')
        params.OMIT_LOW_CORR_BEATS = true;
        pparams.percentile = 90.0;
        pparams.percentile_fraction = 0.3;
        if params.verbose, disp(['params.OMIT_LOW_CORR_BEATS = ', num2str(params.OMIT_LOW_CORR_BEATS)]), end
        if params.verbose, disp(['   pparams.percentile = ', num2str(pparams.percentile)]), end
        if params.verbose, disp(['   pparams.percentile_fraction = ', num2str(pparams.percentile_fraction)]), end
    end
    if params.OMIT_LOW_CORR_BEATS
        [~, peaks_OMIT_LOW_CORR_BEATS] = refine_peaks_waveform_similarity(data_filtered_mn_all_channels, peak_indexes, pparams, 'CORR', params.PLOT_DIAGNOSTIC);
        peaks = cat(1, peaks, peaks_OMIT_LOW_CORR_BEATS);
        refinement_methods = cat(2, refinement_methods, 'OMIT_LOW_CORR_BEATS');
    end

    % remove excess beats based on waveform similarity
    if ~isfield(params, 'OMIT_LOW_CORRCOEF_BEATS')
        params.OMIT_LOW_CORRCOEF_BEATS = true;
        pparams.percentile = 90.0;
        pparams.percentile_fraction = 0.3;
        if params.verbose, disp(['params.OMIT_LOW_CORRCOEF_BEATS = ', num2str(params.OMIT_LOW_CORRCOEF_BEATS)]), end
        if params.verbose, disp(['   pparams.percentile = ', num2str(pparams.percentile)]), end
        if params.verbose, disp(['   pparams.percentile_fraction = ', num2str(pparams.percentile_fraction)]), end
    end
    if params.OMIT_LOW_CORRCOEF_BEATS
        [~, peaks_OMIT_LOW_CORRCOEF_BEATS] = refine_peaks_waveform_similarity(data_filtered_mn_all_channels, peak_indexes, pparams, 'CORRCOEF', params.PLOT_DIAGNOSTIC);
        peaks = cat(1, peaks, peaks_OMIT_LOW_CORRCOEF_BEATS);
        refinement_methods = cat(2, refinement_methods, 'OMIT_LOW_CORRCOEF_BEATS');
    end

    if ~isfield(params, 'OMIT_LOW_AMP_PEAKS_PRCTL_ABS')
        params.OMIT_LOW_AMP_PEAKS_PRCTL_ABS = true;
        if params.verbose, disp('params.OMIT_LOW_AMP_PEAKS_PRCTL_ABS = true'), end
        if ~isfield(params, 'peak_amps_hist_prctile')
            pparams.peak_amps_hist_prctile = 25.0;
            if params.verbose, disp('   pparams.peak_amps_hist_prctile = 25.0'), end
        end
    end
    if params.OMIT_LOW_AMP_PEAKS_PRCTL_ABS
        [~, peaks_OMIT_LOW_AMP_PEAKS_PRCTL_ABS] = refine_peaks_low_amp_peaks_prctile(data_filtered_env, peak_indexes, pparams.peak_amps_hist_prctile, params.PLOT_DIAGNOSTIC);
        peaks = cat(1, peaks, peaks_OMIT_LOW_AMP_PEAKS_PRCTL_ABS);
        refinement_methods = cat(2, refinement_methods, 'OMIT_LOW_AMP_PEAKS_PRCTL_ABS');
    end

    % Omit beats that increase average beat SNR and increase HR variance. not applicable to the first and last beats
    if ~isfield(params, 'OMIT_BEAT_HRV_INCR_BEATS')
        params.OMIT_BEAT_HRV_INCR_BEATS = true;
        if params.verbose, disp(['params.OMIT_BEAT_HRV_INCR_BEATS = ', num2str(params.OMIT_BEAT_HRV_INCR_BEATS)]), end
    end
    if params.OMIT_BEAT_HRV_INCR_BEATS
        mmode = 'HEARTRATE'; % 'MORPHOLOGY' or 'HEARTRATE' or 'MORPHOLOGY-HEARTRATE'
        if params.verbose, disp(['   mmode = ', mmode]), end
        % [~, peaks_OMIT_BEAT_HRV_INCR_BEATS] = refine_peaks_low_snr_beats(data_filtered_mn_all_channels, peak_indexes, mmode, params.PLOT_DIAGNOSTIC);
        [~, peaks_OMIT_BEAT_HRV_INCR_BEATS] = refine_peaks_low_snr_beats(data_filtered_env, peak_indexes, mmode, params.PLOT_DIAGNOSTIC);
        peaks = cat(1, peaks, peaks_OMIT_BEAT_HRV_INCR_BEATS);
        refinement_methods = cat(2, refinement_methods, 'OMIT_BEAT_HRV_INCR_BEATS');
    end

    % Omit beats that increase average beat SNR and increase HR variance. not applicable to the first and last beats
    if ~isfield(params, 'OMIT_BEAT_SNR_REDUC_BEATS')
        params.OMIT_BEAT_SNR_REDUC_BEATS = true;
        if params.verbose, disp(['params.OMIT_BEAT_SNR_REDUC_BEATS = ', num2str(params.OMIT_BEAT_SNR_REDUC_BEATS)]), end
    end
    if params.OMIT_BEAT_SNR_REDUC_BEATS
        mmode = 'MORPHOLOGY'; % 'MORPHOLOGY' or 'HEARTRATE' or 'MORPHOLOGY-HEARTRATE'
        if params.verbose, disp(['   mmode = ', mmode]), end
        [~, peaks_OMIT_BEAT_SNR_REDUC_BEATS] = refine_peaks_low_snr_beats(data_filtered_mn_all_channels, peak_indexes, mmode, params.PLOT_DIAGNOSTIC);
        % [~, peaks_OMIT_BEAT_SNR_REDUC_BEATS] = refine_peaks_low_snr_beats(data_filtered_env, peak_indexes, mmode, params.PLOT_DIAGNOSTIC);
        peaks = cat(1, peaks, peaks_OMIT_BEAT_SNR_REDUC_BEATS);
        refinement_methods = cat(2, refinement_methods, 'OMIT_BEAT_SNR_REDUC_BEATS');
    end

    % remove excess beats based on ampliture thresholding (removes below a fraction of the defined percentile)
    if ~isfield(params, 'OMIT_HIGH_STD_AMP')
        params.OMIT_HIGH_STD_AMP = true;
        if params.verbose, disp(['params.OMIT_HIGH_STD_AMP = ', num2str(params.OMIT_HIGH_STD_AMP)]), end
        pparams.k_sigma = 4.0;
        if params.verbose, disp(['   pparams.k_sigma = ', num2str(pparams.k_sigma)]), end
    end
    if params.OMIT_HIGH_STD_AMP
        [~, peaks_OMIT_HIGH_STD_AMP] = refine_peaks_high_amp_std(data_filtered_env, peak_indexes, pparams.k_sigma, params.PLOT_DIAGNOSTIC);
        peaks = cat(1, peaks, peaks_OMIT_HIGH_STD_AMP);
        refinement_methods = cat(2, refinement_methods, 'OMIT_HIGH_STD_AMP');
    end
end

%% merge the peaks through consensus
peaks_consensus = sum(peaks, 1);
num_peak_refinement_algorithms = size(peaks, 1);

if ~isfield(params, 'NUM_VOTES_TO_KEEP_PEAK')
    params.NUM_VOTES_TO_KEEP_PEAK = ceil(num_peak_refinement_algorithms / 2); % majority vote
    if params.verbose, disp(['num_peak_refinement_algorithms = ', num2str(num_peak_refinement_algorithms), ', params.NUM_VOTES_TO_KEEP_PEAK = ', num2str(params.NUM_VOTES_TO_KEEP_PEAK)]), end
end

if params.NUM_VOTES_TO_KEEP_PEAK > num_peak_refinement_algorithms
    error('Number of required votes to keep a peak exceeds the number of voting algorithms');
end
peak_indexes_consensus = peaks_consensus >= params.NUM_VOTES_TO_KEEP_PEAK;


%% Likelihood-based refinement of the R-peaks (if required)
% search for local peaks with sign of most frequent among previously found peaks
if ~isfield(params, 'CORRECT_PEAKS_LIKELIHOOD_PEAKS')
    params.CORRECT_PEAKS_LIKELIHOOD_PEAKS = false;
    if params.verbose, disp(['params.CORRECT_PEAKS_LIKELIHOOD_PEAKS = ', num2str(params.CORRECT_PEAKS_LIKELIHOOD_PEAKS)]), end
    if ~isfield(params, 'likelihood_power_env_hist_peak_th')
        params.likelihood_power_env_hist_peak_th = 85.0;
        if params.verbose, disp(['params.likelihood_power_env_hist_peak_th = ', num2str(params.likelihood_power_env_hist_peak_th)]), end
    end
end
if params.CORRECT_PEAKS_LIKELIHOOD_PEAKS
    qrs_likelihood = peak_surrounding_likelihood(sig_len, peak_indexes_consensus, fs, params.max_likelihood_span, params.max_likelihood_span);
    % bumps_indexes = refine_peaks_low_amp_peaks_prctile(data_filtered_mn_all_channels.* qrs_likelihood, 1:sig_len, params.likelihood_power_env_hist_peak_th, params.PLOT_DIAGNOSTIC);
    bumps_indexes = refine_peaks_low_amp_peaks_prctile(data_filtered_env.* qrs_likelihood, 1:sig_len, params.likelihood_power_env_hist_peak_th, params.PLOT_DIAGNOSTIC);

% search for all local peaks within a given sliding window length
if ~isfield(params, 'min_peak_distance')
    params.min_peak_distance = 0.18;
    if params.verbose, disp(['params.min_peak_distance = ', num2str(params.min_peak_distance)]), end
end
rpeak_search_half_wlen = floor(fs * params.min_peak_distance);
env_pk_detect_mode = 'POS';
if params.verbose, disp(['env_pk_detect_mode = ', env_pk_detect_mode]), end
peak_indexes_consensus = refine_peaks_too_close_low_amp(data_filtered_env, bumps_indexes, rpeak_search_half_wlen, env_pk_detect_mode, params.PLOT_DIAGNOSTIC);
end

%% replace envelope peaks with original signal peaks, if required
if ~isfield(params, 'RETURN_SIGNAL_PEAKS')
    params.RETURN_SIGNAL_PEAKS = true;
    if params.verbose, disp(['params.RETURN_SIGNAL_PEAKS = ', num2str(params.RETURN_SIGNAL_PEAKS )]), end
end
if params.RETURN_SIGNAL_PEAKS
    if ~isfield(params, 'PEAK_SIGN')
        params.PEAK_SIGN = 'AUTO';
        if params.verbose, disp(['params.PEAK_SIGN = ', params.PEAK_SIGN]), end
    end
    if ~isfield(params, 'envelope_to_peak_search_wlen')
        params.envelope_to_peak_search_wlen = 0.05;
        if params.verbose, disp(['   params.envelope_to_peak_search_wlen = ', num2str(params.envelope_to_peak_search_wlen)]), end
    end
    envelope_to_peak_search_wlen = floor(fs * params.envelope_to_peak_search_wlen / 2);
    peak_indexes = find_closest_peaks(abs(data), peak_indexes, envelope_to_peak_search_wlen, params.PEAK_SIGN, params.PLOT_DIAGNOSTIC);
    % peak_indexes_consensus = find_closest_peaks(data_filtered, peak_indexes_consensus, envelope_to_peak_search_wlen, params.PEAK_SIGN, params.PLOT_DIAGNOSTIC);
end

qrs_likelihood = peak_surrounding_likelihood(sig_len, peak_indexes_consensus, fs, params.max_likelihood_span, params.max_likelihood_span);

% post-extraction peak refinement based on likelihoods
if ~isfield(params, 'POST_EXT_LIKELIHOOD_BASED_IMPROVEMENT')
    params.POST_EXT_LIKELIHOOD_BASED_IMPROVEMENT = false;
    if params.verbose, disp(['params.POST_EXT_LIKELIHOOD_BASED_IMPROVEMENT = ', num2str(params.POST_EXT_LIKELIHOOD_BASED_IMPROVEMENT)]), end
end
if params.POST_EXT_LIKELIHOOD_BASED_IMPROVEMENT
    likelihood_threshold = 0.4;
    if params.verbose, disp(['   likelihood_threshold = ', num2str(likelihood_threshold)]), end
    peak_indexes_consensus = refine_peaks_low_likelihood(data_filtered, peak_indexes_consensus, qrs_likelihood, likelihood_threshold, rpeak_search_half_wlen, params.PLOT_DIAGNOSTIC);
    qrs_likelihood = peak_surrounding_likelihood(sig_len, peak_indexes_consensus, fs, params.max_likelihood_span, params.max_likelihood_span);
end

if isfield(params, 'PLOT_RESULTS') && isequal(params.PLOT_RESULTS, 1)
    lgnds = {};
    tt = (0 : length(data) - 1) / fs;
    figure('units','normalized','outerposition',[0 0.25 1 0.5])
    plot(tt, data); lgnds = cat(2, lgnds, 'data');
    hold on
    plot(tt, data_filtered); lgnds = cat(2, lgnds, 'data_filtered');
    plot(tt, data_filtered_env); lgnds = cat(2, lgnds, 'data_filtered_env');
    plot(tt(bumps_indexes), data_filtered_env(bumps_indexes), 'g.', 'markersize', 14); lgnds = cat(2, lgnds, 'bumps_indexes');
    grid
    all_marks = {'o','+','*','.','x','s','d','>','v','<','^','p','h'};
    lgnds = [lgnds, refinement_methods];
    for ll = 1 : size(peaks, 1)
        pk_indx = find(peaks(ll, :));
        plot(tt(pk_indx), data_filtered_env(pk_indx), all_marks{mod(ll - 1, size(peaks, 1)) + 1}, 'markersize', ll + 10);
    end
    plot(tt(peak_indexes), data(peak_indexes), 'co', 'markersize', 16); lgnds = cat(2, lgnds, 'peak_indexes');
    plot(tt(peak_indexes_consensus), data(peak_indexes_consensus), 'ko', 'MarkerFaceColor','r', 'markersize', 20); lgnds = cat(2, lgnds, 'peak_indexes_consensus');
    legend(lgnds, 'interpreter', 'none', 'Location','eastoutside');
    xlabel('time[s]');
    ylabel('Amplitude');
    set(gca, 'fontsize', 16)
end

% return final peaks
peaks = zeros(1, sig_len);
peaks(peak_indexes) = 1;

% clear consenus-based peaks, if no refinement was requested
if ~params.REFINE_PEAKS
    peak_indexes_consensus = [];
end

end

%% INTERNAL FUNCTIONS

%% remove local peaks, which have nearby peaks with absolute higher amplitudes
function [peak_indexes, peaks] = find_closest_peaks(data, peak_indexes_candidates, peak_search_half_wlen, mode, plot_results)
sig_len = length(data);
polarity = sign(skew(data));
peak_indexes = [];
for jj = 1 : length(peak_indexes_candidates)
    segment = max(1, peak_indexes_candidates(jj) - peak_search_half_wlen) : min(sig_len, peak_indexes_candidates(jj) + peak_search_half_wlen);
    segment_first_index = segment(1);
    switch mode
        case 'AUTO'
            [~, I_max_min] = max(polarity * data(segment));
            pk_indx = I_max_min + segment_first_index - 1;
        case 'POS'
            [~, I_max] = max(data(segment));
            pk_indx = I_max + segment_first_index - 1;
        case 'NEG'
            [~, I_min] = min(data(segment));
            pk_indx = I_min + segment_first_index - 1;
        otherwise
            error('Undefined peak sign detection mode.');
    end
    peak_indexes = cat(2, peak_indexes, pk_indx);
end
peaks = zeros(1, sig_len);
peaks(peak_indexes) = 1;
if plot_results
    n = (1 : length(data));
    figure('units','normalized','outerposition',[0.1 0.25 0.8 0.5])
    plot(n, data)
    hold on
    plot(n(peak_indexes_candidates), data(peak_indexes_candidates), 'gx', 'markersize', 18)
    plot(n(peak_indexes), data(peak_indexes), 'ro', 'markersize', 18)
    legend('data', 'input peaks', 'refined peaks', 'Location','eastoutside')
    title('find_closest_peaks', 'interpreter', 'none');
    set(gca, 'fontsize', 18)
    grid
end
end

%% remove lower-amplitude peaks within a minimal window size
function [peak_indexes, peaks] = refine_peaks_too_close_low_amp(data, peak_indexes_candidates, peak_search_half_wlen, mode, plot_results)
sig_len = length(data);
peak_indexes = [];
for jj = 1 : length(peak_indexes_candidates)
    pk_index_candidate = peak_indexes_candidates(jj);
    segment = max(1, pk_index_candidate - peak_search_half_wlen) : min(sig_len, pk_index_candidate + peak_search_half_wlen);
    switch mode
        case 'POS'
            if max(data(segment)) == data(pk_index_candidate)
                peak_indexes = cat(2, peak_indexes, pk_index_candidate);
            end
        case 'NEG'
            if min(data(segment)) == data(pk_index_candidate)
                peak_indexes = cat(2, peak_indexes, pk_index_candidate);
            end
    end
end
peaks = zeros(1, sig_len);
peaks(peak_indexes) = 1;
if plot_results
    n = (1 : length(data));
    figure('units','normalized','outerposition',[0.1 0.25 0.8 0.5])
    plot(n, data)
    hold on
    plot(n(peak_indexes_candidates), data(peak_indexes_candidates), 'gx', 'markersize', 18)
    plot(n(peak_indexes), data(peak_indexes), 'ro', 'markersize', 18)
    legend('data', 'input peaks', 'refined peaks', 'Location','eastoutside')
    title('refine_peaks_too_close_low_amp', 'interpreter', 'none');
    set(gca, 'fontsize', 18)
    grid
end
end

%% remove excess beats based on ampliture thresholding (removes below the given percentile)
function [peak_indexes_refined, peaks] = refine_peaks_low_amp_peaks_prctile(data_env, peak_indexes, percentile, plot_results)
bumps_amp_threshold = prctile(data_env(peak_indexes), percentile);
peak_indexes_refined = peak_indexes(data_env(peak_indexes) >= bumps_amp_threshold);
peaks = zeros(1, length(data_env));
peaks(peak_indexes_refined) = 1;
if plot_results
    n = (1 : length(data_env));
    figure('units','normalized','outerposition',[0.1 0.25 0.8 0.5])
    plot(n, data_env)
    hold on
    plot(n(peak_indexes), data_env(peak_indexes), 'gx', 'markersize', 18)
    plot(n(peak_indexes_refined), data_env(peak_indexes_refined), 'ro', 'markersize', 18)
    legend('data', 'input peaks', 'refined peaks', 'Location','eastoutside')
    title('refine_peaks_low_amp_peaks_prctile', 'interpreter', 'none');
    set(gca, 'fontsize', 18)
    grid
end

end

%% remove excess beats based on ampliture thresholding (removes below a fraction of the defined percentile)
function [peak_indexes_refined, peaks] = refine_peaks_low_amp_peaks_prctile_fraction(data, peak_indexes, pparams, plot_results)
peak_indexes_refined = peak_indexes;
peak_amps = data(peak_indexes);
I_omit = peak_amps < pparams.percentile_fraction * prctile(peak_amps, pparams.percentile);
peak_indexes_refined(I_omit) = [];
peaks = zeros(1, length(data));
peaks(peak_indexes_refined) = 1;
if plot_results
    n = (1 : length(data));
    figure('units','normalized','outerposition',[0.1 0.25 0.8 0.5])
    plot(n, data)
    hold on
    plot(n(peak_indexes), data(peak_indexes), 'gx', 'markersize', 18)
    plot(n(peak_indexes_refined), data(peak_indexes_refined), 'ro', 'markersize', 18)
    legend('data', 'input peaks', 'refined peaks', 'Location','eastoutside')
    title('refine_peaks_low_amp_peaks_prctile_fraction', 'interpreter', 'none');
    set(gca, 'fontsize', 18)
    grid
end
end

%% remove excess beats based on ampliture thresholding (removes below a fraction of the defined percentile)
function [peak_indexes_refined, peaks] = refine_peaks_high_amp_std(data, peak_indexes, k_sigma, plot_results)
peak_indexes_refined = peak_indexes;
peak_amps = data(peak_indexes);
I_omit = abs(peak_amps - mean(peak_amps)) > k_sigma * std(peak_amps);
peak_indexes_refined(I_omit) = [];
peaks = zeros(1, length(data));
peaks(peak_indexes_refined) = 1;
if plot_results
    n = (1 : length(data));
    figure('units','normalized','outerposition',[0.1 0.25 0.8 0.5])
    plot(n, data)
    hold on
    plot(n(peak_indexes), data(peak_indexes), 'gx', 'markersize', 18)
    plot(n(peak_indexes_refined), data(peak_indexes_refined), 'ro', 'markersize', 18)
    legend('data', 'input peaks', 'refined peaks', 'Location','eastoutside')
    title('refine_peaks_high_amp_std', 'interpreter', 'none');
    set(gca, 'fontsize', 18)
    grid
end
end

%% remove excess beats based on waveform similarity
function [peak_indexes_refined, peaks] = refine_peaks_waveform_similarity(data, peak_indexes, pparams, method, plot_results)
event_width = 2 * round(median(diff(peak_indexes))/2) + 1;
stacked_beats = event_stacker(data, peak_indexes, event_width);
% Below, the average correlation of each beat with the others (excludes the beat itself)
switch method
    case 'CORR'
        rho_beats = stacked_beats * stacked_beats';
        avg_beat_corr_with_others = median(rho_beats + diag(nan(1, size(rho_beats, 1))), "omitmissing");
        I_omit = avg_beat_corr_with_others < pparams.percentile_fraction * prctile(avg_beat_corr_with_others, pparams.percentile);
    case 'CORRCOEF'
        rho_beats = corrcoef(stacked_beats');
        avg_beat_corr_with_others = median(rho_beats + diag(nan(1, size(rho_beats, 1))), "omitmissing");

        % I_omit = avg_beat_corr_with_others < min(pparams.max_corr_coef, pparams.max_corr_coef_fraction * max(avg_beat_corr_with_others));
        I_omit = avg_beat_corr_with_others < pparams.percentile_fraction * prctile(avg_beat_corr_with_others, pparams.percentile);
    case 'ABS-CORRCOEF'
        rho_beats = corrcoef(stacked_beats');
        avg_beat_corr_with_others = mean(rho_beats + diag(nan(1, size(rho_beats, 1))), "omitmissing");
        I_omit = avg_beat_corr_with_others < pparams.beat_corrcoef_th; % omit beats with low correlations
    case 'NEG-CORR'
        rho_beats = corrcoef(stacked_beats');
        avg_beat_corr_with_others = mean(rho_beats + diag(nan(1, size(rho_beats, 1))), "omitmissing");
        I_omit = avg_beat_corr_with_others < 0;
    case 'BEAT-STD'
        rho_beats = corrcoef(stacked_beats');
        avg_beat_corr_with_others = mean(rho_beats + diag(nan(1, size(rho_beats, 1))), "omitmissing");
        I_omit = abs(avg_beat_corr_with_others - mean(avg_beat_corr_with_others)) > pparams.k_sigma * std(avg_beat_corr_with_others); %%%% & avg_beat_corr_with_others < pparams.beat_corrcoef_th;
    otherwise
        error('undefined mode')
end
peak_indexes_refined = peak_indexes;
peak_indexes_refined(I_omit) = [];
peaks = zeros(1, length(data));
peaks(peak_indexes_refined) = 1;
if plot_results
    n = (1 : length(data));
    figure('units','normalized','outerposition',[0.1 0.25 0.8 0.5])
    plot(n, data)
    hold on
    plot(n(peak_indexes), data(peak_indexes), 'gx', 'markersize', 18)
    plot(n(peak_indexes_refined), data(peak_indexes_refined), 'ro', 'markersize', 18)
    legend('data', 'input peaks', 'refined peaks', 'Location','eastoutside')
    title(['refine_peaks_waveform_similarity: ', method], 'interpreter', 'none');
    set(gca, 'fontsize', 18)
    grid
end
end

%% remove low-power beats
function [peak_indexes_refined, peaks] = refine_peaks_low_power_beats(data, peak_indexes, max_amp_prctile, beat_std_med_frac_th, plot_results)
% peak_amps = data(peak_indexes);
event_width = 2 * round(median(diff(peak_indexes))/2) + 1;
[stacked_beats, num_non_zeros] = event_stacker(data, peak_indexes, event_width);
std_beats = (event_width -1) * std(stacked_beats, 0, 2)' ./ (num_non_zeros - 1); % compensates the boundary events STD that were zero-padded in event_stacker
% I_omit = std_beats < beat_std_med_frac_th*median(std_beats) & abs(peak_amps - mean(peak_amps)) > max_amp_k_sigma * std(peak_amps);
I_omit = std_beats < beat_std_med_frac_th * prctile(std_beats, max_amp_prctile); % & abs(peak_amps - mean(peak_amps)) > max_amp_k_sigma * std(peak_amps);
peak_indexes_refined = peak_indexes;
peak_indexes_refined(I_omit) = [];
peaks = zeros(1, length(data));
peaks(peak_indexes_refined) = 1;
if plot_results
    n = (1 : length(data));
    figure('units','normalized','outerposition',[0.1 0.25 0.8 0.5])
    plot(n, data)
    hold on
    plot(n(peak_indexes), data(peak_indexes), 'gx', 'markersize', 18)
    plot(n(peak_indexes_refined), data(peak_indexes_refined), 'ro', 'markersize', 18)
    legend('data', 'input peaks', 'refined peaks', 'Location','eastoutside')
    title('refine_peaks_low_power_beats', 'interpreter', 'none');
    set(gca, 'fontsize', 18)
    grid
end
end

%% peak refinement based on likelihoods
function [peak_indexes_refined, peaks] = refine_peaks_low_likelihood(data, peak_indexes, qrs_likelihood, likelihood_threshold, peak_search_half_wlen, plot_results)
sig_len = length(data);
signal_abs = sqrt(mean(data.^2, 1));
signal_likelihood = qrs_likelihood .* signal_abs/max(signal_abs);
bumps_indexes = find(signal_likelihood >= likelihood_threshold);
peak_indexes_refined = [];
for jj = 1 : length(bumps_indexes)
    index = bumps_indexes(jj);
    segment = max(1, index - peak_search_half_wlen) : min(sig_len, index + peak_search_half_wlen);
    if max(signal_likelihood(segment)) == signal_likelihood(index)
        peak_indexes_refined = cat(2, peak_indexes_refined, index);
    end
end

peak_indexes_refined = intersect(peak_indexes_refined, peak_indexes);

peaks = zeros(1, sig_len);
peaks(peak_indexes_refined) = 1;
if plot_results
    n = (1 : length(data));
    figure('units','normalized','outerposition',[0.1 0.25 0.8 0.5])
    plot(n, data)
    hold on
    plot(n(peak_indexes), data(peak_indexes), 'gx', 'markersize', 18)
    plot(n(peak_indexes_refined), data(peak_indexes_refined), 'ro', 'markersize', 18)
    legend('data', 'input peaks', 'refined peaks', 'Location','eastoutside')
    title('refine_peaks_low_likelihood', 'interpreter', 'none');
    set(gca, 'fontsize', 18)
    grid
end
end

%% omit beats that increase average beat SNR and increase HR variance. not applicable to the first and last beats
function [peak_indexes_refined, peaks] = refine_peaks_low_snr_beats(data, peak_indexes, mmode, plot_results)
max_itr = length(peak_indexes);
event_width = 2 * round(median(diff(peak_indexes))/2) + 1;
stacked_beats = event_stacker(data, peak_indexes, event_width);
%     ECG_robust_mean = RWAverage(stacked_beats);
ECG_robust_mean = mean(stacked_beats, 1);
ECG_robust_mean_replicated = ones(length(peak_indexes), 1) * ECG_robust_mean;
noise = stacked_beats - ECG_robust_mean_replicated;
%snr = 10*log10(trace(ECG_robust_mean_replicated*ECG_robust_mean_replicated')/trace(noise*noise'));
snr_initial = 20*log10(norm(ECG_robust_mean_replicated, 'fro')/norm(noise,'fro'));

% peak_indexes_refined = peak_indexes;
included_indexes = 1 : length(peak_indexes);
for itr = 1 : max_itr
    num_includes_beats = length(peak_indexes(included_indexes));
    rr_std = std(diff(peak_indexes(included_indexes)));
    snr_excluding_this_beat = zeros(1, num_includes_beats);
    rr_std_excluding_this_beat = zeros(1, num_includes_beats);
    for p = 1 : num_includes_beats
        all_included_indexes_but_this_beat = included_indexes;
        % this_beat_index = find(peak_indexes == peak_indexes(included_indexes(p)), 1,'first');
        all_included_indexes_but_this_beat(p) = [];

        % ECG_robust_mean = robust_weighted_average(stacked_beats(all_included_indexes_but_this_beat, :));
        ECG_robust_mean = mean(stacked_beats(all_included_indexes_but_this_beat, :), 1);
        ECG_robust_mean_replicated = ones(length(all_included_indexes_but_this_beat), 1) * ECG_robust_mean;

        % estimate channel quality SNR
        noise = stacked_beats(all_included_indexes_but_this_beat, :) - ECG_robust_mean_replicated;

        %snr_excluding_this_beat(p) = 10*log10(trace(ECG_robust_mean_replicated*ECG_robust_mean_replicated')/trace(noise*noise'));
        snr_excluding_this_beat(p) = 20*log10(norm(ECG_robust_mean_replicated, 'fro')/norm(noise,'fro'));

        % rr_std_excluding_this_beat(p) = max(abs(diff(peak_indexes(all_but_this_beat_index))));
        % rr_std_excluding_this_beat(p) = std(diff([1, peak_indexes(all_but_this_beat_index), sig_len]));
        rr_std_excluding_this_beat(p) = std(diff(peak_indexes(all_included_indexes_but_this_beat)));
    end
    [snr_excluding_worse_beat, I_worse_beat] = max(snr_excluding_this_beat);
    % beat_index = find(peak_indexes == peak_indexes(included_indexes(I_worse_beat)), 1,'first');
    switch mmode
        case 'MORPHOLOGY' % 'MORPHOLOGY' or 'HEARTRATE'
            if snr_excluding_worse_beat > snr_initial
                included_indexes(I_worse_beat) = [];
            else
                break;
            end
        case 'HEARTRATE'
            if rr_std_excluding_this_beat(I_worse_beat) < rr_std
                included_indexes(I_worse_beat) = [];
            else
                break;
            end
        case 'MORPHOLOGY-HEARTRATE'
            if snr_excluding_worse_beat > snr_initial && rr_std_excluding_this_beat(I_worse_beat) < rr_std
                % && data_filtered_env(peak_indexes(I_snr_excluding_this_beat)) < (peak_amps_max + peak_amps_med)/2% && hr_std <= rr_std_excluding_this_beat(I_snr_excluding_this_beat)
                included_indexes(I_worse_beat) = [];
            else
                break;
            end
    end
end

% make sure the first and last indexes are not removed
% if included_indexes(1) > 1
%     included_indexes = cat(2, 1, included_indexes);
% end
% if included_indexes(end) < length(peak_indexes)
%     included_indexes = cat(2, included_indexes, length(peak_indexes));
% end

peak_indexes_refined = peak_indexes(included_indexes);
peaks = zeros(1, length(data));
peaks(peak_indexes_refined) = 1;

if plot_results
    n = (1 : length(data));
    figure('units','normalized','outerposition',[0.1 0.25 0.8 0.5])
    plot(n, data)
    hold on
    plot(n(peak_indexes), data(peak_indexes), 'gx', 'markersize', 18)
    plot(n(peak_indexes_refined), data(peak_indexes_refined), 'ro', 'markersize', 18)
    legend('data', 'input peaks', 'refined peaks', 'Location','eastoutside')
    title('refine_peaks_low_snr_beats', 'interpreter', 'none');
    set(gca, 'fontsize', 18)
    grid
end

end

%% matched filter using average beat shape
function [data_enhanced, data_enhanced_env] = signal_specific_matched_filter(data, peak_indexes)
event_width = 2 * round(median(diff(peak_indexes))/2) + 1;
sig_len = size(data, 2);
data_enhanced = zeros(size(data));
for ch = 1 : size(data, 1)
    stacked_beats = event_stacker(data(ch, :), peak_indexes, event_width);
    robust_mean = robust_weighted_average(stacked_beats);
    matched_filter_out = conv(robust_mean(end : -1 : 1), data(ch, :));
    lag = round(length(robust_mean)/2);
    data_enhanced(ch, :) = matched_filter_out(lag : sig_len + lag - 1);
    data_enhanced(ch, :) = std(data(ch, :)) * data_enhanced(ch, :) / std(data_enhanced(ch, :));
end
data_enhanced_env = sqrt(sum(data_enhanced.^2, 1));
end

%% likelihood of peaks as we move towards or away a peak
function qrs_likelihood = peak_surrounding_likelihood(sig_len, peak_indexes, fs, max_span, likelihood_sigma)
peaks = zeros(1, sig_len);
peaks(peak_indexes) = 1;
sig_len = length(peaks);
tt = -max_span/2 : 1 / fs : max_span/2;
template = exp(-tt.^2/(2 * likelihood_sigma ^ 2));
lag = round(length(template)/2);
qrs_likelihood = conv(template, peaks);
qrs_likelihood = qrs_likelihood(lag : sig_len + lag - 1);
% qrs_likelihood = filtfilt(template, sum(template), peaks);
end
