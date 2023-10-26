function [peaks, peak_indexes, qrs_likelihood] = peak_det_probabilistic(signal, fs, varargin)
% function [peaks, peak_indexes, qrs_likelihood] = peak_det_probabilistic(signal, fs, varargin)
% A probabilistic R-peak detector based on local peaks sorting
%
% Note: Under development
%
% Usage:
%   [peaks, qrs_likelihood, peak_indexes] = peak_det_probabilistic(signal, fs, params)
%
% Inputs:
%   signal: single or multichannel ECG signal with row-wise channels
%   fs: sampling frequency
%   params: a structure containing the algorithm parameters (all parameters have default values if not given as input)
%       params.saturate: saturate (1) the channels before processing or not (0). Default = 0
%       params.filter_type: preprocessing filter type ('MATCHED_FILTER', or 'BANDPASS_FILTER'); Default = 'BANDPASS_FILTER'
%       params.low_cutoff: lower cutoff frequency of R-peak detector bandpass filter (in Hz), when params.filter_type = 'BANDPASS_FILTER'. Default = 1
%       params.up_cutoff: upper cutoff frequency of R-peak detector bandpass filter (in Hz), when params.filter_type = 'BANDPASS_FILTER'. Default = 40
%       params.matched_filter_span: the preprocessing gaussian shape matched filter span, when params.filter_type = 'MATCHED_FILTER'. Default = 0.2
%       params.matched_filter_sigma: the preprocessing gaussian shape matched filter STD, when params.filter_type = 'MATCHED_FILTER'. Default = 0.01
%       params.wlen_power_env: power envelope moving average window length (in seconds). Default = 0.03
%       params.n_bins: number of local peaks histogram bins. Default = max(250, min(500, 10% of signal length))
%       params.hist_search_th: the top percentile of the signal's power envelope, considered for R-peak detection. Default = 0.9 (top 10%)
%       params.rpeak_search_wlen: the R-peak detector search window length (in seconds). Default = 0.2
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

sig_len = size(signal, 2); % signal length

%% saturate the channels at k_sigma times the channel-wise STD
if isfield(params, 'saturate') && isequal(params.saturate, 1)
    if ~isfield(params, 'k_sigma')
        params.k_sigma = 12.0;
    end
    data_sat = tanh_saturation(signal, params.k_sigma, 'ksigma');
else
    data_sat = signal;
end

%% pass the channels through a narrow bandpass filter or a matched filter to remove the baseline and the T-waves
if isfield(params, 'filter_type')
    switch params.filter_type
        case 'MATCHED_FILTER' % Matched filter
            if ~isfield(params, 'matched_filter_span')
                params.matched_filter_span = 0.2;
            end
            if ~isfield(params, 'matched_filter_sigma')
                params.matched_filter_sigma = 0.01;
            end
            t_matched = -params.matched_filter_span/2 : 1 / fs : params.matched_filter_span/2;
            template_matched = exp(-t_matched.^2/(2 * params.matched_filter_sigma ^ 2));
            lag_matched = round(length(template_matched)/2);
            data_filtered = zeros(size(data_sat));
            for kk = 1 : size(data_sat, 1)
                band_passsed = conv(template_matched, data_sat(kk, :))/sum(template_matched.^2);
                data_filtered(kk, :) = band_passsed(lag_matched : sig_len + lag_matched - 1);
            end
        case 'MULT_MATCHED_FILTER' % multi-template matched filter
            span = 0.3;
            sigma1 = 0.005;
            t_matched = -span/2 : 1 / fs : span/2;
            template_matched1 = exp(-t_matched.^2/(2 * sigma1 ^ 2));
            template_matched2 = diff(template_matched1);
            lag1 = round(length(template_matched1)/2);
            lag2 = round(length(template_matched2)/2);
            data_filtered = zeros(size(data_sat));
            for kk = 1 : size(data_sat, 1)
                band_passsed1 = conv(template_matched1, data_sat(kk, :))/sum(template_matched1.^2);
                band_passsed1 = band_passsed1(lag1 : sig_len + lag1 - 1);

                band_passsed2 = conv(template_matched2, data_sat(kk, :))/sum(template_matched2.^2);
                band_passsed2 = band_passsed2(lag2 : sig_len + lag2 - 1);

                data_filtered(kk, :) = sqrt(band_passsed1.^2 + band_passsed2.^2);
            end
        case 'MDMN' % Median + MA filter
            if ~isfield(params, 'wlen_md')
                params.wlen_md = 0.075;
            end
            if ~isfield(params, 'wlen_mn')
                params.wlen_mn = 0.05;
            end
            data_filtered = zeros(size(data_sat));
            wlen1 = round(params.wlen_md * fs);
            wlen2 = round(params.wlen_mn * fs);
            for kk = 1 : size(data_sat, 1)
                bl1 = baseline_sliding_window(data_sat(kk, :), wlen1, 'md');
                bl2 = baseline_sliding_window(bl1, wlen2, 'mn');
                data_filtered(kk, :) = data_sat(kk, :) - bl2;
                % figure
                % plot(data_sat(kk, :));
                % hold on
                % plot(bl1);
                % plot(bl2);
                % plot(data_filtered(kk, :));
                % grid
                % legend('datasat', 'bl1', 'bl2', 'filtered');
            end
        case  'WAVELET' % Wavelet denoiser
            if isfield(params, 'WDtype')
                wtype = params.WDtype;
            else
                wtype = 'sym4'; % mother wavelet used for wavelet denoising
            end

            if isfield(params, 'WDLevel1')
                level1 = params.WDLevel1;
            else
                level1 = floor(log2(fs/22.5)); % Upper frequency range for the wavelet denoiser
            end
            if isfield(params, 'WDLevel2')
                level2 = params.WDLevel2;
            else
                level2 = ceil(log2(fs/6.5)); % Lower frequency range for the wavelet denoiser
            end

            data_filtered = zeros(size(data_sat));
            for kk = 1 : size(data_sat, 1)
                wt = modwt(data_sat(kk, :), level2);
                wtrec = zeros(size(wt));
                wtrec(level1:level2,:) = wt(level1:level2,:);
                data_filtered(kk, :) = imodwt(wtrec, wtype);
            end
        case 'BANDPASS_FILTER' % Bandpass filter
            if ~isfield(params, 'low_cutoff')
                params.low_cutoff = 1.0;
            end
            if ~isfield(params, 'up_cutoff')
                params.up_cutoff = 40.0;
            end
            data_lp = lp_filter_zero_phase(data_sat, params.up_cutoff/fs);
            data_filtered = data_lp - lp_filter_zero_phase(data_lp, params.low_cutoff/fs);

        otherwise
            error('Unknown preprocessing filter');
    end
else % Bandpass filter
    if ~isfield(params, 'low_cutoff')
        params.low_cutoff = 1.0;
    end
    if ~isfield(params, 'up_cutoff')
        params.up_cutoff = 40.0;
    end
    warning(['Preprocessing filter type not specified. Using a generic bandpass filter: (', num2str(params.low_cutoff), ', ', num2str(params.up_cutoff), ') Hz.']);

    data_lp = lp_filter_zero_phase(data_sat, params.up_cutoff/fs);
    data_filtered = data_lp - lp_filter_zero_phase(data_lp, params.low_cutoff/fs);
end

data_filtered_mn_all_channels = mean(data_filtered, 1);

%% calculate the power envelope of one or all channels
if ~isfield(params, 'wlen_power_env')
    params.wlen_power_env = 0.05;
end
wlen_power_env = ceil(params.wlen_power_env * fs);
data_filtered_env = filtfilt(ones(wlen_power_env, 1), wlen_power_env, sqrt(mean(data_filtered.^2, 1)));

%% pick the top percentage of the signal's power envelope for R-peak search
if ~isfield(params, 'hist_search_th')
    params.hist_search_th = 0.9;
end
bumps_indexes = refine_peaks_low_amp_peaks_prctile(data_filtered_env, 1:sig_len, 100*params.hist_search_th);

%% search for all local peaks within a given sliding window length
if ~isfield(params, 'rpeak_search_wlen')
    params.rpeak_search_wlen = 0.2;
end
rpeak_search_half_wlen = floor(fs * params.rpeak_search_wlen);
peak_indexes = refine_peaks_too_close_low_amp(data_filtered_env, bumps_indexes, rpeak_search_half_wlen);

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
end
if ~isfield(params, 'max_likelihood_span')
    params.max_likelihood_span = 0.3;
end

% matched filter using average beat shape
if ~isfield(params, 'ENHANCE_MATCHED_FILTER')
    params.ENHANCE_MATCHED_FILTER = false;
end
if params.ENHANCE_MATCHED_FILTER
    [data_filtered_enhanced, data_filtered_enhanced_env] = signal_specific_matched_filter(data_filtered, peak_indexes);
    figure
    plot(data_filtered);
    hold on
    plot(data_filtered_enhanced);
    plot(data_filtered_enhanced_env);
    grid
end

%% Refine the extracted R-peaks
refinement_methods = {};
peaks = [];
if ~isfield(params, 'REMOVE_LOW_AMP_PEAKS')
    params.REMOVE_LOW_AMP_PEAKS = true;
end
% remove excess beats with extreme amplitude deviations from other peaks
if params.REMOVE_LOW_AMP_PEAKS
    percentile = 75.0;
    percentile_fraction = 0.5;
    [~, peaks1] = refine_peaks_low_amp_peaks_prctile_fraction(data_filtered_mn_all_channels, peak_indexes, percentile, percentile_fraction);
    peaks = cat(1, peaks, peaks1);
    refinement_methods = cat(2, refinement_methods, 'REMOVE_LOW_AMP_PEAKS');
end

% negative correlations
if ~isfield(params, 'REMOVE_NEG_CORR')
    params.REMOVE_NEG_CORR = true;
end
if params.REMOVE_NEG_CORR
    [~, peaks2] = refine_peaks_waveform_similarity(data_filtered_mn_all_channels, peak_indexes, [], 'NEG-CORR');
    peaks = cat(1, peaks, peaks2);
    refinement_methods = cat(2, refinement_methods, 'REMOVE_NEG_CORR');
end

% omit beats with low variance
if ~isfield(params, 'REMOVE_LOW_ENERGY_BEATS')
    params.REMOVE_LOW_ENERGY_BEATS = true;
end
if params.REMOVE_LOW_ENERGY_BEATS
    max_amp_k_sigma = 2.0;
    beat_std_media_frac_th = 0.2;
    % peak_indexes_refined3 = refine_peaks_low_variance(data_filtered_env, peak_indexes, max_amp_k_sigma, beat_std_media_frac_th);
    [~, peaks3] = refine_peaks_low_variance(data_filtered_mn_all_channels, peak_indexes, max_amp_k_sigma, beat_std_media_frac_th);
    peaks = cat(1, peaks, peaks3);
    refinement_methods = cat(2, refinement_methods, 'REMOVE_LOW_ENERGY_BEATS');
end

% omit beats with low correlations
if ~isfield(params, 'REMOVE_LOW_CORR_BEATS')
    params.REMOVE_LOW_CORR_BEATS = true;
end
if params.REMOVE_LOW_CORR_BEATS
    pparams.k_sigma = 3.0;
    pparams.beat_corr_th = 0.5;
    [~, peaks4] = refine_peaks_waveform_similarity(data_filtered_env, peak_indexes, pparams, 'BEAT-STD');
    %[~, peaks4] = refine_peaks_waveform_similarity(data_filtered_mn_all_channels, peak_indexes, pparams, 'BEAT-STD');
    peaks = cat(1, peaks, peaks4);
    refinement_methods = cat(2, refinement_methods, 'REMOVE_LOW_CORR_BEATS');
end

if ~isfield(params, 'REMOVE_LOW_AMP_PRCTILE')
    params.REMOVE_LOW_AMP_PRCTILE = true;
end
if params.REMOVE_LOW_CORR_BEATS
    [~, peaks5] = refine_peaks_low_amp_peaks_prctile(data_filtered_env, peak_indexes, 100*params.hist_search_th);
    peaks = cat(1, peaks, peaks5);
    refinement_methods = cat(2, refinement_methods, 'REMOVE_LOW_AMP_PRCTILE');
end

% remove excess beats based on waveform similarity
if ~isfield(params, 'REMOVE_LOW_WAVEFORM_SIMILARITY')
    params.REMOVE_LOW_WAVEFORM_SIMILARITY = true;
end
if params.REMOVE_LOW_WAVEFORM_SIMILARITY
    pparams.max_corr_coef = 0.9;
    pparams.max_corr_coef_fraction = 0.7;
    [~, peaks6] = refine_peaks_waveform_similarity(data_filtered_env, peak_indexes, pparams, 'LOW-CORR');
    peaks = cat(1, peaks, peaks6);
    refinement_methods = cat(2, refinement_methods, 'REMOVE_LOW_WAVEFORM_SIMILARITY');
end

% Omit beats that increase average beat SNR and increase HR variance. not applicable to the first and last beats
if ~isfield(params, 'REMOVE_HRV_REDUCING_BEATS')
    params.REMOVE_HRV_REDUCING_BEATS = true;
end
if params.REMOVE_HRV_REDUCING_BEATS
    [~, peaks7] = refine_peaks_low_snr_beats(data_filtered_mn_all_channels, peak_indexes);
    peaks = cat(1, peaks, peaks7);
    refinement_methods = cat(2, refinement_methods, 'REMOVE_HRV_REDUCING_BEATS');
end

% remove excess beats based on ampliture thresholding (removes below a fraction of the defined percentile)
if ~isfield(params, 'REMOVE_HIGH_AMP_STD')
    params.REMOVE_HIGH_AMP_STD = true;
    pparams.k_sigma = 4.0;
end
if params.REMOVE_HIGH_AMP_STD
    [~, peaks8] = refine_peaks_high_amp_std(data_filtered_mn_all_channels, peak_indexes, pparams.k_sigma);
    peaks = cat(1, peaks, peaks8);
    refinement_methods = cat(2, refinement_methods, 'REMOVE_HIGH_AMP_STD');
end

% search for local peaks with sign of most frequent among previously found peaks
% Likelihood-based refinement of the R-peaks
if ~isfield(params, 'LIKELIHOOD_BASED_REFINEMENT')
    params.LIKELIHOOD_BASED_REFINEMENT = true;
end
if params.LIKELIHOOD_BASED_REFINEMENT
    qrs_likelihood = peak_surrounding_likelihood(sig_len, peak_indexes, fs, params.max_likelihood_span, params.max_likelihood_span);
    [peaks9, ~] = peak_det_local_search(mode(sign(data_filtered_mn_all_channels(peak_indexes))) * data_filtered_mn_all_channels .* qrs_likelihood, 1.0/median(diff(peak_indexes)), 1);
    % [~, peak_indexes_refined9] = peak_det_local_search(data_filtered_env .* qrs_likelihood, 1.0/median(diff(peak_indexes)));
    % % [~, peak_indexes_refined9] = peak_det_local_search(mean(data_filtered, 1).*qrs_likelihood, 1.0/median(diff(peak_indexes)));
    % % [~, peak_indexes_refined9] = peak_det_local_search(mean(signal, 1) .* qrs_likelihood, 1.0/median(diff(peak_indexes)));
    % % [~, peak_indexes_refined9] = peak_det_local_search(mean(data_filtered_enhanced, 1) .* qrs_likelihood, 1.0/median(diff(peak_indexes)));
    % % [~, peak_indexes_refined9] = peak_det_local_search(sign(skew(mean(data_filtered, 1))) * mean(data_filtered, 1).*qrs_likelihood, 1.0/median(diff(peak_indexes)), 1);
    % nn = (0 : length(data_filtered_env) - 1)/fs;
    % figure
    % plot(nn, data_filtered_env)
    % hold on
    % plot(nn, qrs_likelihood .* data_filtered_env)
    % plot(nn(peak_indexes), data_filtered_env(peak_indexes), 'ro')
    % grid
    peaks = cat(1, peaks, peaks9);
    refinement_methods = cat(2, refinement_methods, 'LIKELIHOOD_BASED_REFINEMENT');
end

peaks_consensus = sum(peaks, 1);
peak_indexes_consensus = peaks_consensus >= ceil(size(peaks, 1) / 2);

% figure
% hold on
% grid
% plot(signal)
% plot(data_sat)
% plot(data_filtered)
% plot(data_filtered_env)
% plot(peak_indexes, data_filtered_env(peak_indexes), 'ko')
% plot(peak_indexes_refined1, data_filtered_env(peak_indexes_refined1), 'gx')
% plot(peak_indexes_refined2, data_filtered_env(peak_indexes_refined2), 'ro', 'markersize', 10)
% plot(peak_indexes_refined3, data_filtered_env(peak_indexes_refined3), 'bo', 'markersize', 12)
% plot(peak_indexes_refined4, data_filtered_env(peak_indexes_refined4), 'co', 'markersize', 18, 'linewidth', 3)
% legend({'signal', 'data_sat', 'data_filtered', 'data_filtered_env', 'peaks', 'peak_indexes_refined1', 'peak_indexes_refined2', 'peak_indexes_refined3', 'peak_indexes_refined4'},'interpreter', 'none');

if ~isfield(params, 'RETURN_SIGNAL_PEAKS')
    params.RETURN_SIGNAL_PEAKS = true;
end
if params.RETURN_SIGNAL_PEAKS
    if ~isfield(params, 'PEAK_SIGN')
        params.PEAK_SIGN = 'AUTO';
    end
    if ~isfield(params, 'envelope_to_peak_search_wlen')
        params.envelope_to_peak_search_wlen = 0.1;
    end
    envelope_to_peak_search_wlen = floor(fs * params.envelope_to_peak_search_wlen / 2);
    peak_indexes = find_closest_peaks(data_filtered, peak_indexes, envelope_to_peak_search_wlen, params.PEAK_SIGN);
end

qrs_likelihood = peak_surrounding_likelihood(sig_len, peak_indexes, fs, params.max_likelihood_span, params.max_likelihood_span);

% peak refinement based on likelihoods
if ~isfield(params, 'LIKELIHOOD_BASED_IMPROVEMENT')
    params.LIKELIHOOD_BASED_IMPROVEMENT = true;
end
if params.LIKELIHOOD_BASED_IMPROVEMENT
    likelihood_threshold = 0.4;
    peak_indexes = refine_peaks_low_likelihood(data_filtered, qrs_likelihood, likelihood_threshold, rpeak_search_half_wlen);
end

if isfield(params, 'PLOT_RESULTS') && isequal(params.PLOT_RESULTS, 1)
    % time = (0 : sig_len - 1) / fs;
    % figure
    % scatter(time, signal', 28, qrs_likelihood(:)*[1 0 0], 'filled');
    % hold on
    % plot(time, data_filtered_env);
    % % plot(time, signal_likelihood);
    % grid

    tt = (0 : length(signal) - 1) / fs;
    figure
    plot(tt, signal)
    hold on
    plot(tt, data_filtered_env)
    plot(tt(bumps_indexes), data_filtered_env(bumps_indexes), 'g.', 'markersize', 14)
    plot(tt(peak_indexes), data_filtered_env(peak_indexes), 'co', 'markersize', 16)
    plot(tt(peak_indexes_consensus), signal(peak_indexes_consensus), 'ko', 'MarkerFaceColor','r', 'markersize', 20)
    grid
    all_marks = {'o','+','*','.','x','s','d','^','v','>','<','p','h'};
    for ll = 1 : size(peaks, 1)
        pk_indx = find(peaks(ll, :));
        plot(tt(pk_indx), data_filtered_env(pk_indx), all_marks{mod(ll - 1, size(peaks, 1)) + 1}, 'markersize', ll + 10);
    end
    % refinement_methods
end


peaks = zeros(1, sig_len);
peaks(peak_indexes) = 1;
end

% find actual beats (finds the closest local peaks of the actual signal
% from energy envelopes of the signsl
function [peak_indexes, peaks] = find_closest_peaks(data, peak_indexes, rpeak_search_half_wlen, mode)
sig_len = length(data);
polarity = sign(skew(data));
peak_indexes_refined = [];
for jj = 1 : length(peak_indexes)
    segment = max(1, peak_indexes(jj) - rpeak_search_half_wlen) : min(sig_len, peak_indexes(jj) + rpeak_search_half_wlen);
    switch mode
        case 'AUTO'
            [~, I_max_min] = max(polarity * data(segment));
            pk_indx = I_max_min + segment(1) - 1;
        case 'POS'
            [~, I_max] = max(data(segment));
            pk_indx = I_max + segment(1) - 1;
        case 'NEG'
            [~, I_min] = min(data(segment));
            pk_indx = I_min + segment(1) - 1;
        otherwise
            error('Undefined peak sign detection mode.');
    end
    peak_indexes_refined = cat(2, peak_indexes_refined, pk_indx);
end
peak_indexes = peak_indexes_refined;
peaks = zeros(1, sig_len);
peaks(peak_indexes) = 1;
end

% remove lower-amplitude peaks within a minimal window size
function [peak_indexes, peaks] = refine_peaks_too_close_low_amp(data, bumps_indexes, rpeak_search_half_wlen)
sig_len = length(data);
peak_indexes = [];
for jj = 1 : length(bumps_indexes)
    bump_index = bumps_indexes(jj);
    segment = max(1, bump_index - rpeak_search_half_wlen) : min(sig_len, bump_index + rpeak_search_half_wlen);
    if max(data(segment)) == data(bump_index)
        peak_indexes = cat(2, peak_indexes, bump_index);
    end
end
peaks = zeros(1, sig_len);
peaks(peak_indexes) = 1;
end

% remove excess beats based on ampliture thresholding (removes below the given percentile)
function [peak_indexes_refined, peaks] = refine_peaks_low_amp_peaks_prctile(data_env, peak_indexes, percentile)
bumps_amp_threshold = prctile(data_env(peak_indexes), percentile);
peak_indexes_refined = peak_indexes(data_env(peak_indexes) >= bumps_amp_threshold);
peaks = zeros(1, length(data_env));
peaks(peak_indexes_refined) = 1;
end

% remove excess beats based on ampliture thresholding (removes below a fraction of the defined percentile)
function [peak_indexes_refined, peaks] = refine_peaks_low_amp_peaks_prctile_fraction(data, peak_indexes, percentile, percentile_fraction)
peak_indexes_refined = peak_indexes;
peak_amps = data(peak_indexes);
I_omit = abs(peak_amps) < percentile_fraction * prctile(abs(peak_amps), percentile);
peak_indexes_refined(I_omit) = [];
peaks = zeros(1, length(data));
peaks(peak_indexes_refined) = 1;
end

% remove excess beats based on ampliture thresholding (removes below a fraction of the defined percentile)
function [peak_indexes_refined, peaks] = refine_peaks_high_amp_std(data, peak_indexes, k_sigma)
peak_indexes_refined = peak_indexes;
peak_amps = data(peak_indexes);
I_omit = abs(peak_amps - mean(peak_amps)) > k_sigma * std(peak_amps);
peak_indexes_refined(I_omit) = [];
peaks = zeros(1, length(data));
peaks(peak_indexes_refined) = 1;
end

% remove excess beats based on waveform similarity
function [peak_indexes_refined, peaks] = refine_peaks_waveform_similarity(data, peak_indexes, params, method)
event_width = 2 * round(median(diff(peak_indexes))/2) + 1;
stacked_beats = event_stacker(data, peak_indexes, event_width);
rho_beats = corrcoef(stacked_beats');
C_beats = (sum(rho_beats, 1) - 1.0) / (length(peak_indexes) - 1); % Average correlation of each beat with the others
switch method
    case 'LOW-CORR'
        I_omit = C_beats < min(params.max_corr_coef, params.max_corr_coef_fraction * max(C_beats));
    case 'NEG-CORR'
        I_omit = C_beats < 0;
    case 'BEAT-STD'
        I_omit = abs(C_beats - mean(C_beats)) > params.k_sigma * std(C_beats) & C_beats < params.beat_corr_th; % omit beats with low correlations
end
peak_indexes_refined = peak_indexes;
peak_indexes_refined(I_omit) = [];
peaks = zeros(1, length(data));
peaks(peak_indexes_refined) = 1;
end

% remove low-variance beats
function [peak_indexes_refined, peaks] = refine_peaks_low_variance(data, peak_indexes, max_amp_k_sigma, beat_std_media_frac_th)
peak_amps = data(peak_indexes);
event_width = 2 * round(median(diff(peak_indexes))/2) + 1;
[stacked_beats, num_non_zeros] = event_stacker(data, peak_indexes, event_width);
std_beats = (event_width -1) * std(stacked_beats, 0, 2)' ./ (num_non_zeros - 1); % compensates the boundary events STD that were zero-padded in event_stacker
I_omit = std_beats < beat_std_media_frac_th*median(std_beats) & abs(peak_amps - mean(peak_amps)) > max_amp_k_sigma * std(peak_amps);
peak_indexes_refined = peak_indexes;
peak_indexes_refined(I_omit) = [];
peaks = zeros(1, length(data));
peaks(peak_indexes_refined) = 1;
end

% peak refinement based on likelihoods
function [peak_indexes, peaks] = refine_peaks_low_likelihood(data, qrs_likelihood, likelihood_threshold, rpeak_search_half_wlen)
sig_len = length(data);
signal_abs = sqrt(mean(data.^2, 1));
signal_likelihood = qrs_likelihood .* signal_abs/max(signal_abs);
bumps_indexes = find(signal_likelihood >= likelihood_threshold);
peak_indexes = [];
for jj = 1 : length(bumps_indexes)
    index = bumps_indexes(jj);
    segment = max(1, index - rpeak_search_half_wlen) : min(sig_len, index + rpeak_search_half_wlen);
    if max(signal_likelihood(segment)) == signal_likelihood(index)
        peak_indexes = cat(2, peak_indexes, index);
    end
end
peaks = zeros(1, sig_len);
peaks(peak_indexes) = 1;
end

% matched filter using average beat shape
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

% likelihood of peaks as we move towards or away a peak
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

% omit beats that increase average beat SNR and increase HR variance. not applicable to the first and last beats
function [peak_indexes, peaks] = refine_peaks_low_snr_beats(data, peak_indexes)
max_itr = length(peak_indexes);
for itr = 1 : max_itr
    num_beats = length(peak_indexes);
    event_width = 2 * round(median(diff(peak_indexes))/2) + 1;
    stacked_beats = event_stacker(data, peak_indexes, event_width);

    %     ECG_robust_mean = RWAverage(stacked_beats);
    ECG_robust_mean = mean(stacked_beats, 1);
    ECG_robust_mean_replicated = ones(num_beats, 1) * ECG_robust_mean;
    noise = stacked_beats - ECG_robust_mean_replicated;

    %snr = 10*log10(trace(ECG_robust_mean_replicated*ECG_robust_mean_replicated')/trace(noise*noise'));
    snr = 20*log10(norm(ECG_robust_mean_replicated, 'fro')/norm(noise,'fro'));

    % hr_std = max(abs(diff(peak_indexes)));
    % hr_std = std(diff([1, peak_indexes, sig_len]));
    hr_std = std(diff(peak_indexes));

    snr_excluding_this_beat = zeros(1, num_beats);
    hr_std_excluding_this_beat = zeros(1, num_beats);
    for p = 1 : num_beats
        all_but_this_beat_index = 1 : num_beats;
        all_but_this_beat_index(p) = [];

        % ECG_robust_mean = robust_weighted_average(stacked_beats(all_but_this_beat_index, :));
        ECG_robust_mean = mean(stacked_beats(all_but_this_beat_index, :), 1);
        ECG_robust_mean_replicated = ones(num_beats - 1, 1) * ECG_robust_mean;

        % estimate channel quality SNR
        noise = stacked_beats(all_but_this_beat_index, :) - ECG_robust_mean_replicated;

        %snr_excluding_this_beat(p) = 10*log10(trace(ECG_robust_mean_replicated*ECG_robust_mean_replicated')/trace(noise*noise'));
        snr_excluding_this_beat(p) = 20*log10(norm(ECG_robust_mean_replicated, 'fro')/norm(noise,'fro'));

        % hr_std_excluding_this_beat(p) = max(abs(diff(peak_indexes(all_but_this_beat_index))));
        % hr_std_excluding_this_beat(p) = std(diff([1, peak_indexes(all_but_this_beat_index), sig_len]));
        hr_std_excluding_this_beat(p) = std(diff(peak_indexes(all_but_this_beat_index)));
    end
    [snr_max_imporvement, I_max_snr_imporvement] = max(snr_excluding_this_beat);
    if snr_max_imporvement > snr && hr_std_excluding_this_beat(I_max_snr_imporvement) < hr_std ...
            && I_max_snr_imporvement > 1 && I_max_snr_imporvement < num_beats % && data_filtered_env(peak_indexes(I_max_snr_imporvement)) < (peak_amps_max + peak_amps_med)/2% && hr_std <= hr_std_excluding_this_beat(I_max_snr_imporvement)
        peak_indexes(I_max_snr_imporvement) = [];
    else
        break;
    end
end
peaks = zeros(1, length(data));
peaks(peak_indexes) = 1;
end


% if 1
%     f0 = 1/params.rpeak_search_wlen;
%     peaks = PeakDetection(data_filtered_env, f0/fs);
%     peak_indexes = find(peaks);
%     peak_amps = data_filtered_env(peak_indexes);
%     peak_amps_max = max(peak_amps);
%     peak_amps_min = min(peak_amps);
%     peak_amps_med = median(peak_amps);
%
%     figure
%     hold on
%     grid
%     plot(signal)
%     plot(data_sat)
%     plot(data_filtered)
%     plot(data_filtered_env)
%     plot(peak_indexes, data_filtered_env(peak_indexes), 'ko', 'markersize', 14, 'linewidth', 3)
%     for itr = 1 : 150
%         num_beats = length(peak_indexes);
%         event_width = 2 * round(0.75    *           median(diff(peak_indexes))/2) + 1;
%         %     event_width = 2 * round(fs * params.rpeak_search_wlen/2) + 1;
%
%
%
%         %     stacked_beats = EventStacker(data_filtered, peak_indexes, event_width);
%         stacked_beats = EventStacker(data_filtered_env, peak_indexes, event_width);
%
%
%
%         %     ECG_robust_mean = RWAverage(stacked_beats);
%         ECG_robust_mean = mean(stacked_beats, 1);
%         ECG_robust_mean_replicated = ones(num_beats, 1) * ECG_robust_mean;
%         noise = stacked_beats - ECG_robust_mean_replicated;
%         %snr = 10*log10(trace(ECG_robust_mean_replicated*ECG_robust_mean_replicated')/trace(noise*noise'));
%         snr = 20*log10(norm(ECG_robust_mean_replicated, 'fro')/norm(noise,'fro'));
%
%         %     hr_std = max(abs(diff(peak_indexes)));
%         hr_std = std(diff([1, peak_indexes, sig_len]));
%         snr_excluding_this_beat = zeros(1, num_beats);
%         hr_std_excluding_this_beat = zeros(1, num_beats);
%         for p = 1 : num_beats
%             all_but_this_beat_index = 1 : num_beats;
%             all_but_this_beat_index(p) = [];
%             %         ECG_robust_mean = RWAverage(stacked_beats(all_but_this_beat_index, :));
%             ECG_robust_mean = mean(stacked_beats(all_but_this_beat_index, :), 1);
%             ECG_robust_mean_replicated = ones(num_beats - 1, 1) * ECG_robust_mean;
%             % estimate channel quality SNR
%             noise = stacked_beats(all_but_this_beat_index, :) - ECG_robust_mean_replicated;
%             %snr_excluding_this_beat(p) = 10*log10(trace(ECG_robust_mean_replicated*ECG_robust_mean_replicated')/trace(noise*noise'));
%             snr_excluding_this_beat(p) = 20*log10(norm(ECG_robust_mean_replicated, 'fro')/norm(noise,'fro'));
%
%             %         hr_std_excluding_this_beat(p) = max(abs(diff(peak_indexes(all_but_this_beat_index))));
%             hr_std_excluding_this_beat(p) = std(diff([1, peak_indexes(all_but_this_beat_index), sig_len]));
%             %         if snr < snr_excluding_this_beat(p)
%             %             peaks(peak_indexes(p)) = 0;
%             %         end
%         end
%         [snr_max_imporvement, I_max_snr_imporvement] = max(snr_excluding_this_beat);
%         % % %     if snr >= snr_max_imporvement% && hr_std <= hr_std_excluding_this_beat(I_max_snr_imporvement)
%         % % %         continue
%         % % %         %     elseif data_filtered_env(peak_indexes(I_max_snr_imporvement)) < (peak_amps_max + peak_amps_min)/2 %if hr_std > hr_std_excluding_this_beat(I_max_snr_imporvement)
%         % % %     elseif data_filtered_env(peak_indexes(I_max_snr_imporvement)) < (peak_amps_max + peak_amps_med)/2 %if hr_std > hr_std_excluding_this_beat(I_max_snr_imporvement)
%         % % %         peaks(peak_indexes(I_max_snr_imporvement)) = 0;
%         % % %         peak_indexes = find(peaks);
%         % % %     end
%         if snr_max_imporvement > snr% && data_filtered_env(peak_indexes(I_max_snr_imporvement)) < (peak_amps_max + peak_amps_med)/2% && hr_std <= hr_std_excluding_this_beat(I_max_snr_imporvement)
%             peaks(peak_indexes(I_max_snr_imporvement)) = 0;
%             peak_indexes = find(peaks);
%             %     elseif data_filtered_env(peak_indexes(I_max_snr_imporvement)) < (peak_amps_max + peak_amps_min)/2 %if hr_std > hr_std_excluding_this_beat(I_max_snr_imporvement)
%             %     elseif data_filtered_env(peak_indexes(I_max_snr_imporvement)) < (peak_amps_max + peak_amps_med)/2 %if hr_std > hr_std_excluding_this_beat(I_max_snr_imporvement)
%             %
%         end
%         %     peak_indexes = find(peaks);
%     end
%     plot(peak_indexes, data_filtered_env(peak_indexes), 'ro', 'markersize', 16, 'linewidth', 3)
%
%     figure
%     hold on
%     plot((0:size(stacked_beats, 2)-1)/fs, stacked_beats');
%     plot((0:size(stacked_beats, 2)-1)/fs, mean(stacked_beats, 1), 'k', 'linewidth', 3);
%     grid
%
% end
%

