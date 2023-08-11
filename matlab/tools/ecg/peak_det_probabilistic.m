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

SigLen = size(signal, 2); % signal length

% Saturate the channels at k_sigma times the channel-wise STD
if isfield(params, 'saturate') && isequal(params.saturate, 1)
    if ~isfield(params, 'k_sigma')
        params.k_sigma = 12.0;
    end
    alpha = params.k_sigma * std(signal, [], 2);
    alpha = alpha(:);
    data_sat = diag(alpha) * tanh(diag(1./alpha) * signal);
else
    data_sat = signal;
end

% pass the channels through a narrow bandpass filter or a matched filter
if isfield(params, 'filter_type') && isequal(params.filter_type, 'MATCHED_FILTER')
    % Matched filter
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
        data_filtered(kk, :) = band_passsed(lag_matched : SigLen + lag_matched - 1);
    end
elseif isfield(params, 'filter_type') && isequal(params.filter_type, 'MULT_MATCHED_FILTER')
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
        band_passsed1 = band_passsed1(lag1 : SigLen + lag1 - 1);
        
        band_passsed2 = conv(template_matched2, data_sat(kk, :))/sum(template_matched2.^2);
        band_passsed2 = band_passsed2(lag2 : SigLen + lag2 - 1);
        
        data_filtered(kk, :) = sqrt(band_passsed1.^2 + band_passsed2.^2);
    end
elseif isfield(params, 'filter_type') && isequal(params.filter_type, 'MDMN')
    % Bandpass filter
    if ~isfield(params, 'wlen_md')
        params.wlen_md = 0.11;
    end
    if ~isfield(params, 'wlen_mn')
        params.wlen_mn = 0.01;
    end
    data_filtered = zeros(size(data_sat));
    wlen1 = round(params.wlen_md * fs);
    wlen2 = round(params.wlen_mn * fs);
    for kk = 1 : size(data_sat, 1)
        bl1 = BaseLine1(data_sat(kk, :), wlen1, 'md');
        bl2 = BaseLine1(bl1, wlen2, 'mn');
        data_filtered(kk, :) = data_sat(kk, :) - bl2;
        %         figure
        %         plot(data_sat(kk, :));
        %         hold on
        %         plot(bl1);
        %         plot(bl2);
        %         plot(data_filtered(kk, :));
        %         grid
        %         legend('datasat', 'bl1', 'bl2', 'filtered');
    end
elseif isfield(params, 'filter_type') && isequal(params.filter_type, 'WAVELET')
    % Wavelet denoiser
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
elseif ~isfield(params, 'filter_type') || (isfield(params, 'filter_type') && isequal(params.filter_type, 'BANDPASS_FILTER'))
    % Bandpass filter
    if ~isfield(params, 'low_cutoff')
        params.low_cutoff = 1.0;
    end
    if ~isfield(params, 'up_cutoff')
        params.up_cutoff = 40.0;
    end
    
    %     load('BPFilter15_40_300Hz.mat', 'h');
    %     data_filtered = filtfilt(h, 1, data_sat);
    data_lp = LPFilter(data_sat, params.up_cutoff/fs);
    data_filtered = data_lp - LPFilter(data_lp, params.low_cutoff/fs);
else
    error('Unknown preprocessing filter');
end

% calculate the power envelope of one or all channels
if ~isfield(params, 'wlen_power_env')
    params.wlen_power_env = 0.03;
end
wlen_power_env = ceil(params.wlen_power_env * fs);
data_filtered_env = filtfilt(ones(1, wlen_power_env), wlen_power_env, sqrt(mean(data_filtered.^2, 1)));

% pick the top percentage of the signal's power envelope for R-peak search
if ~isfield(params, 'hist_search_th')
    params.hist_search_th = 0.9;
end
bumps_amp_threshold = prctile(data_filtered_env, 100*params.hist_search_th);
bumps_indexes = find(data_filtered_env >= bumps_amp_threshold);

% search for the local peaks within a given sliding window length
if ~isfield(params, 'rpeak_search_wlen')
    params.rpeak_search_wlen = 0.03;
end
rpeak_search_half_wlen = floor(fs * params.rpeak_search_wlen / 2);
peak_indexes = [];
for jj = 1 : length(bumps_indexes)
    bump_index = bumps_indexes(jj);
    segment = max(1, bump_index - rpeak_search_half_wlen) : min(SigLen, bump_index + rpeak_search_half_wlen);
    if max(data_filtered_env(segment)) == data_filtered_env(bump_index)
        peak_indexes = cat(2, peak_indexes, bump_index);
    end
end
% peaks = zeros(1, SigLen);
% peaks(peak_indexes) = 1;


% bumps_amp_threshold = prctile(data_filtered_env(peak_indexes), 100*params.hist_search_th);
% peak_indexes_refined1 = peak_indexes(data_filtered_env(peak_indexes) >= bumps_amp_threshold);
%
% % remove excess beats based on waveform similarity
% event_width = 2 * round(median(diff(peak_indexes_refined1))/2) + 1;
% stacked_beats = EventStacker(data_filtered_env, peak_indexes_refined1, event_width);
% Rho_beats = corrcoef(stacked_beats');
% C_beats = (sum(Rho_beats, 1) - 1.0) / (length(peak_indexes_refined1) - 1); % Average correlation of each beat with the others
% I_omit = C_beats < min(0.9, 0.7 * max(C_beats));
% peak_indexes_refined2 = peak_indexes_refined1;
% peak_indexes_refined2(I_omit) = [];
% peaks = zeros(1, SigLen);
% peaks(peak_indexes_refined2) = 1;
%
% event_width = 2 * round(median(diff(peak_indexes_refined2))/2) + 1;
% stacked_beats = EventStacker(data_filtered_env, peak_indexes_refined2, event_width);
% % Rho_beats = stacked_beats * stacked_beats';
% Rho_beats = corrcoef(stacked_beats');
% C_beats = (sum(Rho_beats, 1) - 1.0) / (length(peak_indexes_refined2) - 1); % Average correlation of each beat with the others
% I_omit = C_beats < min(0.9, 0.7 * max(C_beats));
% peak_indexes_refined3 = peak_indexes_refined2;
% peak_indexes_refined3(I_omit) = [];
% peaks = zeros(1, SigLen);
% peaks(peak_indexes_refined3) = 1;
%
% event_width = 2 * round(median(diff(peak_indexes_refined3))/2) + 1;
% stacked_beats = EventStacker(data_filtered_env, peak_indexes_refined3, event_width);
% Rho_beats = stacked_beats * stacked_beats';
% C_beats = mean(Rho_beats, 1);
% I_omit = C_beats < 0.5 * max(C_beats);
% % Rho_beats = corrcoef(stacked_beats');
% % C_beats = (sum(Rho_beats, 1) - 1.0) / (length(peak_indexes_refined3) - 1); % Average correlation of each beat with the others
% % I_omit = C_beats < min(0.9, 0.85 * max(C_beats));
% peak_indexes_refined4 = peak_indexes_refined3;
% peak_indexes_refined4(I_omit) = [];
% peaks = zeros(1, SigLen);
% peaks(peak_indexes_refined4) = 1;

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

% % calculate the signal's power envelope histogram
% if ~isfield(params, 'n_bins')
%     params.n_bins = min(500, max(250, round(SigLen/10))); % somewhere between 250 and 500, if not specified
% end
% [data_filtered_env_hist.Values, data_filtered_env_hist.BinEdges] = histcounts(data_filtered_env, params.n_bins);
% data_filtered_env_PDF = data_filtered_env_hist.Values/sum(data_filtered_env_hist.Values); % calculate the PDF
% data_filtered_env_CDF = cumsum(data_filtered_env_PDF); % calculate the CDF

% bumps_threshold_index = find(data_filtered_env_CDF >= params.hist_search_th, 1, 'first');
% bumps_amp_threshold = data_filtered_env_hist.BinEdges(bumps_threshold_index);
% bumps_indexes = find(data_filtered_env >= bumps_amp_threshold);

% refine the R-peaks
if 0
    % matched filter using average beat shape
    event_width = 2 * round(median(diff(peak_indexes))/2) + 1;
    data_filtered_enhanced = zeros(size(data_filtered));
    for ch = 1 : size(data_filtered, 1)
        stacked_beats = EventStacker(data_filtered(ch, :), peak_indexes, event_width);
        robust_mean = RWAverage(stacked_beats);
        matched_filter_out = conv(robust_mean(end : -1 : 1), data_filtered(ch, :));
        lag = round(length(robust_mean)/2);
        data_filtered_enhanced(ch, :) = matched_filter_out(lag : SigLen + lag - 1);
    end
    data_filtered_env = sqrt(sum(data_filtered_enhanced.^2, 1));
end

% peaks = PeakDetection(mean(data_filtered, 1), 1.0/median(diff(peak_indexes)));
peaks = PeakDetection(data_filtered_env, 1.0/median(diff(peak_indexes)), 1);
peak_indexes = find(peaks);

% calculate a likelihood function for the R-peaks (useful for classification and scoring purposes)
if ~isfield(params, 'likelihood_sigma')
    %     params.likelihood_sigma = 0.01;
    params.likelihood_sigma = 0.3;
end
if ~isfield(params, 'max_likelihood_span')
    params.max_likelihood_span = 0.3;
end
tt = -params.max_likelihood_span/2 : 1 / fs : params.max_likelihood_span/2;
template = exp(-tt.^2/(2 * params.likelihood_sigma ^ 2));
lag = round(length(template)/2);
qrs_likelihood = conv(template, peaks);
qrs_likelihood = qrs_likelihood(lag : SigLen + lag - 1);

% Likelihood-based refinement of the R-peaks:
% peaks = PeakDetection(data_filtered_env .* qrs_likelihood, 1.0/median(diff(peak_indexes)));
% peaks = PeakDetection(mean(data_filtered, 1).*qrs_likelihood, 1.0/median(diff(peak_indexes)));
% peaks = PeakDetection(mean(signal, 1) .* qrs_likelihood, 1.0/median(diff(peak_indexes)));
% peaks = PeakDetection(mean(data_filtered_enhanced, 1) .* qrs_likelihood, 1.0/median(diff(peak_indexes)));
% peaks = PeakDetection(sign(skew(mean(data_filtered, 1))) * mean(data_filtered, 1).*qrs_likelihood, 1.0/median(diff(peak_indexes)), 1);

% search for local peaks with sign of most frequent among previously found peaks
data_filtered_mn_all_channels = mean(data_filtered, 1);
peaks = PeakDetection(mode(sign(data_filtered_mn_all_channels(peak_indexes))) * data_filtered_mn_all_channels .* qrs_likelihood, 1.0/median(diff(peak_indexes)), 1);

peak_indexes = find(peaks);

% remove excess beats based on waveform similarity
if ~isfield(params, 'RemoveUnsimilarBeats')
        params.RemoveUnsimilarBeats = true;
end
if params.RemoveUnsimilarBeats
    % Zero round (remove excess beats with extreme amplitude deviations from other peaks)
    peak_amps = data_filtered_mn_all_channels(peak_indexes);
    %     I_omit = abs(peak_amps - mean(peak_amps)) > 4.0 * std(peak_amps);
    I_omit = abs(peak_amps) < 0.5 * prctile(abs(peak_amps), 70.0);
    peak_indexes(I_omit) = [];
    
    % First round (negative correlations)
    event_width = 2 * round(median(diff(peak_indexes))/2) + 1;
    %     stacked_beats = EventStacker(data_filtered_env, peak_indexes, event_width);
    stacked_beats = EventStacker(data_filtered_mn_all_channels, peak_indexes, event_width);
    Rho_beats = corrcoef(stacked_beats');
    C_beats = (sum(Rho_beats, 1) - 1.0) / (length(peak_indexes) - 1); % Average correlation of each beat with the others
    I_omit = C_beats < 0;
    peak_indexes(I_omit) = [];
    
    % Second round (omit beats with low variance)
    peak_amps = data_filtered_mn_all_channels(peak_indexes);
    event_width = 2 * round(median(diff(peak_indexes))/2) + 1;
    %     stacked_beats = EventStacker(data_filtered_env, peak_indexes, event_width);
    [stacked_beats, num_non_zeros] = EventStacker(data_filtered_mn_all_channels, peak_indexes, event_width);
    std_beats = (event_width -1) * std(stacked_beats, 0, 2)' ./ (num_non_zeros - 1); % compensates the boundary events STD that were zero-padded in EventStacker
    I_omit = std_beats < 0.2*median(std_beats) & abs(peak_amps - mean(peak_amps)) > 2.0 * std(peak_amps);
    if(length(I_omit) < peak_indexes) % don't remove all peaks!
        peak_indexes(I_omit) = [];
    end
    
    % Third round (omit beats with low correlations)
    event_width = 2 * round(median(diff(peak_indexes))/2) + 1;
    %     stacked_beats = EventStacker(data_filtered_env, peak_indexes, event_width);
    stacked_beats = EventStacker(data_filtered_mn_all_channels, peak_indexes, event_width);
    Rho_beats = corrcoef(stacked_beats');
    C_beats = (sum(Rho_beats, 1) - 1.0) / (length(peak_indexes) - 1); % Average correlation of each beat with the others
    %     I_omit = C_beats < max(0.3, min(0.8, 0.5 * max(C_beats))); % omit beats with low correlations
    I_omit = abs(C_beats - mean(C_beats)) > 3.0 * std(C_beats) & C_beats < 0.5; % omit beats with low correlations
    if(length(I_omit) < peak_indexes) % don't remove all peaks!
        peak_indexes(I_omit) = [];
    end
    
    % Fourth round (omit beats that increase average beat SNR and reduce HR variance. not applicable to the first and last beats)
    max_itr = length(peak_indexes);
    for itr = 1 : max_itr
        NumBeats = length(peak_indexes);
        event_width = 2 * round(median(diff(peak_indexes))/2) + 1;
        stacked_beats = EventStacker(data_filtered_mn_all_channels, peak_indexes, event_width);
        
        %     ECG_robust_mean = RWAverage(stacked_beats);
        ECG_robust_mean = mean(stacked_beats, 1);
        ECG_robust_mean_replicated = ones(NumBeats, 1) * ECG_robust_mean;
        noise = stacked_beats - ECG_robust_mean_replicated;
        %snr = 10*log10(trace(ECG_robust_mean_replicated*ECG_robust_mean_replicated')/trace(noise*noise'));
        snr = 20*log10(norm(ECG_robust_mean_replicated, 'fro')/norm(noise,'fro'));
        %     hr_std = max(abs(diff(peak_indexes)));
        %         hr_std = std(diff([1, peak_indexes, SigLen]));
        hr_std = std(diff(peak_indexes));
        snr_excluding_this_beat = zeros(1, NumBeats);
        hr_std_excluding_this_beat = zeros(1, NumBeats);
        for p = 1 : NumBeats
            all_but_this_beat_index = 1 : NumBeats;
            all_but_this_beat_index(p) = [];
            %         ECG_robust_mean = RWAverage(stacked_beats(all_but_this_beat_index, :));
            ECG_robust_mean = mean(stacked_beats(all_but_this_beat_index, :), 1);
            ECG_robust_mean_replicated = ones(NumBeats - 1, 1) * ECG_robust_mean;
            % estimate channel quality SNR
            noise = stacked_beats(all_but_this_beat_index, :) - ECG_robust_mean_replicated;
            %snr_excluding_this_beat(p) = 10*log10(trace(ECG_robust_mean_replicated*ECG_robust_mean_replicated')/trace(noise*noise'));
            snr_excluding_this_beat(p) = 20*log10(norm(ECG_robust_mean_replicated, 'fro')/norm(noise,'fro'));
            %         hr_std_excluding_this_beat(p) = max(abs(diff(peak_indexes(all_but_this_beat_index))));
            %             hr_std_excluding_this_beat(p) = std(diff([1, peak_indexes(all_but_this_beat_index), SigLen]));
            hr_std_excluding_this_beat(p) = std(diff(peak_indexes(all_but_this_beat_index)));
        end
        [snr_max_imporvement, I_max_snr_imporvement] = max(snr_excluding_this_beat);
        if snr_max_imporvement > snr && hr_std > hr_std_excluding_this_beat(I_max_snr_imporvement) && I_max_snr_imporvement > 1 && I_max_snr_imporvement < NumBeats % && data_filtered_env(peak_indexes(I_max_snr_imporvement)) < (peak_amps_max + peak_amps_med)/2% && hr_std <= hr_std_excluding_this_beat(I_max_snr_imporvement)
            peak_indexes(I_max_snr_imporvement) = [];
        else
            break;
        end
    end
end

peaks = zeros(1, SigLen);
peaks(peak_indexes) = 1;
lag = round(length(template)/2);
qrs_likelihood = conv(template, peaks);
qrs_likelihood = qrs_likelihood(lag : SigLen + lag - 1);

if 0 % peak refinement based on likelihoods
    signal_abs = sqrt(mean(data_filtered.^2, 1));
    signal_likelihood = qrs_likelihood .* signal_abs/max(signal_abs);
    bumps_indexes2 = find(signal_likelihood >= 0.4);
    peak_indexes = [];
    for jj = 1 : length(bumps_indexes2)
        index = bumps_indexes2(jj);
        segment = max(1, index - rpeak_search_half_wlen) : min(SigLen, index + rpeak_search_half_wlen);
        if max(signal_likelihood(segment)) == signal_likelihood(index)
            peak_indexes = cat(2, peak_indexes, index);
        end
    end
    peaks = zeros(1, SigLen);
    peaks(peak_indexes) = 1;
end

if 0
    test plot results
    time = (0 : SigLen - 1) / fs;
    figure
    scatter(time, signal', 28, qrs_likelihood(:)*[1 0 0], 'filled');
    hold on
    plot(time, data_filtered_env);
    plot(time, signal_likelihood);
    grid
end

if 0
    f0 = 1/params.rpeak_search_wlen;
    peaks = PeakDetection(data_filtered_env, f0/fs);
    peak_indexes = find(peaks);
    peak_amps = data_filtered_env(peak_indexes);
    peak_amps_max = max(peak_amps);
    peak_amps_min = min(peak_amps);
    peak_amps_med = median(peak_amps);
    
    figure
    hold on
    grid
    plot(signal)
    plot(data_sat)
    plot(data_filtered)
    plot(data_filtered_env)
    plot(peak_indexes, data_filtered_env(peak_indexes), 'ko', 'markersize', 14, 'linewidth', 3)
    for itr = 1 : 150
        NumBeats = length(peak_indexes);
        event_width = 2 * round(0.75    *           median(diff(peak_indexes))/2) + 1;
        %     event_width = 2 * round(fs * params.rpeak_search_wlen/2) + 1;
        
        
        
        %     stacked_beats = EventStacker(data_filtered, peak_indexes, event_width);
        stacked_beats = EventStacker(data_filtered_env, peak_indexes, event_width);
        
        
        
        %     ECG_robust_mean = RWAverage(stacked_beats);
        ECG_robust_mean = mean(stacked_beats, 1);
        ECG_robust_mean_replicated = ones(NumBeats, 1) * ECG_robust_mean;
        noise = stacked_beats - ECG_robust_mean_replicated;
        %snr = 10*log10(trace(ECG_robust_mean_replicated*ECG_robust_mean_replicated')/trace(noise*noise'));
        snr = 20*log10(norm(ECG_robust_mean_replicated, 'fro')/norm(noise,'fro'));
        
        %     hr_std = max(abs(diff(peak_indexes)));
        hr_std = std(diff([1, peak_indexes, SigLen]));
        snr_excluding_this_beat = zeros(1, NumBeats);
        hr_std_excluding_this_beat = zeros(1, NumBeats);
        for p = 1 : NumBeats
            all_but_this_beat_index = 1 : NumBeats;
            all_but_this_beat_index(p) = [];
            %         ECG_robust_mean = RWAverage(stacked_beats(all_but_this_beat_index, :));
            ECG_robust_mean = mean(stacked_beats(all_but_this_beat_index, :), 1);
            ECG_robust_mean_replicated = ones(NumBeats - 1, 1) * ECG_robust_mean;
            % estimate channel quality SNR
            noise = stacked_beats(all_but_this_beat_index, :) - ECG_robust_mean_replicated;
            %snr_excluding_this_beat(p) = 10*log10(trace(ECG_robust_mean_replicated*ECG_robust_mean_replicated')/trace(noise*noise'));
            snr_excluding_this_beat(p) = 20*log10(norm(ECG_robust_mean_replicated, 'fro')/norm(noise,'fro'));

            %         hr_std_excluding_this_beat(p) = max(abs(diff(peak_indexes(all_but_this_beat_index))));
            hr_std_excluding_this_beat(p) = std(diff([1, peak_indexes(all_but_this_beat_index), SigLen]));
            %         if snr < snr_excluding_this_beat(p)
            %             peaks(peak_indexes(p)) = 0;
            %         end
        end
        [snr_max_imporvement, I_max_snr_imporvement] = max(snr_excluding_this_beat);
        % % %     if snr >= snr_max_imporvement% && hr_std <= hr_std_excluding_this_beat(I_max_snr_imporvement)
        % % %         continue
        % % %         %     elseif data_filtered_env(peak_indexes(I_max_snr_imporvement)) < (peak_amps_max + peak_amps_min)/2 %if hr_std > hr_std_excluding_this_beat(I_max_snr_imporvement)
        % % %     elseif data_filtered_env(peak_indexes(I_max_snr_imporvement)) < (peak_amps_max + peak_amps_med)/2 %if hr_std > hr_std_excluding_this_beat(I_max_snr_imporvement)
        % % %         peaks(peak_indexes(I_max_snr_imporvement)) = 0;
        % % %         peak_indexes = find(peaks);
        % % %     end
        if snr_max_imporvement > snr% && data_filtered_env(peak_indexes(I_max_snr_imporvement)) < (peak_amps_max + peak_amps_med)/2% && hr_std <= hr_std_excluding_this_beat(I_max_snr_imporvement)
            peaks(peak_indexes(I_max_snr_imporvement)) = 0;
            peak_indexes = find(peaks);
            %     elseif data_filtered_env(peak_indexes(I_max_snr_imporvement)) < (peak_amps_max + peak_amps_min)/2 %if hr_std > hr_std_excluding_this_beat(I_max_snr_imporvement)
            %     elseif data_filtered_env(peak_indexes(I_max_snr_imporvement)) < (peak_amps_max + peak_amps_med)/2 %if hr_std > hr_std_excluding_this_beat(I_max_snr_imporvement)
            %
        end
        %     peak_indexes = find(peaks);
    end
    plot(peak_indexes, data_filtered_env(peak_indexes), 'ro', 'markersize', 16, 'linewidth', 3)
    
    figure
    hold on
    plot((0:size(stacked_beats, 2)-1)/fs, stacked_beats');
    plot((0:size(stacked_beats, 2)-1)/fs, mean(stacked_beats, 1), 'k', 'linewidth', 3);
    grid
    
end

