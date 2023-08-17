function [x_denoised, baseline_reference, x_mains_cancelled] = ecg_strip_viewer_waterfall(x_raw, fs, f_mains, ttl, varargin)
% 
% ecg_strip_viewer_waterfall: A function to process and visualize ECG signals.
%
% USAGE:
%   [xd, bl, xm] = ecg_strip_viewer_waterfall(x, fs, f_mains, title, params);
%
% INPUTS:
%   x_raw: Raw ECG signal(s). Can be a single-channel or multi-channel signal.
%   fs: Sampling frequency of the ECG signal in Hz.
%   f_mains: Mains interference frequency (typically 50 or 60 Hz).
%   title: Title string for the visualization and naming of outputs.
%   params: Optional input, contains the parameters for processing.
%
% OUTPUTS:
%   x_denoised: Denoised ECG signal.
%   baseline_reference: Extracted baseline from the ECG signal.
%   x_mains_cancelled: ECG signal after mains interference cancellation.
%
% EXAMPLE:
%   [xd, bl, xm] = ecg_strip_viewer_waterfall(x, 500, 50, 'Test ECG');
%
% PARAMETER DESCRIPTION: The last input params is a structure with these
% fields:
%   - plot_results: Boolean, decides if results should be plotted. (Default: 1)
%   - save_output_files: Boolean, decides if results should be saved. (Default: 0)
%   - save_output_images: Boolean, decides if images should be saved. (Default: 0)
%   - aply_source_separation: Boolean, decides if Blind Source Separation (BSS) should be applied. (Default: 0)
%   - build_lead_III: Boolean, to build lead III from leads I and II. (Default: 0)
%   - correct_lead_I_polarity: Boolean, to correct the polarity of lead I. (Default: 0)
%   - file_out_name: String, naming for the output files. (Default: [ttl , '_processed'])
%   - image_name: String, naming for the output images. (Default: [ttl , '_processed'])
%   - Q_factor: Numeric, Q-factor for the mains notch filter. (Default: 40)
%   - w1, w2, w3, w4: Numerics, defining window sizes in seconds for different stages of baseline wander removal. (Defaults: 0.29, 0.32, 0.35, 0.2 respectively)
%   - baseline_removal_method: String, method for removing baseline wander. (Default: 'MDMDMDMN')
%   - mains_notch_method: String, method for removing mains interference. (Default: 'IIRNOTCH')
%   - remove_2nd_mains_harmonic: Boolean, decides if the second harmonic of the mains frequency should be removed. (Default: 0)
%   - denoising_method: String, method for denoising the ECG signal. (Default: 'BYPASS')
%   - DiffOrder: Integer, smoothness constraint order for denoising. (Default: 2)
%   - beat_avg_method: String, method for beat averaging. (Default: 'TIME')
%
%   Revision History:
%       2023: First release.
% 
%   Reza Sameni, 2006-2023
%   The Open-Source Electrophysiological Toolbox
%   https://github.com/alphanumericslab/OSET

if nargin > 4 && ~isempty(varargin{1})
    params_in = varargin{1};
else
    params_in = [];
end

if isfield(params_in, 'plot_results')
    params.plot_results = params_in.plot_results;
else
    params.plot_results = 1; % plot the results or not
end

if isfield(params_in, 'save_output_files')
    params.save_output_files = params_in.save_output_files;
else
    params.save_output_files = 0; % save the results or not
end

if isfield(params_in, 'save_output_images')
    params.save_output_images = params_in.save_output_images;
else
    params.save_output_images = 0; % save the images or not
end

if isfield(params_in, 'aply_source_separation')
    params.aply_source_separation = params_in.aply_source_separation;
else
    params.aply_source_separation = 0; % Apply BSS or not
end

if isfield(params_in, 'build_lead_III')
    params.build_lead_III = params_in.build_lead_III;
else
    params.build_lead_III = 0;
end

if isfield(params_in, 'correct_lead_I_polarity')
    params.correct_lead_I_polarity = params_in.correct_lead_I_polarity;
else
    params.correct_lead_I_polarity = 0;
end

if isfield(params_in, 'file_out_name')
    params.file_out_name = params_in.file_out_name;
else
    params.file_out_name = [ttl , '_porcessed'];
end

if isfield(params_in, 'image_name')
    params.image_name = params_in.image_name;
else
    params.image_name = [ttl , '_porcessed'];
end

if isfield(params_in, 'w1')
    params.w1 = params_in.w1;
else
    params.w1 = 0.29; % First stage baseline wander removal window size in seconds
    % params.w1 = 0.25; % First stage baseline wander removal window size in seconds
end

if isfield(params_in, 'w2')
    params.w2 = params_in.w1;
else
    params.w2 = 0.32; % Second stage baseline wander removal window size in seconds
    % params.w2 = 0.52; % Second stage baseline wander removal window size in seconds
end

if isfield(params_in, 'w3')
    params.w3 = params_in.w1;
else
    params.w3 = 0.35; % Third stage baseline wander removal window size in seconds
    % params.w3 = 0.54; % Third stage baseline wander removal window size in seconds
end

if isfield(params_in, 'w4')
    params.w4 = params_in.w1;
else
    params.w4 = 0.2; % Fourth stage baseline wander removal window size in seconds
    % params.w4 = 0.49; % Fourth stage baseline wander removal window size in seconds
end

% Set baseline wander removal method:
if isfield(params_in, 'baseline_removal_method')
    params.baseline_removal_method = params_in.baseline_removal_method;
else
    params.baseline_removal_method = 'MDMDMDMN'; % 'BYPASS', 'MNMN', 'MDMN', 'TIKHONOV', 'MDMDMDMN'
end

% Set mains removal method:
if isfield(params_in, 'mains_notch_method')
    params.mains_notch_method = params_in.mains_notch_method;
else
    params.mains_notch_method = 'IIRNOTCH';% 'BYPASS', 'IIRNOTCH', 'KFNOTCH', 'KFNOTCH2'
end

if isfield(params_in, 'Q_factor')
    params.Q_factor = params_in.Q_factor;
else
    params.Q_factor = 40; % Q-factor of the mains notch filter
end

if isfield(params_in, 'remove_2nd_mains_harmonic')
    params.remove_2nd_mains_harmonic = params_in.remove_2nd_mains_harmonic;
else
    params.remove_2nd_mains_harmonic = 0;
end

if isfield(params_in, 'denoising_method')
    params.denoising_method = params_in.denoising_method;
else
    params.denoising_method = 'BYPASS';%'BYPASS';%'SEGMENTWISETIKHONOV'; %'TIKHONOV';%'WAVELET';
end

if isfield(params_in, 'DiffOrder')
    params.DiffOrder = params_in.DiffOrder;
else
    params.DiffOrder = 2; % smoothness constraint order greater than 1
    %     params.DiffImpulseResponse = [1 2 1]; % smoothness filter impulse (can replace DiffOrder in some modes of TikhonovRegularization)
end

if isfield(params_in, 'beat_avg_method')
    params.beat_avg_method = params_in.beat_avg_method;
else
    params.beat_avg_method = 'TIME';
end


% NOTE: Correct the dimensions of single channel signals
if size(x_raw, 2) < size(x_raw, 1)
    x_raw = x_raw';
end

% Only porcess the first two independent channels (the other data are built from the first two channels)
alpha = 1;
if size(x_raw, 1) == 2 && params.build_lead_III
    if params.correct_lead_I_polarity % Find Lead I polarity
        alpha = sign(sum(x_raw(1,:).*x_raw(2,:)));
    end
    lead3_corrected = x_raw(2,:) - alpha*x_raw(1,:);
    x_raw = cat(1, x_raw, lead3_corrected);
end
% select baseline wander removal method
BLFilterParams.wlen1 = round(params.w1 * fs);
BLFilterParams.wlen2 = round(params.w2 * fs);
BLFilterParams.wlen3 = round(params.w3 * fs);
BLFilterParams.wlen4 = round(params.w4 * fs);
BLFilterParams.DiffOrder = 2; % smoothness constraint order greater than 1
BLFilterParams.lambda = (4.8e3)^BLFilterParams.DiffOrder; % The smoothness parameter for Tikhonov regularization (higher values penalize the roughness of the signal)
baseline_reference = baseline_filter(x_raw, params.baseline_removal_method, BLFilterParams);

% select the mains removal method
switch params.mains_notch_method
    case 'BYPASS' % bypass the notch filter
        x_mains_cancelled = x_raw;
    case 'IIRNOTCH' % second-order IIR notch filter
        W1 = f_mains/(fs/2);
        BW1 = W1/params.Q_factor;
        [b1,a1] = iirnotch(W1,BW1);
        x_mains_cancelled = filtfilt(b1, a1, x_raw')';
        if params.remove_2nd_mains_harmonic % remove second harmonic
            W2 = (2.0*f_mains)/(fs/2);
            BW2 = W2/Q_factor;
            [b2,a2] = iirnotch(W2,BW2);
            x_mains_cancelled = filtfilt(b2, a2, x_mains_cancelled')';
        end
    case 'KFNOTCH' % standard KF implementation
        W1 = f_mains/(fs/2);
        BW1 = W1/Q_factor;
        [b1,a1] = iirnotch(W1,BW1);
        x_mains_cancelled = filtfilt(b1, a1, x_raw')';

        Q = 0.01*var(x_raw);
        R = var(x_raw - x_mains_cancelled);
        gamma = 0.99;
        x_mains_cancelled = zeros(size(x_raw));
        for ch = 1 : size(x_raw, 1)
            [~, x_mains_cancelled(ch, :)] = KFNotch(x_raw(ch, :), f_mains, fs, Q, R, gamma);
        end
    case 'KFNOTCH2' % simplified KF implementation (works with the Kalman filter Q/R ratio)
        W1 = f_mains/(fs/2);
        BW1 = W1/Q_factor;
        [b1,a1] = iirnotch(W1,BW1);
        x_mains_cancelled = filtfilt(b1, a1, x_raw')';

        Q = 0.01*var(x_raw);
        R = var(x_raw - x_mains_cancelled);
        wlen = 10;
        x_mains_cancelled = zeros(size(x_raw));
        for ch = 1 : size(x_raw, 1)
            [~, x_mains_cancelled(ch, :)] = KFNotch2(x_raw(ch, :), f_mains, fs, 1000*Q/R, wlen);
        end
    otherwise
        error('Notch filter option undefined');
end

% Apply source separation
if params.aply_source_separation
    x_diff = x_mains_cancelled(1, :) - x_mains_cancelled(2, :);
    x_aug = [x_mains_cancelled ; x_diff ; baseline_reference];

    W = real(jadeR(x_aug, 2));
    s = W * x_aug;
    baseline_s = baseline_removal_method(s, baseline_removal_method, BLFilterParams);
    s_den = s - baseline_s;
    x_common_mode_den = real(pinv(W) * s_den);
    x_common_mode_den = x_common_mode_den(1 : 3, :); % Just pick up the first two original channels plus their difference
    if params.plot_bss
        plot_multichannel_data(x_aug, size(x_aug, 1), 'b', fs, 'BSS stage signal xaug');
        plot_multichannel_data(s, size(s, 1), 'r', fs, 'BSS stage source signal s');
        plot_multichannel_data(baseline_s, size(baseline_s, 1), 'g', fs, 'BSS stage baselines');
        plot_multichannel_data(s_den, size(s_den, 1), 'm', fs, 'BSS stage baseline removed');
        plot_multichannel_data(x_common_mode_den, size(x_den, 1), 'k', fs, 'BSS stage x-common-mode-den');
    end
else
    x_common_mode_den = x_mains_cancelled - baseline_reference;
end


% select the mains removal method
switch params.denoising_method
    case 'BYPASS' % bypass the denoiser
        x_denoised = x_common_mode_den;
    case 'TIKHONOV'
        % A Tikhonov regularizer baseline detector. See TikhonovRegularization for parameter choices
        lambda = (9.47)^params.DiffOrder; % The smoothness parameter for Tikhonov regularization (higher values penalize the roughness of the signal)
        x_denoised = TikhonovRegularization(x_common_mode_den, params.DiffOrder, lambda);
    case 'SEGMENTWISETIKHONOV'
        % Segment-wise Tikhonov regularizer. See ECGSmoothnessPriorsDenoiserBW for help and references
        wlen = 15e-3;%100e-3; % window length (s)
        % withwindow = 0; % with(1) or without(0) windowing
        mode = 6; % 0..7
        nvar = 0.8e-3;
        [~, x_denoised] = ECGSmoothnessPriorsDenoiserBW(x_common_mode_den, nvar, mode, params.DiffOrder, round(wlen*fs), 1e-8, 250, 1, 500);
    case 'WAVELET'
        % A wavelet denoiser
        x_denoised = zeros(size(x_common_mode_den));
        for ch = 1 : size(x_common_mode_den, 1)
            x_denoised(ch, :) = wden(x_common_mode_den(ch, :), 'rigrsure', 's', 'sln', 4, 'sym5'); % use SLN since we seek beat-wise performance
        end
    case 'GP_FILTER'
        % % %     % Apply MAP filter
        peaks = peak_det_local_search(x_common_mode_den(1, :), 1.0/fs, [], 3);
        GPfilterparams.bins = 300; % number of phase domain bins
        GPfilterparams.BEAT_AVG_METHOD = 'MEDIAN'; % 'MEAN' or 'MEDIAN'
        GPfilterparams.NOISE_VAR_EST_METHOD = 'AVGLOWER'; %'MIN', 'AVGLOWER', 'MEDLOWER', 'PERCENTILE'
        GPfilterparams.p = 0.5;
        GPfilterparams.avg_bins = 10;
        GPfilterparams.SMOOTH_PHASE = 'GAUSSIAN';
        GPfilterparams.gaussianstd = 1.0;
        GPfilterparams.plotresults = 0;
        GPfilterparams.nvar_factor = 3.0; % noise variance over/under estimation factor (1 by default)
        %     [x_denoised, data_prior_est_phase_based, n_var_est_phase_based] = ECGPhaseDomainMAPFilter(x_common_mode_den(max_snr_index, :), peaks, GPfilterparams);
        x_denoised = ECGTimeDomainMAPFilter(x_common_mode_den, peaks, GPfilterparams);

    otherwise
        error('Denoising filter option undefined');
end

% Detect the R-peaks
SignalLen = size(x_denoised, 2);

peaks_all_channels = zeros(size(x_denoised));
snr = -inf(1, size(x_denoised, 1));
for ch = 1 : size(x_denoised, 1)
    peak_detector_params.saturate = 1;
    peak_detector_params.k_sigma = 4;
    peak_detector_params.hist_search_th = 0.9;
    peak_detector_params.rpeak_search_wlen = 0.4; % MAX detectable HR (in BPM) = 60/rpeak_search_wlen
    peak_detector_params.filter_type = 'MDMN';%'MULT_MATCHED_FILTER';%'BANDPASS_FILTER', 'MATCHED_FILTER', 'MULT_MATCHED_FILTER', 'MDMN', 'WAVELET'
    [peaks_all_channels(ch, :), peak_indexes, ~] = peak_det_probabilistic(x_denoised(ch, :), fs, peak_detector_params); % qrs_likelihood
    NumBeats = length(peak_indexes);
    if NumBeats == 0
        continue;
    end
    event_width = round(1.2 * median(diff([1, peak_indexes, SignalLen]))); % Add the first and last indexes
    if(mod(event_width, 2) == 0)
        event_width = event_width + 1;
    end

    % time domain average beats
    stacked_beats = event_stacker(x_denoised(ch, :), peak_indexes, event_width);
    %         ECG_mean = mean(stacked_beats, 1);
    %         ECG_median = median(stacked_beats, 1);
    %         [ECG_robust_mean, ECG_var_mn, ECG_robust_median, ECG_var_md] = RWAverage(stacked_beats);
    ECG_robust_mean = robust_weighted_average(stacked_beats);
    ECG_robust_mean_replicated = ones(NumBeats, 1) * ECG_robust_mean;

    % estimate channel quality SNR
    noise = stacked_beats - ECG_robust_mean_replicated;
    snr(ch) = 10*log10(trace(ECG_robust_mean_replicated*ECG_robust_mean_replicated')/trace(noise*noise'));
end

[max_snr, max_snr_index] = max(snr);
peaks = peaks_all_channels(max_snr_index, :);
peak_indexes = find(peaks);
NumBeats = length(peak_indexes);
event_width = round(1.2 * median(diff([1, peak_indexes, SignalLen]))); % Add the first and last indexes
if(mod(event_width, 2) == 0)
    event_width = event_width + 1;
end

switch params.beat_avg_method
    case 'TIME'
        % time domain average beats
        stacked_beats = event_stacker(x_denoised(max_snr_index, :), peak_indexes, event_width);
        ECG_mean = mean(stacked_beats, 1);
        ECG_median = median(stacked_beats, 1);
        [ECG_robust_mean, ECG_var_mn, ECG_robust_median, ECG_var_md] = robust_weighted_average(stacked_beats);

    case 'PHASE'
        % ECG_robust_mean_replicated = ones(NumBeats, 1) * ECG_robust_mean;
        % ECG phase
        [phase, ~] = phase_calculator(peaks);     % phase calculation
        teta = 0;                                       % desired phase shift
        pphase = PhaseShifting(phase,teta);             % phase shifting
        % phase domain average beats
        bins = round(0.95 * median(diff([1, peak_indexes, SignalLen])));                                     % number of phase bins
        [ECG_robust_mean, ECG_std_phase, ~, ECG_robust_median] = MeanECGExtraction(x_denoised(max_snr_index, :), pphase, bins, 0); % mean ECG extraction
        ECG_var_mn = ECG_std_phase .^2;
        ECG_var_md = ECG_var_mn;
    otherwise
        error('Undefined beat averaging method');
end

% Plot the signals
if params.plot_results
    tm = 1000.0 * ((-event_width/2 : event_width/2 - 1) + 0.5)/fs;

    figure('Position', [100, 100, 1250, 900],'PaperUnits', 'points', 'PaperOrientation', 'landscape');%, 'visible','off');
    %     h3 = subplot(421);
    subplot('Position', [0.05, 0.63, 0.92, 0.29]);
    lbl = {};
    t = (0 : length(x_mains_cancelled) - 1)/fs;
    for ch = 1 : size(x_denoised, 1)
        offset = (ch - 1) * (2.0);

        plot(t, x_raw(ch, :) - offset, 'linewidth', 1); lbl = cat(2, lbl, {'Raw signal'});
        hold on
        plot(t, baseline_reference(ch, :) - offset, 'linewidth', 2); lbl = cat(2, lbl, {'Baseline'});
        %         plot(t(peak_indexes), x_raw(ch, peak_indexes) + offset, 'co', 'linewidth', 3, 'markersize', 12); lbl = cat(2, lbl, {'R-peaks'});
        mid_time_index = round(0.95 * length(t));
        text(t(mid_time_index), baseline_reference(ch, mid_time_index) - offset - 0.3, ['Lead ', num2str(ch)], 'fontsize', 14);
    end
    grid
    ylabel('Amplitude (mv)');
    %             xlabel('time (s)');
    set(gca, 'fontsize', 14);
    hold off
    %legend(lbl, 'interpreter', 'none', 'Location', 'Best', 'Orientation', 'horizontal');
    %title(['Baseline removal for subject:' subjectcode, ' Ch #' num2str(ch)])
    if size(x_denoised, 1) > 1
        title([ttl, ', highest SNR lead: #', num2str(max_snr_index), ', SNR=', num2str(snr(max_snr_index), 3), 'dB'], 'fontsize', 14);
        if alpha < 0
            note = ', Note: Lead 1 inverted';
        else
            note = ' ';
        end
        subtitle(['Channel SNRs (dB)= ', mat2str(snr, 3), note], 'fontsize', 14);
    else
        title([ttl, ', lead #', num2str(max_snr_index), ', SNR=', num2str(snr(max_snr_index), 3), 'dB'], 'fontsize', 14);
    end
    %     set(gca,'Xticklabel',[])
    %     h3.Position = [h3.Position(1), max(0, h3.Position(2) - 0.08), h3.Position(3) + 0.5, h3.Position(4) + 0.08];
    %     h3.Position = [h3.Position(1), h3.Position(2) - 0.06, h3.Position(3) + 0.5, h3.Position(4) + 0.08];
    aa = axis;
    aa(1) = t(1);
    aa(2) = t(end);
    axis(aa);
    %     h3.Position = [0.05, 0.52, 0.92, 0.4];
    %     h3_copy = h3.Position;

    %     h4 = subplot(423);
    subplot('Position', [0.05, 0.48, 0.92, 0.10]);
    plot(t, x_denoised(max_snr_index, :), 'linewidth', 1);
    hold on
    % % %     plot(t, x_denoised_map(max_snr_index, :), 'linewidth', 1);
    plot(t(peak_indexes), x_denoised(max_snr_index, peak_indexes), 'ro', 'linewidth', 1, 'markersize', 14);

    %     legend({['Highest SNR lead: ', num2str(max_snr_index)], 'R-peaks'}, 'Orientation', 'horizontal', 'Location', 'southeast')
    %     title(h4, 'Heart rate');
    title(['Highest SNR lead #', num2str(max_snr_index), ' and R-peaks'], 'fontsize', 14)

    ylabel('Amplitude (mv)');
    %     xlabel('time (s)');
    set(gca, 'fontsize', 14);
    set(gca,'Xticklabel',[])
    grid
    aa = axis;
    aa(1) = t(1);
    aa(2) = t(end);
    %             aa(1) = t(peak_indexes(2));
    %             aa(2) = t(peak_indexes(end));
    axis(aa);
    %     h4.Position = [h4.Position(1), h4.Position(2), h3_copy(3), max(0, h4.Position(4) - 0.08)];
    %     h4.Position = [h4.Position(1), h4.Position(2)- 0.06, h4.Position(3) + 0.5, h4.Position(4)-0.03];
    %     h4.Position = [0.05, 0.32, 0.92, 0.2];

    %     h5 = subplot(425);
    subplot('Position', [0.05, 0.36, 0.92, 0.10]);
    %     h5.Position = [h5.Position(1), h5.Position(2), h3_copy(3), max(0, h5.Position(4) - 0.08)];
    %     h5.Position = [h5.Position(1), h5.Position(2), h5.Position(3) + 0.5, h5.Position(4)-0.07];
    plot(t(peak_indexes(2:end)), 60.0*fs./diff(peak_indexes), 'linewidth', 1);
    hold on
    plot(t(peak_indexes(2:end)), 60.0*fs./diff(peak_indexes), 'ro', 'linewidth', 1);
    %     title(h4, 'Heart rate');
    ylabel('HR (BPM)');
    xlabel('time (s)');
    set(gca, 'fontsize', 14);
    grid
    aa = axis;
    aa(1) = t(1);
    aa(2) = t(end);
    %             aa(1) = t(peak_indexes(2));
    %             aa(2) = t(peak_indexes(end));
    axis(aa);

    %     h1 = subplot(427);
    subplot('Position', [0.05, 0.05, 0.41, 0.27]);
    lgnd = {};
    plot(tm, ECG_robust_median, 'linewidth', 3); lgnd = cat(2, lgnd, 'Robust Median Beat');
    hold on
    %plot(tm, ECG_mean, 'linewidth', 3); lgnd = cat(2, lgnd, 'Mean Beat');
    plot(tm, ECG_median, 'linewidth', 3); lgnd = cat(2, lgnd, 'Median Beat');
    %plot(tm, ECG_robust_mean, 'linewidth', 3); lgnd = cat(2, lgnd, 'Robust Mean Beat');
    plot(tm, stacked_beats', 'color', 0.7*ones(1, 3)); lgnd = cat(2, lgnd, 'All beats');
    legend(lgnd, 'Location', 'Best');
    grid
    chi = get(gca, 'Children'); %Returns handles to the patch and line objects
    set(gca, 'Children',flipud(chi)); %Reverse the stacking order so that the averages overlay the beats
    %     title(h1, 'Time-domain averaging');
    xlabel('time (ms)');
    ylabel('Amplitude (mV)');
    set(gca, 'fontsize', 14)
    set(gca, 'box', 'on')
    aa = axis;
    aa(1) = tm(1);
    aa(2) = tm(end);
    aa(3) = 1.1 * min([ECG_robust_median, ECG_median]);
    aa(4) = 1.1 * max([ECG_robust_median, ECG_median]);
    axis(aa);

    %     h2 = subplot(428);
    subplot('Position', [0.55, 0.05, 0.41, 0.27]);
    if false
        lgnd = {};
        tmp_copy = x_denoised(ch, :);
        pphase_copy = pphase;
        phase_jumps = abs(diff(pphase)) > pi;
        tmp_copy(phase_jumps) = nan;
        pphase_copy(phase_jumps) = nan;
        plot(pphase_copy, tmp_copy, 'color', 0.7*ones(1, 3)); lgnd = cat(2, lgnd, 'All beats');
        %             plot(pphase, x_denoised(ch, :), '.', 'color', 0.7*ones(1, 3)); lgnd = cat(2, lgnd, 'All beats');
        hold on
        plot(mean_phase, ECG_mean_phase, 'linewidth', 3); lgnd = cat(2, lgnd, 'Mean Beat');
        plot(mean_phase, ECG_median_phase, 'linewidth', 3); lgnd = cat(2, lgnd, 'Median Beat');
        legend(lgnd);
        grid
        xlabel('phase (rad.)');
        ylabel('Amplitude (mV)');
        set(gca, 'fontsize', 14)
        set(gca, 'box', 'on')
        title(h2, 'Phase-domain averaging');
        aa = axis;
        aa(1) = mean_phase(1);
        aa(2) = mean_phase(end);
        axis(aa);
        hold off
    end

    %             z_offset = 1;
    %             figure
    %             hold on
    %             for kk = 1 : size(stacked_beats, 1)
    %                 plot3(tm, (kk-1)*z_offset*ones(1, event_width), stacked_beats(kk, :), 'color', 0.7*ones(1, 3));
    %                 %                 scatter(time, data_preprocessed_env, 28, qrs_likelihood(:)*[1 0 0], 'filled');
    %
    %             end
    %             grid
    %             xlabel('time (ms)')
    %             ylabel('beat index')
    %             zlabel('Amplitude (mV)')
    waterfall(ones(NumBeats, 1)*tm, cumsum(ones(NumBeats, event_width), 1), stacked_beats);
    %             waterfall((ones(NumBeats, 1)*tm)', cumsum(ones(NumBeats, event_width), 1)', stacked_beats')
    xlabel('time (ms)')
    ylabel('beat index')
    zlabel('Amplitude (mV)')
    set(gca, 'fontsize', 16)
    %     title(h2, 'Beat waterfall');
    view(25, 30);

    if params.save_output_images
        saveas(gcf, [params.image_name, '.png']);
    end
end
%                 sgtitle(['ID: ' subjectcode, ', Ch #', num2str(ch), ', SNR=', num2str(snr), 'dB'], 'fontsize', 18);
%             CD = get (h0, 'CData');
%             CD(1,:) = nan;
%             CD(end-2:end,:) = nan;
%             set (h0, 'CData', CD)

% make the results in JSON format
if params.save_output_files
    s = [];
    s.x_raw = x_raw;
    s.baseline_reference = baseline_reference;
    s.x_mains_cancelled = x_mains_cancelled;
    s.x_common_mode_den = x_common_mode_den;
    s.x_denoised = x_denoised;
    s.peaks_all_channels = peaks_all_channels;
    s.peaks = peaks;
    s.peak_indexes = peak_indexes;
    s.snr = snr;
    s.max_snr = max_snr;
    s.max_snr_index = max_snr_index;
    s.lead_I_polarity = alpha;

    s.ECG_mean = ECG_mean;
    s.ECG_median = ECG_median;
    s.ECG_robust_mean = ECG_robust_mean;
    s.ECG_var_mn = ECG_var_mn;
    s.ECG_robust_median = ECG_robust_median;
    s.ECG_var_md = ECG_var_md;

    JSONformatted = jsonencode(s);
    fname_results = [params.file_out_name, '.json'];
    fid = fopen(fname_results, 'w');
    fprintf(fid, '%s', JSONformatted);
    fclose(fid);
end