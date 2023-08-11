% Test script for removing ventricular activity

function testVentricularWaveSuppression
close all
clear
clc

% correct data sign or not, to have an upwards peak (positive skew)?
correct_sign = true;

% Baseline wander removal filter
baseline_removal_method = 'MDMDMDMN'; % 'LP' or 'MDMN' or 'MDMDMDMN'; % baseline wander removal method

% mains removal
mains_notch_method = 'BYPASS'; % 'IIRNOTCH', 'KFNOTCH, 'KFNOTCH2'
f_mains = 60.0; % mains frequency
Q_factor = 50; % Q-factor of the mains notch filter

% Denoising method
denoising_method = 'BYPASS';%'MDMN', 'LP', 'BYPASS', 'SEGMENTWISETIKHONOV', 'TIKHONOV', 'WAVELET', 'LP';

% R-peak detector
peak_det_method = 'PeakDetectionMultiStage'; % 'PeakDetection', 'PeakDetectionMultiStage', 'PeakDetectionProbabilistic'

% Plot results or not?
plot_results = true;

% Save results or not?
save_results = false;

% Save the spectrogram inages
save_the_output_images = false;


% Load data
% base_input_data_file_path = '../../../DataFiles/Physionet.org/files/ptbdb/1.0.0/'; % Change this path to where you have the .mat data files
% base_output_data_file_path = '../../../DataFiles/ECGFeatureResults/Physionet.org/files/ptbdb/1.0.0/'; % Change this path to where you have the .mat data files
% fs = 1000.0; % Sampling frequency of the data (put it in the loop and read it from the data if not fixed across all records)
% sub_folders = {''};

% base_input_data_file_path = '../../../DataFiles/Physionet.org/files/qtdb/1.0.0/'; % Change this path to where you have the .mat data files
% base_output_data_file_path = '../../../DataFiles/ECGFeatureResults/Physionet.org/files/qtdb/1.0.0/'; % Change this path to where you have the .mat data files
% fs = 250.0; % Sampling frequency of the data (put it in the loop and read it from the data if not fixed across all records)
% sub_folders = {''};

base_input_data_file_path = '../../../DataFiles/Physionet.org/files/challenge-2017/1.0.0/training/'; % Change this path to where you have the .mat data files
base_output_data_file_path = '../../../DataFiles/ECGFeatureResults/Physionet.org/files/challenge-2017/1.0.0/training/'; % Change this path to where you have the .mat data files
base_output_image_file_path = '../../../DataFiles/ECGImageFeatureResults/Physionet.org/files/challenge-2017/1.0.0/training/'; % Change this path to where you have the .mat data files
fs = 300.0; % Sampling frequency of the data (put it in the loop and read it from the data if not fixed across all records)
sub_folders = {'A00', 'A01', 'A02', 'A03', 'A04', 'A05', 'A06', 'A07', 'A08'};

for mm = 1 : length(sub_folders)
    input_data_file_path = [fullfile(base_input_data_file_path, sub_folders{mm}), '/']; % Change this path to where you have the .mat data files
    output_data_file_path = [fullfile(base_output_data_file_path, sub_folders{mm}), '/']; % Change this path to where you have the .mat data files
    if exist('base_output_image_file_path', 'var')
        output_image_file_path = [fullfile(base_output_image_file_path, sub_folders{mm}), '/']; % Change this path to where you have the .mat data files
    end

    input_file_list = dir(fullfile(input_data_file_path, '*.mat'));  % get list of all mat files of interest
    if ~exist(output_data_file_path, 'dir') && save_results % create output directory if it does not exist
        mkdir(output_data_file_path)
    end

    if exist('output_image_file_path', 'var')
        if ~exist(output_image_file_path, 'dir') && save_the_output_images % create output directory if it does not exist
            mkdir(output_image_file_path)
        end
    end

    for k = 1 : length(input_file_list) % Sweep over all or specific records
        input_data_file_name = fullfile(input_file_list(k).folder, input_file_list(k).name);
        data = load(input_data_file_name);
        data = data.val;

        % Sign correction (if selected)
        if correct_sign
            sign_correct_params.fc1 = 3.0;
            sign_correct_params.fc2 = 40.0;
            sign_correct_params.fs = fs;
            data = CorrectSign_(data, sign_correct_params);
        end

        % Baseline wander estimation
        baseline_params.baseline_removal_method = baseline_removal_method;
        baseline_params.fs = fs;
        baseline = ExtractBaseline_(data, baseline_params);

        data_baseline_removed = data - baseline;

        % Remove powerline noise
        powerline_params.mains_notch_method = mains_notch_method;
        powerline_params.f_mains = f_mains;
        powerline_params.fs = fs;
        powerline_params.Q_factor = Q_factor;
        data_mains_removed = RemovePowerline_(data_baseline_removed, powerline_params);

        % Remove other in-band noise
        denoiser_params.fs = fs;
        denoiser_params.denoising_method = denoising_method;
        data_denoised = RemoveNoise_(data_mains_removed, denoiser_params);

        % Detect the R-peaks
        r_peak_detector_params.ref_ch = 1;
        r_peak_detector_params.peak_det_method = peak_det_method;
        r_peak_detector_params.fs = fs;
        [peaks, peaks_indexes, rr_intervals, hr, hr_max, hr_min, hr_mean, hr_median, hr_std] = DetectRPeaks_(data_denoised, r_peak_detector_params);

        % Ventricular activity cancellation by Gaussian Process filtering (requries R-peaks)
        GPfilterparams.bins = 100;              %300; % number of phase domain bins
        GPfilterparams.BEAT_AVG_METHOD = 'MEAN'; % 'MEAN' or 'MEDIAN'
        GPfilterparams.NOISE_VAR_EST_METHOD = 'AVGLOWER'; %'MIN', 'AVGLOWER', 'MEDLOWER', 'PERCENTILE'
        GPfilterparams.p = 0.5;
        GPfilterparams.avg_bins = 10;
        GPfilterparams.SMOOTH_PHASE = 'GAUSSIAN';
        GPfilterparams.gaussianstd = 1.0;
        %  GPfilterparams.SMOOTH_PHASE = 'MA';
        %  GPfilterparams.wlen_phase = 3;
        %  GPfilterparams.wlen_time = 3;
        GPfilterparams.plotresults = 0;
        GPfilterparams.nvar_factor = 1.0; % noise variance over/under estimation factor (1 by default)
        % GPfilterparams.nvar = var(noise);

        %         GPfilterparams.avg_beat_shaping_window = round(GPfilterparams.bins/5);
        %         GPfilterparams.avg_beat_shaping_window = round(GPfilterparams.bins/20);
        GPfilterparams.Q_onset_lead = round(GPfilterparams.bins/20);
        GPfilterparams.T_offset_lag = round(GPfilterparams.bins/2.0);
        [data_posterior_est_phase_based, data_prior_est_phase_based] = ECGPhaseDomainMAPFilter(data_mains_removed, peaks, GPfilterparams);
        %         GPfilterparams.avg_beat_shaping_window = round(0.1*fs);
        %         GPfilterparams.avg_beat_shaping_window = round(0.01*fs);
        GPfilterparams.Q_onset_lead = round(0.0375*fs);
        GPfilterparams.T_offset_lag = round(0.35*fs);
        [data_posterior_est_time_based, data_prior_est_time_based] = ECGTimeDomainMAPFilter(data_mains_removed, peaks, GPfilterparams);

        if save_results
            output_data_file_name = fullfile(output_data_file_path, ['Features_', input_file_list(k).name]);
            save(output_data_file_name, 'ECG_mean_ph', 'ECG_median_ph', 'ECG_cov_ph', 'ECG_mean_robust_ph', 'ECG_med_robust_ph', 'beats_sqi_ph',...
                'ECG_mean_tm', 'ECG_median_tm', 'ECG_cov_tm', 'ECG_mean_robust_tm', 'ECG_med_robust_tm', 'beats_sqi_tm',...
                'eigen_beats_ph', 'eigen_beats_tm',...
                'peaks_indexes', 'hr', 'hr_mean', 'hr_median', 'hr_max', 'hr_min', 'hr_std',...
                '-v4');
        end

        if plot_results
            % {
            close all
            for ch = 1 : size(data, 1)
                lgnd = {};
                t = (0:size(data, 2)-1)/fs;
                figure
                plot(t, data_mains_removed(ch, :)); lgnd = cat(2, lgnd, {'Mains Removed'});
                hold on
                %                 plot(t, data_mains_removed(ch, :) - data_prior_est_time_based(ch, :)); lgnd = cat(2, lgnd, {'Prior Time Domain Residual'});
                %                 plot(t, data_mains_removed(ch, :) - data_prior_est_phase_based(ch, :)); lgnd = cat(2, lgnd, {'Prior Phase Domain Residual'});
                plot(t, data_mains_removed(ch, :) - data_posterior_est_time_based(ch, :)); lgnd = cat(2, lgnd, {'Posterior Time Domain Residual'});
                plot(t, data_mains_removed(ch, :) - data_posterior_est_phase_based(ch, :)); lgnd = cat(2, lgnd, {'Posterior Phase Domain Residual'});
                grid
                legend(lgnd)

                figure
                nfft = 512;
                subplot(311)
                spectrogram(data_mains_removed(ch, :), hamming(round(0.25*fs)), round(0.20*fs), nfft, fs, 'yaxis', 'psd');
                subplot(312)
                spectrogram(data_mains_removed(ch, :) - data_posterior_est_time_based(ch, :), hamming(round(0.25*fs)), round(0.20*fs), nfft, fs, 'yaxis', 'psd');
                subplot(313)
                spectrogram(data_mains_removed(ch, :) - data_posterior_est_phase_based(ch, :), hamming(round(0.25*fs)), round(0.20*fs), nfft, fs, 'yaxis', 'psd');

                figure
                lgnd = {};
                hold on
                periodogram(data_mains_removed(ch, :),hamming(length(data_mains_removed(ch, :))),nfft,fs); lgnd = cat(2, lgnd, {'Mains Removed'});
                periodogram(data_mains_removed(ch, :) - data_posterior_est_time_based(ch, :),hamming(length(data_mains_removed(ch, :))),nfft,fs); lgnd = cat(2, lgnd, {'Posterior Time Domain Residual'});
                periodogram(data_mains_removed(ch, :) - data_posterior_est_phase_based(ch, :),hamming(length(data_mains_removed(ch, :))),nfft,fs); lgnd = cat(2, lgnd, {'Posterior Phase Domain Residual'});
                legend(lgnd)
                %                 graphs = get(gca, 'Children');
                %                 set(graphs(2), 'Color', 'r');
            end
            % }

            t = (0:size(data, 2)-1)/fs;
            for ch = 1 : size(data, 1)
                lgnd = {};
                figure
                hold on
                plot(t, data(ch, :)); lgnd = cat(2, lgnd, {'Raw'});
                plot(t, data_baseline_removed(ch, :)); lgnd = cat(2, lgnd, {'Baseline removed'});
                plot(t, data_mains_removed(ch, :)); lgnd = cat(2, lgnd, {'Mains removed'});
                plot(t, data_denoised(ch, :)); lgnd = cat(2, lgnd, {'Denoised'});
                plot(t, data_posterior_est_time_based(ch, :)); lgnd = cat(2, lgnd, {'Denoised with Time-based GP'});
                plot(t, data_posterior_est_phase_based(ch, :)); lgnd = cat(2, lgnd, {'Denoised with Phase-based GP'});
                plot(t, baseline(ch, :), 'linewidth', 2); lgnd = cat(2, lgnd, {'Baseline'});
                legend(lgnd);
                set(gca, 'fontsize', 18);
                grid
            end


        end

        if(save_the_output_images)
            t = (0:size(data, 2)-1)/fs;
            for ch = 1 : size(data, 1)
                peak_width = round(0.01*fs);
                nfft = 1024;
                %             [S_baseline, F, T] = spectrogram(baseline(ch, :), hamming(round(0.3*fs)), round(0.25*fs), nfft, fs, 'yaxis', 'psd');
                % [S_raw, F, T] = spectrogram(data(ch, :), hamming(round(0.25*fs)), round(0.20*fs), nfft, fs, 'yaxis', 'psd');
                [S_denoised, F, T] = spectrogram(data_denoised(ch, :), hamming(round(0.25*fs)), round(0.20*fs), nfft, fs, 'yaxis', 'psd');
                %                 S = abs([S_denoised ; S_raw(1:nfft/4, :)]);
                S = abs(S_denoised);
                standard_length = length(data(ch, :))/(30.*fs);
                fh = figure('Position', [10, 10, round(1900*standard_length), 900],'PaperUnits', 'points', 'PaperOrientation', 'landscape', 'visible','off');
                ax = axes('position',[0 0 1 1]);
                s = mesh(gca, T, (1:size(S, 1)), 10*log10(S));
                s.FaceColor = 'flat';
                view(2)

                baseline_scaled = data(ch, :) - min(data(ch, :));
                baseline_scaled = (length(F)/2)*baseline_scaled / max(baseline_scaled(:)) + 1.0*length(F);

                data_scaled = data_denoised(ch, :) - min(data_denoised(ch, :));
                data_scaled = (length(F)/2)*data_scaled / max(data_scaled(:)) + length(F)/2;

                peaks_scaled = (length(F)/2)*peaks - 5;
                peaks_scaled(peaks == 0) = nan;

                hold on
                plot3(t, data_scaled, max(S(:))*ones(1, length(t)), 'r', 'linewidth', 1)
                %                 plot3(t, baseline_scaled, max(S(:))*ones(1, length(t)), 'r', 'linewidth', 1)
                plot3(t, peaks_scaled, max(S(:))*ones(1, length(t)), 's', 'MarkerFaceColor', 'm', 'MarkerEdgeColor', 'm')
                axis tight
                set(gca,'XTick',[], 'YTick', [])
                set(gca,'XColor','none','YColor','none')
                set(gca, 'color', 'none');
                ax.Toolbar.Visible = 'off';

                output_data_file_name = fullfile(output_image_file_path, ['Images_', input_file_list(k).name]);
                saveas(gcf, [output_data_file_name(1:end-4) '.png']);
                close(fh);
            end
        end

        disp(input_file_list(k).name);
    end
end
end

% A local function for sign correction
function data = CorrectSign_(data, params)
for ch = 1 : size(data, 1)
    data_hp = data(ch, :) - LPFilter(data(ch, :), params.fc1/params.fs);
    data_lp = LPFilter(data_hp, params.fc2/params.fs);
    polarity = sign(skew(data_lp));
    data(ch, :) = polarity * data(ch, :);
end
end

% A local function for baseline extraction
function baseline = ExtractBaseline_(data, params)
switch params.baseline_removal_method
    case 'BYPASS'
        baseline = zeros(size(data));
    case 'LP'
        baseline_params.fc = 0.5;
        baseline_params.fs = params.fs;
        baseline = BaselineEstimator(data, 'LP', baseline_params);
    case 'MDMN'
        w1 = 0.42; % First stage baseline wander removal window size (in seconds)
        w2 = 0.57; % Second stage baseline wander removal window size (in seconds)
        baseline_params.wlen1 = round(w1 * params.fs);
        baseline_params.wlen2 = round(w2 * params.fs);
        baseline = BaselineEstimator(data, 'MDMN', baseline_params);
    case 'MDMDMDMN'
        w1 = 0.57; % First stage baseline wander removal window size (in seconds)
        w2 = 0.60; % Second stage baseline wander removal window size (in seconds)
        w3 = 0.63; % First stage baseline wander removal window size (in seconds)
        w4 = 0.60; % Second stage baseline wander removal window size (in seconds)
        baseline_params.wlen1 = round(w1 * params.fs);
        baseline_params.wlen2 = round(w2 * params.fs);
        baseline_params.wlen3 = round(w3 * params.fs);
        baseline_params.wlen4 = round(w4 * params.fs);
        baseline = BaselineEstimator(data, 'MDMDMDMN', baseline_params);
    otherwise
        warning('Unknown method. Bypassing baseline wander removal.');
end
end

% A local function for powerline noise removal
function data_mains_removed = RemovePowerline_(data, params)
switch params.mains_notch_method % mains removal method
    case 'BYPASS' % bypass the notch filter
        data_mains_removed = data;
    case 'IIRNOTCH' % second-order IIR notch filter
        W1 = params.f_mains/(params.fs/2);
        BW1 = W1/params.Q_factor;
        [b1,a1] = iirnotch(W1,BW1);
        data_mains_removed = filtfilt(b1, a1, data')';
        % remove second harmonic
        if 0
            W2 = (2.0*params.f_mains)/(params.fs/2);
            BW2 = W2/Q_factor;
            [b2,a2] = iirnotch(W2,BW2);
            data_mains_removed = filtfilt(b2, a2, data_mains_removed')';
        end
    case 'KFNOTCH' % standard KF implementation
        W1 = params.f_mains/(params.fs/2);
        BW1 = W1/Q_factor;
        [b1,a1] = iirnotch(W1,BW1);
        data_mains_removed1 = filtfilt(b1, a1, data')';

        gamma = 0.999;
        data_mains_removed = zeros(size(data));
        for jj = 1 : size(data, 1)
            Q = 0.01*var(data(jj, :));
            R = var(data(jj, :) - data_mains_removed1(jj, :));
            [~, data_mains_removed(jj, :)] = KFNotch(data(jj, :), params.f_mains, params.fs, Q, R, gamma);
        end
    case 'KFNOTCH2' % simplified KF implementation (works with the Kalman filter Q/R ratio)
        W1 = params.f_mains/(params.fs/2);
        BW1 = W1/Q_factor;
        [b1,a1] = iirnotch(W1,BW1);
        data_mains_removed1 = filtfilt(b1, a1, data')';

        wlen = round(0.1*params.fs);
        data_mains_removed = zeros(size(data));
        for jj = 1 : size(data, 1)
            Q = 0.01*var(data(jj, :));
            R = var(data(jj, :) - data_mains_removed1(jj, :));
            [~, data_mains_removed(jj, :)] = KFNotch2(data(jj, :), params.f_mains, params.fs, 1000*Q/R, wlen);
        end
    otherwise
        error('Notch filter option undefined');
end
end

% A local function for in-band noise removal
function data_denoised = RemoveNoise_(data, params)
switch params.denoising_method
    case 'BYPASS' % bypass the denoiser
        data_denoised = data;
    case 'TIKHONOV'
        % A Tikhonov regularizer baseline detector. See TikhonovRegularization for parameter choices
        DiffOrder = 2; % smoothness constraint order greater than 1
        %     DiffImpulseResponse = [1 2 1]; % smoothness filter impulse (can replace DiffOrder in some modes of TikhonovRegularization)
        lambda = (9.47)^DiffOrder; % The smoothness parameter for Tikhonov regularization (higher values penalize the roughness of the signal)
        data_denoised = TikhonovRegularization(data, DiffOrder, lambda);
    case 'SEGMENTWISETIKHONOV'
        % Segment-wise Tikhonov regularizer. See ECGSmoothnessPriorsDenoiserBW for help and references
        DiffOrder = 2; % smoothness constraint order greater than 1
        wlen = 15e-3;%100e-3; % window length (s)
        %         withwindow = 0; % with(1) or without(0) windowing
        mode = 6; % 0..7
        nvar = 0.8e-3;
        [~, data_denoised] = ECGSmoothnessPriorsDenoiserBW(data, nvar, mode, DiffOrder, round(wlen*fs), 1e-8, 250, 1, 500);
    case 'WAVELET'
        % A wavelet denoiser
        data_denoised = zeros(size(data));
        for jj = 1 : size(data, 1)
            %             data_denoised(jj, :) = wden(data_mains_removed(jj, :), 'rigrsure', 's', 'sln', 4, 'sym5'); % use SLN since we seek beat-wise performance
            data_denoised(jj, :) = wden(data(jj, :), 'sqtwolog', 'h', 'mln', 6, 'sym7'); % use SLN since we seek beat-wise performance
        end
    case 'LP'
        fc = 25.0;
        data_denoised = LPFilter(data, fc/params.fs);
    case 'MDMN'
        w1 = 0.01; % First stage baseline wander removal window size (in seconds)
        w2 = 0.01; % Second stage baseline wander removal window size (in seconds)
        denoiser_params.wlen1 = ceil(w1 * params.fs);
        denoiser_params.wlen2 = ceil(w2 * params.fs);
        data_denoised = BaselineEstimator(data, 'MDMN', denoiser_params);
    otherwise
        error('Denoiser option undefined');
end
end

function [peaks, peaks_indexes, rr_intervals, hr, hr_max, hr_min, hr_mean, hr_median, hr_std] = DetectRPeaks_(data, params)
switch params.peak_det_method
    case 'PeakDetection'
        f0 = 1.0; % approximate heart rate (in Hz) used for R-peak detection
        [peaks, peaks_indexes] = PeakDetection(data(params.ref_ch, :), f0/params.fs, 1); % peak detection
    case 'PeakDetectionMultiStage'
        fc1 = 8.0;
        fc2 = 40.0;
        ref_ch_hp = data(params.ref_ch, :) - LPFilter(data(params.ref_ch, :), fc1/params.fs);
        ref_ch_lp = LPFilter(ref_ch_hp, fc2/params.fs); % data_denoised(ref_ch, :)

        f0 = 1.2; % approximate heart rate (in Hz) used for R-peak detection
        [peaks, peaks_indexes] = PeakDetection(ref_ch_lp, f0/params.fs, 1, 3); % peak detection
    case 'PeakDetectionProbabilistic'
        peak_detector_params.saturate = 1;
        peak_detector_params.hist_search_th = 0.9;
        peak_detector_params.rpeak_search_wlen = 0.3; % MAX detectable HR (in BPM) = 60/rpeak_search_wlen
        peak_detector_params.filter_type = 'MULT_MATCHED_FILTER';%'BANDPASS_FILTER', 'MATCHED_FILTER', 'MULT_MATCHED_FILTER'
        [peaks, peaks_indexes] = PeakDetectionProbabilistic(data(params.ref_ch, :), params.fs, peak_detector_params);
    otherwise
        error('Undefined R-peak detector');
end
rr_intervals = diff(peaks_indexes)/params.fs;
hr = 60 ./ rr_intervals; % BPM
hr_mean = mean(hr);
hr_median = median(hr);
hr_std = std(hr);
hr_max = max(hr);
hr_min = min(hr);
end