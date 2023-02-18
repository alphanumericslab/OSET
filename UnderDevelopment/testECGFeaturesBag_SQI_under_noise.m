% Test script for an ECG denoiser based on a data driven MAP estimator in
% the phase and time domains

close all
clear
clc

% Load data
% datafilepath = '../../../../DataFiles/Physionet.org/files/ptbdb/1.0.0/'; % Change this path to where you have the .mat data files
% fs = 1000.0; % Sampling frequency of the data (put it in the loop and read it from the data if not fixed across all records)
% ofname = 'PTBAnalysisResults.csv';

datafilepath = '../../../DataFiles/Physionet.org/files/qtdb/1.0.0/'; % Change this path to where you have the .mat data files
fs = 250.0; % Sampling frequency of the data (put it in the loop and read it from the data if not fixed across all records)
ofname = 'QTAnalysisResults.csv';

filelist = dir(fullfile([datafilepath, '**/*.mat']));  % get list of all mat files of interest
% Make N by 2 matrix of fieldname + value type
variable_names_types = [["FileID", "string"]; ...
    ["Channel", "int32"]; ...
    ["SignalLen", "int32"]; ...
    ["NumRPeaks", "int32"]; ...
    ["Run", "int32"]; ...
    ["SetInSNR", "double"]; ...
    ["ActInSNR", "double"]; ...
    ["SignalVar", "double"]; ...
    ["NoiseVar", "double"]; ...
    ["OutSNR_PriorPhaseBased", "double"]; ...
    ["OutSNR_PostPhaseBased", "double"]; ...
    ["OutSNR_PriorTimeBased", "double"]; ...
    ["OutSNR_PostTimeBased", "double"]];

% Make table using fieldnames & value types from above
ResultsTableEmpty = table('Size',[1,size(variable_names_types,1)],...
    'VariableNames', variable_names_types(:,1),...
    'VariableTypes', variable_names_types(:,2));

ResultsTable = ResultsTableEmpty;
% Baseline wander removal filter
w1 = 0.72; % First stage baseline wander removal window size (in seconds)
w2 = 0.87; % Second stage baseline wander removal window size (in seconds)
BASELINE_REMOVAL_APPROACH = 'MDMN'; % 'BP' or 'MDMN'; % baseline wander removal method
OneRowSet = false;
for k = 15 : length(filelist) % Sweep over all or specific records
    datafilename = [filelist(k).folder '/' filelist(k).name];
    data = load(datafilename);
    data = data.val;

    switch(BASELINE_REMOVAL_APPROACH)
        case 'BP'
            data = data - LPFilter(data, 5.0/fs);
            data = LPFilter(data, 80.0/fs);
        case 'MDMN'
            wlen1 = round(w1 * fs);
            wlen2 = round(w2 * fs);
            for jj = 1 : size(data, 1)
                bl1 = BaseLine1(data(jj, :), wlen1, 'md');
                data(jj, :) = data(jj, :) - BaseLine1(bl1, wlen2, 'mn');
            end
        otherwise
            warning('Unknown method. Bypassing baseline wander removal.');
    end

    % for SNR_pre_set = 55 : 5 : 60 % the desired input SNR
    for ch = 1 : size(data, 1) % sweep over all or a single desired channel
        sig = data(ch, :);
%         sd = sqrt(var(sig) / 10^(SNR_pre_set/10));
%         noise = sd * randn(size(sig));
%         % x = sig + noise;


        noise0 = zeros(size(sig));
        noisy_indexes0 = round(0.05*length(noise0)) : round(0.2*length(noise0));
        noise0(noisy_indexes0) = randn(1, length(noisy_indexes0));

        noise1 = zeros(size(sig));
        noisy_indexes1 = round(0.45*length(noise1)) : round(0.65*length(noise1));
        noise1(noisy_indexes1) = randn(1, length(noisy_indexes1));

        noise2 = zeros(size(sig));
        noisy_indexes2 = round(0.75*length(noise2)) : round(0.95*length(noise2));
        noise2(noisy_indexes2) = randn(1, length(noisy_indexes2));

        x = sig + (0.5*std(sig))*noise0 + (0.7*std(sig))*noise1 + (0.9*std(sig))*noise2;

        if 0
            f0 = 1.0; % approximate heart rate (in Hz) used for R-peak detection
            [peaks, peaks_indexes] = PeakDetection(x, f0/fs);                  % peak detection
        else
            peak_detector_params.saturate = 1;
            peak_detector_params.hist_search_th = 0.9;
            peak_detector_params.rpeak_search_wlen = 0.3; % MAX detectable HR (in BPM) = 60/rpeak_search_wlen
            peak_detector_params.filter_type = 'MULT_MATCHED_FILTER';%'BANDPASS_FILTER', 'MATCHED_FILTER', 'MULT_MATCHED_FILTER'
            [peaks, peaks_indexes] = PeakDetectionProbabilistic(x, fs, peak_detector_params);
        end

        FeatureParams.window = true;
        FeatureParams.BEAT_AVG_METHOD = 'ROBUST_MEAN'; % 'MEAN' or 'MEDIAN' or 'ROBUST_MEAN' or 'ROBUST_MEDIAN'
        FeatureParams.bins = 100;

        [ECG_mean_ph, ECG_median_ph, ECG_cov_ph, ECG_mean_robust_ph, ECG_median_robust_ph, stacked_beats_ph, beats_sqi_ph] = PhaseDomainECGFeatures(x, peaks, FeatureParams);
        FeatureParams.width = 401;
        [ECG_mean_tm, ECG_median_tm, ECG_cov_tm, ECG_mean_robust_tm, ECG_median_robust_tm, stacked_beats_tm, beats_sqi_tm] = TimeDomainECGFeatures(x, peaks, FeatureParams);

        t = (0:length(x)-1)/fs;
        figure
        subplot(211);
        rectangle('Position',[t(noisy_indexes0(1)), min(x(:)), t(noisy_indexes0(end))- t(noisy_indexes0(1)), max(x(:)) - min(x(:))], 'EdgeColor', 'm')
        hold on
        rectangle('Position',[t(noisy_indexes1(1)), min(x(:)), t(noisy_indexes1(end))- t(noisy_indexes1(1)), max(x(:)) - min(x(:))], 'EdgeColor', 'm')
        rectangle('Position',[t(noisy_indexes2(1)), min(x(:)), t(noisy_indexes2(end))- t(noisy_indexes2(1)), max(x(:)) - min(x(:))], 'EdgeColor', 'm')
        plot(t, x);
        plot(t(peaks_indexes), x(peaks_indexes), 'ro');
        grid
        xlabel('time(s)');
        ylabel('Amplitude(uV)');
        set(gca, 'fontsize', 16);
        legend('ECG', 'R-peaks');

        subplot(212)
        rectangle('Position',[t(noisy_indexes0(1)), 0, t(noisy_indexes0(end))- t(noisy_indexes0(1)), 1], 'EdgeColor', 'm')
        hold on
        rectangle('Position',[t(noisy_indexes1(1)), 0, t(noisy_indexes1(end))- t(noisy_indexes1(1)), 1], 'EdgeColor', 'm')
        rectangle('Position',[t(noisy_indexes2(1)), 0, t(noisy_indexes2(end))- t(noisy_indexes2(1)), 1], 'EdgeColor', 'm')
        % plot(t(peaks_indexes), LPFilter(beats_sqi_tm', 150.0/fs));
        % plot(t(peaks_indexes(1:end-1)), LPFilter(beats_sqi_ph', 150.0/fs));
        plot(t(peaks_indexes), beats_sqi_tm);
        plot(t(peaks_indexes(1:end-1)), beats_sqi_ph);
        grid
        xlabel('time(s)');
        ylabel('Beatwise SQI');
        set(gca, 'fontsize', 16);
        legend('Time-based algorithm', 'Phase-based algorithm');

    end
end
