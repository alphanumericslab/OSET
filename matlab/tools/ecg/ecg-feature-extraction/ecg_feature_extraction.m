function [ecg_features_vector, ecg_feature_info, ecg_features_units, ecg_fiducial_position, exit_flag] = ecg_feature_extraction(data, fs, n_eigenvalues, num_samples, norm_flag, lead_names)
% all_features_vector = extract_ecg_features_one_record(signal_name, data, fs, n_eigenvalues, n_features, f_notch)
% Description: Extract features from one record of ECG signal
%
% INPUT:
% data - ECG signal as matrix CxT which C indicates number of ECG channels (in microvolts).
% fs   - Sampling frequency of the ECG signal (in Hz).
% f_notch - Frequency for the notch filter (in Hz).
% n_eigenvalues - Number of eigenvalues for SVD features (scaler).
%
% DEPENDENCIES:
% 1. peak_det_likelihood function from the OSET package
% 2. peak_det_likelihood_long_recs function from the OSET package
% 3. fiducial_det_lsim function from the OSET package
%
% OUTPUT:
% ecg_features_vector - includes the following features per each channel:
%   1. SNR features (2)
%   2. HRV features (7)
%   3. Features related to time interval (29)
%   4. Features related to amplitude and area under curve (32)
%   5. Features related to angle (12)
%   6. Complexity analysis features (3)
%   7. SVD features (default 5)
%   7. Average morphology (default 25)
% ecg_feature_names - name of features for all channels
% ecg_features_units - unit of all extracted features
% ecg_fiducial_position - cell array of detected fiducial points of ECG
% exit_flag - with 0 indicates successful, and minus values (-K)
% indicates error in processing K channels out of C channels

% Author: Seyedeh Somayyeh Mousavi
% Date: Jan 11, 2025
% Author: Sajjad Karimi
% Date: Mar 16, 2025
% Location: Emory University, Georgia, USA
% Email: bmemousavi@gmail.com


C = size(data, 1); % number of channels

if nargin<3
    n_eigenvalues = 5;
end

if nargin < 4
    num_samples = 25;
end

if nargin < 5
    norm_flag = true;
end


if nargin<6
    for c = 1:C
        lead_names{c} = ['ch',num2str(c)];
    end
end


% Define ECG-related features
ecg_feature_names = {
    'snr_median', 'snr_mean',...
    'n_beats', 'rmssd', 'sdnn', 'hr_median', 'hr_mean', 'hr_upper_5', 'hr_lower_5', ...
    'mean_pr_interval', 'std_pr_interval', 'median_pr_interval', 'mean_qt_interval', 'std_qt_interval', 'median_qt_interval', 'qtc_b','qtc_f' ...
    'mean_st_segment', 'std_st_segment', 'median_st_segment', 'mean_pr_segment', 'std_pr_segment', 'median_pr_segment', ...
    'mean_pr_peaks_interval', 'std_pr_peaks_interval', 'median_pr_peaks_interval', 'mean_qr_peaks_interval', 'std_qr_peaks_interval', 'median_qr_peaks_interval', ...
    'mean_qrs_complex_interval', 'std_qrs_complex_interval', 'median_qrs_complex_interval', ...
    'mean_rs_peaks_interval', 'std_rs_peaks_interval', 'median_rs_peaks_interval', 'mean_rt_peaks_interval', 'std_rt_peaks_interval', 'median_rt_peaks_interval', ...
    'mean_qrs_amp', 'std_qrs_amp', 'median_qrs_amp', 'mean_qrs_area', 'std_qrs_area', 'median_qrs_area', 'mean_qrs_abs_area', 'std_qrs_abs_area', 'median_qrs_abs_area', ...
    'mean_t_amp', 'std_t_amp', 'median_t_amp', 'ratio_qrs_t_amp', 'mean_t_area', 'std_t_area', 'median_t_area',  'mean_t_abs_area', 'std_t_abs_area', 'median_t_abs_area', ...
    'mean_p_amp', 'std_p_amp', 'median_p_amp', 'ratio_qrs_p_amp', 'mean_p_area', 'std_p_area', 'median_p_area', 'mean_p_abs_area', 'std_p_abs_area', 'median_p_abs_area', ...
    'mean_st_amp', 'std_st_amp', 'median_st_amp', ...
    'mean_pr_angles', 'std_pr_angles', 'median_pr_angles', ...
    'mean_qr_angles', 'std_qr_angles', 'median_qr_angles', ...
    'mean_rs_angles', 'std_rs_angles', 'median_rs_angles', ...
    'mean_rt_angles', 'std_rt_angles', 'median_rt_angles' , ...
    'activity','mobility', 'complexity'};

% Define corresponding units for ECG features
ecg_features_units = {
    'dB', 'dB', 'scaler', 'ms', 'ms', 'bpm', 'bpm', 'bpm', 'bpm',...
    'ms', 'ms', 'ms', 'ms', 'ms', 'ms','ms', 'ms', ...
    'ms', 'ms', 'ms', 'ms', 'ms', 'ms', ...
    'ms', 'ms', 'ms', 'ms', 'ms', 'ms', ...
    'ms', 'ms', 'ms', ...
    'ms', 'ms', 'ms', 'ms', 'ms', 'ms', ...
    'mv', 'mv', 'mv', 'mv*ms', 'mv*ms', 'mv*ms', 'mv*ms', 'mv*ms', 'mv*ms', ...
    'mv', 'mv', 'mv', 'scaler', 'mv*ms', 'mv*ms', 'mv*ms', 'mv*ms', 'mv*ms', 'mv*ms',...
    'mv', 'mv', 'mv', 'scaler', 'mv*ms', 'mv*ms', 'mv*ms', 'mv*ms', 'mv*ms', 'mv*ms',...
    'mv', 'mv', 'mv',...
    'deg', 'deg', 'deg', 'deg', 'deg', 'deg', ...
    'deg', 'deg', 'deg', 'deg', 'deg', 'deg', ...
    'mv*ms','Hz', 'Hz'};

% Define features description
ecg_feature_description = {
    "Median beat SNR", "Mean beat SNR", ...
    "Number of ECG beat", "Root Mean Square of Successive Differences for HRV", "Standard Deviation of Normal-to-Normal Intervals", "Median heart-rate", "Mean heart-rate", "95% heart-rate", "5% heart-rate", ...
    "Mean PR-interval", "Standard deviation PR-interval", "Median PR-interval", "Mean QT-interval", "Standard deviation QT-interval", "Median QT-interval", "QTc (Bazett)","QTc (Fridericia)" ...
    "Mean ST-segment (QRSoff-Ton)", "Standard deviation ST-segment (QRSoff-Ton)", "Median ST-segment (QRSoff-Ton)", "Mean PR-segment (Poff-QRSon)", "Standard deviation PR-segment (Poff-QRSon)", "Median PR-segment (Poff-QRSon)", ...
    "Mean P to R peaks interval", "Standard deviation  P to R peaks interval", "Median  P to R peaks interval", "Mean  QRSon to R peak interval", "Standard deviation  QRSon to R peak interval", "Median  QRSon to R peak interval", ...
    "Mean  QRS-complex interval", "Standard deviation QRS-complex interval", "Median  QRS-complex interval", ...
    "Mean R peak to QRSoff interval", "Standard deviation R peak to QRSoff interval", "Median R peak to QRSoff interval", "Mean R peak to T peak interval", "Standard deviation R peak to T peak interval", "Median R peak to T peak interval", ...
    "Mean R-peak amplitude", "Standard deviation R-peak amplitude", "Median R-peak amplitude", 'mean_qrs_area', 'std_qrs_area', 'median_qrs_area', 'mean_qrs_abs_area', 'std_qrs_abs_area', 'median_qrs_abs_area', ...
    'mean_t_amp', 'std_t_amp', 'median_t_amp', 'ratio_qrs_t_amp', 'mean_t_area', 'std_t_area', 'median_t_area',  'mean_t_abs_area', 'std_t_abs_area', 'median_t_abs_area', ...
    'mean_p_amp', 'std_p_amp', 'median_p_amp', 'ratio_qrs_p_amp', 'mean_p_area', 'std_p_area', 'median_p_area', 'mean_p_abs_area', 'std_p_abs_area', 'median_p_abs_area', ...
    'mean_st_amp', 'std_st_amp', 'median_st_amp', ...
    'mean_pr_angles', 'std_pr_angles', 'median_pr_angles', ...
    'mean_qr_angles', 'std_qr_angles', 'median_qr_angles', ...
    'mean_rs_angles', 'std_rs_angles', 'median_rs_angles', ...
    'mean_rt_angles', 'std_rt_angles', 'median_rt_angles' , ...
    'activity','mobility', 'complexity'};


exit_flag = 0;
% n_features = 2 + 7 + 29 + 32 + 12 + 3 + n_eigenvalues + num_samples;
n_features = length(ecg_feature_names) + n_eigenvalues + num_samples; % per channel

for n = 1:n_eigenvalues
    ecg_feature_names{end+1} = ['svd', num2str(n)];
    ecg_features_units{end+1} = 'scaler';
end

for n = 1:num_samples

    if norm_flag
        ecg_feature_names{end+1} = ['norm_morphology_sample_', num2str(n)];
        ecg_features_units{end+1} = 'scaler';
    else
        ecg_feature_names{end+1} = ['morphology_sample_', num2str(n)];
        ecg_features_units{end+1} = 'mv';
    end

end



%% Feature Extraction
ecg_features_vector = [];
ecg_fiducial_position = cell(C,1);
ecg_feature_info = cell(1,C*n_features);

for c = 1:C

    data_channel = data(c, :);
    try
        % Detect peaks based on recording length
        if (size(data, 2)/fs) < 90.0
            [~, R_peaks_indexes] = peak_det_likelihood(data_channel, fs);
        else
            [~, R_peaks_indexes] = peak_det_likelihood_long_recs(data_channel, fs);
        end

        % Run ECG fiducial points detector
        % position = fiducial_det_lsim(data_channel, R_peaks_indexes, fs);
        % Run wavedet_3D
        % heasig.nsig = 1;
        % heasig.freq = fs;
        % heasig.nsamp = length(data_channel);
        % [position, ~, ~] = wavedet_3D(data_channel', R_peaks_indexes, heasig);

        flag_post_processing = 1;
        position = fiducial_det_lsim( data_channel', R_peaks_indexes, fs, flag_post_processing);

        ecg_fiducial_position{c} = position;
        R_peaks_indexes = position.R;

        % Extract all features
        [snr_features, mean_beat] = ecg_snr_features(data_channel, R_peaks_indexes);
        hrv_features = ecg_hrv_features(R_peaks_indexes, fs);
        ti_features = ecg_time_intervals_features(R_peaks_indexes, position, fs);
        amps_areas_features = ecg_area_amp_features(data_channel, position, fs);
        angles_features = ecg_angles_features(data_channel, position, fs);
        complexity_analysis_features = complexity_features(mean_beat,fs);
        svd_features = ecg_svd_features(data_channel, R_peaks_indexes, n_eigenvalues);
        sampled_values = morphology_features(mean_beat, num_samples, norm_flag);

        % Convert feature structures to arrays
        snr_features = struct2array(snr_features);
        hrv_features = struct2array(hrv_features);
        ti_features = struct2array(ti_features);
        amps_areas_features = struct2array(amps_areas_features);
        angles_features = struct2array(angles_features);
        complexity_features= struct2array(complexity_analysis_features);
        % svd_features = struct2array(svd_features);

        % Concatenate extracted features into a single vector
        all_features = [snr_features, hrv_features, ti_features, amps_areas_features, angles_features, complexity_features, svd_features, sampled_values];

        % Append to the overall feature vector
        ecg_features_vector = cat(2, ecg_features_vector, all_features);

    catch ME
        % Error handling: display message and set NaN values for problematic channel
        disp(ME.message);
        fprintf("Error in processing signal_channel: %s, %d\n", c);
        ecg_features_vector = cat(2, ecg_features_vector, nan(1, n_features));
        exit_flag = exit_flag -1;

    end

    for n = 1:length(ecg_feature_names)
        ecg_feature_info{(c-1)*n_feature+n} = [ecg_feature_names{n}, lead_names{c}];
    end


end

% Repeat cell array C times
ecg_features_units = repmat(ecg_features_units, 1, C);

end
