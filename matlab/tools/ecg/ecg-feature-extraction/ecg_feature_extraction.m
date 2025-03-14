function [ecg_features_vector, ecg_feature_info, ecg_fiducial_position, exit_flag] = ...
    ecg_feature_extraction(ecg_data, fs, lead_names, n_svd, num_morpho_samples, norm_flag, feature_list, flatten_flag)

% Description: Extract features from one record of multi-channel ECG signal
%
% INPUT:
% ecg_data - ECG signal as matrix CxT which C indicates number of ECG channels (in microvolts).
% fs   - Sampling frequency of the ECG signal (in Hz).
% lead_names - a cell array containing each ECG lead name
% n_svd - Number of eigenvalues for SVD features (default 5).
% num_morpho_samples - Number of samples for ECG average beat (default 25).
% norm_flag - A boolian value for normalizing ECG average beat before sampling (default true).
% feature_list - A cell array contaning a subset of {'snr', 'hrv','interval', 'amplitude', 'angle', 'svd', 'morpho'}
%                to return selected subset of all ECG features (default all features)

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
% vecg_feature_info - includes following information feilds
%   ecg_feature_names - name of features for all channels
%   ecg_features_units - unit of all extracted features
%   ecg_fiducial_position - cell array of detected fiducial points of ECG
% exit_flag - with 0 indicates successful, and minus values (-K)
%   indicates error in processing K channels out of C channels

% Author: Seyedeh Somayyeh Mousavi
% Date: Jan 11, 2025
% Location: Emory University, Georgia, USA
% Email: bmemousavi@gmail.com
% Author:
%   Sajjad Karimi
%   Emory University, Georgia, USA
%   Email: sajjadkarimi91@gmail.com
%   Date: Mar 14, 2025


[C, T] = size(ecg_data); % number of channels

if nargin<3  || isempty(lead_names)
    for c = 1:C
        lead_names{c} = ['ch',num2str(c)];
    end
end

if nargin<4 || isempty(n_svd)
    n_svd = 5;
end

if nargin < 5 || isempty(num_morpho_samples)
    num_morpho_samples = 25;
end

if nargin < 6 || isempty(norm_flag)
    norm_flag = true;
end


if nargin < 7 || isempty(feature_list)
    feature_list = { 'amplitude', 'angle', 'hjorth' };
end

if nargin<8
    flatten_flag = true;
end

exit_flag = 0;

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
    'degree', 'degree', 'degree', 'degree', 'degree', 'degree', ...
    'degree', 'degree', 'degree', 'degree', 'degree', 'degree', ...
    'mv*ms','Hz', 'Hz'};

% Define features description
ecg_feature_description = {
    "Median beat SNR", "Mean beat SNR", ...
    "Number of ECG beat", "Root Mean Square of Successive Differences for HRV", "Standard Deviation of Normal-to-Normal Intervals", "Median heart-rate", "Mean heart-rate", "95% heart-rate", "5% heart-rate", ...
    "Mean PR-interval", "Standard deviation PR-interval", "Median PR-interval", "Mean QT-interval", "Standard deviation QT-interval", "Median QT-interval", "QTc (Bazett)","QTc (Fridericia)" ...
    "Mean ST-segment (QRSoff-Ton)", "Standard deviation ST-segment (QRSoff-Ton)", "Median ST-segment (QRSoff-Ton)", "Mean PR-segment (Poff-QRSon)", "Standard deviation PR-segment (Poff-QRSon)", "Median PR-segment (Poff-QRSon)", ...
    "Mean P to R peaks interval", "Standard deviation P to R peaks interval", "Median P to R peaks interval", "Mean QRSon to R peak interval", "Standard deviation QRSon to R peak interval", "Median QRSon to R peak interval", ...
    "Mean  QRS-complex interval", "Standard deviation QRS-complex interval", "Median  QRS-complex interval", ...
    "Mean R peak to QRSoff interval", "Standard deviation R peak to QRSoff interval", "Median R peak to QRSoff interval", "Mean R peak to T peak interval", "Standard deviation R peak to T peak interval", "Median R peak to T peak interval", ...
    "Mean R-peak amplitude", "Standard deviation R-peak amplitude", "Median R-peak amplitude", ...
    "Mean QRS area", "Standard deviation QRS area", "Median QRS area", "Mean area of absolute QRS", "Standard deviation area of absolute QRS", "Median area of absolute QRS", ...
    "Mean T-peak amplitude", "Standard deviation T-peak amplitude", "Median T-peak amplitude", "R peak to T peak ratio", ...
    "Mean T wave area", "Standard deviation T wave area", "Median T wave area", "Mean area of absolute T wave", "Standard deviation area of absolute T wave", "Median area of absolute T wave", ...
    "Mean P-peak amplitude", "Standard deviation P-peak amplitude", "Median P-peak amplitude", "R peak to P peak ratio", ...
    "Mean P wave area", "Standard deviation P wave area", "Median P wave area", "Mean area of absolute P wave", "Standard deviation area of absolute P wave", "Median area of absolute P wave", ...
    "Mean ST-segment amplitude", "Standard deviation ST-segment amplitude", "Median ST-segment amplitude", ...
    "Mean P and R peaks angle", "Standard deviation P and R peaks angle", "Median P and R peaks angle", ...
    "Mean QRSon and R peak angle", "Standard deviation QRSon and R peak angle", "Median QRSon and R peak angle", ...
    "Mean R peak and QRSoff angle", "Standard deviation R peak and QRSoff peak angle", "Median R peak and QRSoff angle", ...
    "Mean R and T peaks angle", "Standard deviation R and T peaks angle", "Median R and T peaks angle", ...
    "Hjorth parameters: Activity","Hjorth parameters: Mobility","Hjorth parameters: Complexity"};




% n_features = length(ecg_feature_names) + n_svd + num_morpho_samples; % per channel

for n = 1:n_svd
    ecg_feature_names{end+1} = ['svd', num2str(n)];
    ecg_features_units{end+1} = 'scaler';
    ecg_feature_description{end+1} = ['Normalized SVD value ', num2str(n)];
end

for n = 1:num_morpho_samples

    if norm_flag
        ecg_feature_names{end+1} = ['norm_morphology_sample_', num2str(n)];
        ecg_features_units{end+1} = 'scaler';
        ecg_feature_description{end+1} = ['Normalized average beat sample ', num2str(n)];
    else
        ecg_feature_names{end+1} = ['morphology_sample_', num2str(n)];
        ecg_features_units{end+1} = 'mv';
        ecg_feature_description{end+1} = ['Average beat sample ', num2str(n)];
    end

end



%% Feature Extraction

feature_list_all = {'snr', 'hrv','interval', 'amplitude', 'angle', 'hjorth', 'svd', 'morpho'};
[~ , ilocb] = ismember(feature_list,feature_list_all);

num_features =         [2     7      29           32          12        3      n_svd   num_morpho_samples]; % per channel

n_features = sum(num_features(ilocb));

ecg_features_vector = [];
ecg_fiducial_position = cell(C,1);
ecg_feature_names_ch = cell(1,C*n_features);

seg_len_time = min(10,T/fs); % 10 seconds or signal length for data with shorter than 10s
pad_len_time = 1;

ecg_feature_names = [];

for c = 1:C

    data_channel = ecg_data(c, :);
    try

        % Detect peaks based on recording length
        [~, rpeak_indexes] = peak_det_likelihood_long_recs(data_channel, fs, seg_len_time, pad_len_time);

        % Run ECG fiducial points detector
        % position = fiducial_det_lsim(data_channel, R_peaks_indexes, fs);
        % Run wavedet_3D
        % heasig.nsig = 1;
        % heasig.freq = fs;
        % heasig.nsamp = length(data_channel);
        % [position, ~, ~] = wavedet_3D(data_channel', R_peaks_indexes, heasig);

        flag_post_processing = 1;
        position = fiducial_det_lsim( data_channel', rpeak_indexes, fs, flag_post_processing);

        ecg_fiducial_position{c} = position;
        rpeak_indexes = position.R;

        % Extract all features
        % [snr_features, snr_feature_info, mean_beat] = ecg_snr_features(data_channel, rpeak_indexes);
        % [hrv_features, hrv_feature_info] = ecg_hrv_features(rpeak_indexes, fs);
        % [ti_features, ti_feature_info] = ecg_time_intervals_features(rpeak_indexes, position, fs);
        % [amps_areas_features, amp_feature_info] = ecg_area_amp_features(data_channel, position, fs);
        % [angles_features, angle_feature_info] = ecg_angles_features(data_channel, position, fs);
        % [hjorth_features, hjorth_feature_info] = hjorth_time_features(mean_beat,fs);
        % [svd_features, svd_feature_info] = ecg_svd_features(data_channel, rpeak_indexes, n_svd);
        % [sample_morpho, morpho_feature_info] = morphology_features(mean_beat, num_morpho_samples, norm_flag);

        % % Concatenate extracted features into a single vector
        % all_features = [snr_features, hrv_features, ti_features, amps_areas_features, angles_features, hjorth_features, svd_features, sample_morpho];

        [featureset{1}, this_feature_info{1}, mean_beat] = ecg_snr_features(data_channel, rpeak_indexes);
        [featureset{2}, this_feature_info{2}] = ecg_hrv_features(rpeak_indexes, fs);
        [featureset{3}, this_feature_info{3}] = ecg_time_intervals_features(rpeak_indexes, position, fs);
        [featureset{4}, this_feature_info{4}] = ecg_area_amp_features(data_channel, position, fs);
        [featureset{5}, this_feature_info{5}] = ecg_angles_features(data_channel, position, fs);
        [featureset{6}, this_feature_info{6}] = hjorth_time_features(mean_beat,fs);
        [featureset{7}, this_feature_info{7}] = ecg_svd_features(data_channel, rpeak_indexes, n_svd);
        [featureset{8}, this_feature_info{8}] = morphology_features(mean_beat, num_morpho_samples, norm_flag);



        all_features = cell2mat(featureset(ilocb));

        % Append to the overall feature vector
        ecg_features_vector = cat(2, ecg_features_vector, all_features);

        this_feature_info = this_feature_info(ilocb);

    catch ME
        % Error handling: display message and set NaN values for problematic channel
        disp(ME.message);
        fprintf("Error in processing signal_channel: %s, %d\n", c);
        ecg_features_vector = cat(2, ecg_features_vector, nan(1, n_features));
        exit_flag = exit_flag -1;

    end

    if c==1 || isempty(ecg_feature_names)
        ecg_features_units = [];
        ecg_feature_description = [];

        for f = 1:length(this_feature_info)
            ecg_feature_names = [ecg_feature_names, this_feature_info{f}.names];
            ecg_features_units = [ecg_features_units, this_feature_info{f}.units];
            ecg_feature_description = [ecg_feature_description, this_feature_info{f}.description];
        end

    end

    for n = 1:length(ecg_feature_names)
        ecg_feature_names_ch{(c-1)*n_features+n} = [ecg_feature_names{n},'_', lead_names{c}];
    end


end

% Repeat cell array C times
ecg_features_units = repmat(ecg_features_units, 1, C);
ecg_feature_description = repmat(ecg_feature_description, 1, C);

ecg_feature_description = ecg_feature_description';
ecg_feature_names_ch = ecg_feature_names_ch';
ecg_features_units = ecg_features_units';

ecg_feature_info.names = ecg_feature_names_ch;
ecg_feature_info.units = ecg_features_units;
ecg_feature_info.description = ecg_feature_description;

end
