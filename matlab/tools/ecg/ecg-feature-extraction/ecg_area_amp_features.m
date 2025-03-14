function [feature_vec, feature_info] = ecg_area_amp_features(data, position, fs)
% featureset = ecg_area_amp_features(data, position, fs)
% Extract features related to the ECG amplitude and area under the
% curve.
%
% Inputs:
%   data: ECG signal (1D array).
%   position: Fiducial points of the ECG signal.
%   fs: Sampling frequency (Hz).
%
% Output:
%   featureset: Structure containing features related to the ECG
%   amplitude and area under the curve of the QRS complex, T wave and P wave,
%   amplitude of the ST segment, and amplitude ratio of the R peak to the S
%   and T waves.
%
% Author:
%   Seyedeh Somayyeh Mousavi
%   Emory University, Georgia, USA
%   Email: bmemousavi@gmail.com
%   Date: SEP 24, 2024
% Author: Sajjad Karimi
% Location: Emory University, Georgia, USA
% Email: sajjadkarimi91@gmail.com
% Date: Mar 14, 2025

%% Constant value
convert_s_ms =1000;
% Calculate the duration of each time step
time_step = 1 / fs;

%% Define ECG points
p_onset = position.Pon;
p = position.P;
p_offset = position.Poff;
qrs_onset = position.QRSon;
rpeak = position.R;
qrs_offset = position.QRSoff;
t_onset = position.Ton;
t = position.T;
t_offset = position.Toff;

%% Calculate amplitude
% Extract T-peak amplitude from the ECG signal
valid_indices_t = ~isnan(t) & ~isnan(t_onset) & ~isnan(t_offset);

% Check if any valid amplitudes are present for T-wave
if ~any(valid_indices_t)
    mean_t_amplitude = NaN;
    std_t_amplitude = NaN;
    median_t_amplitude = NaN;
    mean_t_auc = NaN;
    std_t_auc = NaN;
    median_t_auc = NaN;
    mean_t_abs_auc = NaN;
    std_t_abs_auc = NaN;
    median_t_abs_auc = NaN;

else
    t_valid = t(valid_indices_t);
    t_onset_valid = t_onset(valid_indices_t);
    t_offset_valid = t_offset(valid_indices_t);

    % Calculate statistics for T-wave
    amplitude_t_peak = data(t_valid) - (data(t_onset_valid) + data(t_offset_valid))/2;
    mean_t_amplitude = mean(amplitude_t_peak);
    std_t_amplitude = std(amplitude_t_peak);
    median_t_amplitude = median(amplitude_t_peak);

    % Compute AUC
    t_auc = nan(size(t_valid));
    t_abs_auc = nan(size(t_valid));
    for i = 1:length(t_valid)

        t_auc(i) = sum(data(t_onset_valid(i):t_offset_valid(i))) * time_step * convert_s_ms;
        t_abs_auc(i) = sum(abs(data(t_onset_valid(i):t_offset_valid(i)))) * time_step * convert_s_ms;

    end
    mean_t_auc = mean(t_auc);
    std_t_auc = std(t_auc);
    median_t_auc = median(t_auc);
    mean_t_abs_auc = mean(t_abs_auc);
    std_t_abs_auc = std(t_abs_auc);
    median_t_abs_auc = median(t_abs_auc);
end

%% Calculate amplitude
% Extract P-peak amplitude from the ECG signal
valid_indices_p = ~isnan(p) & ~isnan(p_onset) & ~isnan(p_offset);

% Check if any valid amplitudes are present for P-wave
if ~any(valid_indices_p)
    mean_p_amplitude = NaN;
    std_p_amplitude = NaN;
    median_p_amplitude = NaN;
    mean_p_auc = NaN;
    std_p_auc = NaN;
    median_p_auc = NaN;
    mean_p_abs_auc = NaN;
    std_p_abs_auc = NaN;
    median_p_abs_auc = NaN;
else
    valid_indices_p = ~isnan(p) & ~isnan(p_onset) & ~isnan(p_offset);
    p_valid = p(valid_indices_p);
    p_onset_valid = p_onset(valid_indices_p);
    p_offset_valid = p_offset(valid_indices_p);

    % Calculate statistics for P-wave
    amplitude_p_peak = data(p_valid) - (data(p_onset_valid) + data(p_offset_valid))/2;
    mean_p_amplitude = mean(amplitude_p_peak);
    std_p_amplitude = std(amplitude_p_peak);
    median_p_amplitude = median(amplitude_p_peak);

    % Compute AUC
    p_auc = nan(size(p_valid));
    p_abs_auc = nan(size(p_valid));

    for i = 1:length(p_valid)

        p_auc(i) = sum(data(p_onset_valid(i):p_offset_valid(i)))* time_step * convert_s_ms;
        p_abs_auc(i) = sum(abs(data(p_onset_valid(i):p_offset_valid(i))))* time_step * convert_s_ms;

    end

    mean_p_auc = mean(p_auc);
    std_p_auc = std(p_auc);
    median_p_auc = median(p_auc);
    mean_p_abs_auc = mean(p_abs_auc);
    std_p_abs_auc = std(p_abs_auc);
    median_p_abs_auc = median(p_abs_auc);

end
%% Calculate amplitude
% Extract ST segment amplitude from the ECG signal
valid_indices_st =  ~isnan(t_onset) & ~isnan(qrs_offset) ;

% Check if any valid amplitudes are present for the ST segment
if ~any(valid_indices_st)
    mean_st_amplitude = NaN;
    std_st_amplitude = NaN;
    median_st_amplitude = NaN;
else
    st_onset_valid = t_onset(valid_indices_st);
    qrs_offset_valid = qrs_offset(valid_indices_st);
    amplitude_st = nan(size(st_onset_valid));

    % Extract amplitudes for the valid ST segment range
    for idx = 1:length(st_onset_valid)

        amplitude_st(idx) =  mean(data(qrs_offset_valid(idx):st_onset_valid(idx)));

    end

    % Calculate statistics
    mean_st_amplitude = mean(amplitude_st);
    std_st_amplitude = std(amplitude_st);
    median_st_amplitude = median(amplitude_st);
end

%% Extract QRS-peak amplitude from the ECG signal
valid_indices_qrs = ~isnan(rpeak) & ~isnan(qrs_onset) & ~isnan(qrs_offset);

% Check if any valid amplitudes are present for QRS complex
if ~any(valid_indices_qrs)
    mean_qrs_amplitude = NaN;
    std_qrs_amplitude = NaN;
    median_qrs_amplitude = NaN;
    ratio_qrs_t_amplitude = NaN;
    ratio_qrs_p_amplitude = NaN;
    mean_qrs_auc = NaN;
    std_qrs_auc = NaN;
    median_qrs_auc = NaN;
    mean_qrs_abs_auc = NaN;
    std_qrs_abs_auc = NaN;
    median_qrs_abs_auc = NaN;
else
    qrs_valid = rpeak(valid_indices_qrs);
    qrs_onset_valid = qrs_onset(valid_indices_qrs);
    qrs_offset_valid = qrs_offset(valid_indices_qrs);

    % Calculate statistics for QRS complex
    amplitude_qrs_peak = data(qrs_valid) - (data(qrs_onset_valid) + data(qrs_offset_valid))/2;
    mean_qrs_amplitude = mean(amplitude_qrs_peak);
    std_qrs_amplitude = std(amplitude_qrs_peak);
    median_qrs_amplitude = median(amplitude_qrs_peak);
    ratio_qrs_t_amplitude = median_t_amplitude / median_qrs_amplitude;
    ratio_qrs_p_amplitude = median_p_amplitude/ median_qrs_amplitude;

    % Compute AUC
    qrs_auc = nan(size(qrs_valid));
    qrs_abs_auc = nan(size(qrs_valid));

    for i = 1:length(qrs_valid)

        qrs_auc(i) = sum(data(qrs_onset_valid(i):qrs_offset_valid(i))) * time_step * convert_s_ms;
        qrs_abs_auc(i) = sum(abs(data(qrs_onset_valid(i):qrs_offset_valid(i)))) * time_step * convert_s_ms;

    end

    mean_qrs_auc = mean(qrs_auc);
    std_qrs_auc = std(qrs_auc);
    median_qrs_auc = median(qrs_auc);
    mean_qrs_abs_auc = mean(qrs_abs_auc);
    std_qrs_abs_auc = std(qrs_abs_auc);
    median_qrs_abs_auc= median(qrs_abs_auc);
end

%%
% Store QRS complex amplitude and area results
featureset.mean_qrs_amplitude = mean_qrs_amplitude;
featureset.std_qrs_amplitude = std_qrs_amplitude;
featureset.median_qrs_amplitude = median_qrs_amplitude;
featureset.mean_qrs_auc = mean_qrs_auc;
featureset.std_qrs_auc = std_qrs_auc;
featureset.median_qrs_auc = median_qrs_auc;
featureset.mean_qrs_abs_auc = mean_qrs_abs_auc;
featureset.std_qrs_abs_auc = std_qrs_abs_auc;
featureset.median_qrs_abs_auc = median_qrs_abs_auc;

% Store T-wave amplitude results
featureset.mean_t_amplitude = mean_t_amplitude;
featureset.std_t_amplitude = std_t_amplitude;
featureset.median_t_amplitude = median_t_amplitude;
featureset.ratio_qrs_t_amplitude = ratio_qrs_t_amplitude;
featureset.mean_t_auc = mean_t_auc;
featureset.std_t_auc = std_t_auc;
featureset.median_t_auc = median_t_auc;
featureset.mean_t_abs_auc = mean_t_abs_auc;
featureset.std_t_abs_auc = std_t_abs_auc;
featureset.median_t_abs_auc = median_t_abs_auc;

% Store P-wave amplitude results
featureset.mean_p_amplitude = mean_p_amplitude;
featureset.std_p_amplitude = std_p_amplitude;
featureset.median_p_amplitude = median_p_amplitude;
featureset.ratio_qrs_p_amplitude = ratio_qrs_p_amplitude;
featureset.mean_p_auc = mean_p_auc;
featureset.std_p_auc = std_p_auc;
featureset.median_p_auc = median_p_auc;
featureset.mean_p_abs_auc = mean_p_abs_auc;
featureset.std_p_abs_auc = std_p_abs_auc;
featureset.median_p_abs_auc = median_p_abs_auc;

% Store ST-segment amplitude results
featureset.mean_st_amplitude = mean_st_amplitude;
featureset.std_st_amplitude = std_st_amplitude;
featureset.median_st_amplitude = median_st_amplitude;

feature_vec = [mean_qrs_amplitude, std_qrs_amplitude, median_qrs_amplitude, mean_qrs_auc, std_qrs_auc , median_qrs_auc, mean_qrs_abs_auc, std_qrs_abs_auc, median_qrs_abs_auc,...
    mean_t_amplitude, std_t_amplitude, median_t_amplitude, ratio_qrs_t_amplitude, mean_t_auc, std_t_auc, median_t_auc, mean_t_abs_auc, std_t_abs_auc, median_t_abs_auc, ...
    mean_p_amplitude, std_p_amplitude, median_p_amplitude, ratio_qrs_p_amplitude, mean_p_auc, std_p_auc, median_p_auc, mean_p_abs_auc, std_p_abs_auc, median_p_abs_auc, ...
    mean_st_amplitude, std_st_amplitude, median_st_amplitude];
% Define feature info
feature_info.names = {    'mean_qrs_amp', 'std_qrs_amp', 'median_qrs_amp', 'mean_qrs_area', 'std_qrs_area', 'median_qrs_area', 'mean_qrs_abs_area', 'std_qrs_abs_area', 'median_qrs_abs_area', ...
    'mean_t_amp', 'std_t_amp', 'median_t_amp', 'ratio_t_r_amp', 'mean_t_area', 'std_t_area', 'median_t_area',  'mean_t_abs_area', 'std_t_abs_area', 'median_t_abs_area', ...
    'mean_p_amp', 'std_p_amp', 'median_p_amp', 'ratio_p_r_amp', 'mean_p_area', 'std_p_area', 'median_p_area', 'mean_p_abs_area', 'std_p_abs_area', 'median_p_abs_area', ...
    'mean_st_amp', 'std_st_amp', 'median_st_amp'};

feature_info.units = {'mv', 'mv', 'mv', 'mv*ms', 'mv*ms', 'mv*ms', 'mv*ms', 'mv*ms', 'mv*ms', ...
    'mv', 'mv', 'mv', 'scaler', 'mv*ms', 'mv*ms', 'mv*ms', 'mv*ms', 'mv*ms', 'mv*ms',...
    'mv', 'mv', 'mv', 'scaler', 'mv*ms', 'mv*ms', 'mv*ms', 'mv*ms', 'mv*ms', 'mv*ms',...
    'mv', 'mv', 'mv'};

feature_info.description = {"Mean R-peak amplitude", "Standard deviation R-peak amplitude", "Median R-peak amplitude", ...
    "Mean QRS area", "Standard deviation QRS area", "Median QRS area", "Mean area of absolute QRS", "Standard deviation area of absolute QRS", "Median area of absolute QRS", ...
    "Mean T-peak amplitude", "Standard deviation T-peak amplitude", "Median T-peak amplitude", "T peak to R peak ratio", ...
    "Mean T wave area", "Standard deviation T wave area", "Median T wave area", "Mean area of absolute T wave", "Standard deviation area of absolute T wave", "Median area of absolute T wave", ...
    "Mean P-peak amplitude", "Standard deviation P-peak amplitude", "Median P-peak amplitude", "P peak to R peak ratio", ...
    "Mean P wave area", "Standard deviation P wave area", "Median P wave area", "Mean area of absolute P wave", "Standard deviation area of absolute P wave", "Median area of absolute P wave", ...
    "Mean ST-segment amplitude", "Standard deviation ST-segment amplitude", "Median ST-segment amplitude"};


end
