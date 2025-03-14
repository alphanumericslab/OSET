function [feature_vec, feature_info] = ecg_time_intervals_features(rpeak_indexes, position, fs)
% featureset = ecg_time_intervals_features(R_peaks_indexes, position, fs)
% Extract features related to the ECG time interval
%
% Inputs:
%   rpeak_indexes - A vector containing the R-peak indices of the ECG signal (expressed as sample points).
%   position: Fiducial points of the ECG signal (expressed as sample points).
%   fs: Sampling frequency (Hz).
%
% Output:
%   featureset: Structure containing features related to the ECG time
%   intervals, including:
%       - QRS complex
%       - QT interval
%       - PR interval
%       - ST interval
%       - PR segment
%       - ST segment
%       - Time intervals between peaks, such as:
%           * P to R peaks
%           * Q to R peaks
%           * S to R peaks
%           * T to R peaks
%       - Corrected QT intervals (QTC):
%           * QTC_F (Fridericia correction)
%           * QTC_B (Bazett correction)
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
convert_ms_s = 1/ convert_s_ms ;

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

%% Calculate intervals
rr_interval = convert_s_ms * diff(rpeak_indexes) / fs;
mean_rr_interval = mean(rr_interval);

qrs_complex = convert_s_ms * (qrs_offset - qrs_onset) / fs;
qt_interval = convert_s_ms * (t_offset - qrs_onset) / fs;
pr_interval = convert_s_ms * (qrs_onset - p_onset) / fs ;
st_segment = convert_s_ms * (t_onset - qrs_offset) / fs ;
pr_segment = convert_s_ms * (qrs_onset - p_offset) / fs ;
rt_peaks = convert_s_ms * ( t - rpeak) / fs;
pr_peaks = convert_s_ms * (rpeak - p) / fs;
qr_peaks = convert_s_ms * (rpeak - qrs_onset) / fs;
rs_peaks = convert_s_ms * (qrs_offset - rpeak) / fs;

% Omit NaN values
valid_qrs_complex = qrs_complex(~isnan(qrs_complex));
valid_qt_interval = qt_interval(~isnan(qt_interval));
valid_pr_interval = pr_interval(~isnan(pr_interval));
valid_st_segment = st_segment(~isnan(st_segment));
valid_pr_segment = pr_segment(~isnan(pr_segment));
valid_rt_peaks = rt_peaks(~isnan(rt_peaks));
valid_pr_peaks = pr_peaks(~isnan(pr_peaks));
valid_qr_peaks = qr_peaks(~isnan(qr_peaks));
valid_rs_peaks = rs_peaks(~isnan(rs_peaks));

%% Check if any valid intervals are present for qrs_complex
if isempty(valid_qrs_complex)
    mean_qrs = NaN;
    std_qrs = NaN;
    median_qrs = NaN;
else
    % Calculate statistics
    mean_qrs = mean(valid_qrs_complex);
    std_qrs = std(valid_qrs_complex);
    median_qrs = median(valid_qrs_complex);
end

%% Check if any valid intervals are present for qt_interval
if isempty(valid_qt_interval)
    mean_qt_interval = NaN;
    std_qt_interval = NaN;
    median_qt_interval = NaN;
else
    % Calculate statistics
    mean_qt_interval = mean(valid_qt_interval);
    std_qt_interval = std(valid_qt_interval);
    median_qt_interval = median(valid_qt_interval);
end

QTc_b = mean_qt_interval/ ((mean_rr_interval*convert_ms_s )^(1/2));
QTc_f = mean_qt_interval / ((mean_rr_interval*convert_ms_s )^(1/3));

%% Check if any valid intervals are present for pr_interval
if isempty(valid_pr_interval)
    mean_pr_interval = NaN;
    std_pr_interval = NaN;
    median_pr_interval = NaN;
else

    % Calculate statistics for pr_interval
    mean_pr_interval = mean(valid_pr_interval);
    std_pr_interval = std(valid_pr_interval);
    median_pr_interval = median(valid_pr_interval);
end

%% Check if any valid intervals are present for st_segment
if isempty(valid_st_segment)
    mean_st_segment= NaN;
    std_st_segment = NaN;
    median_st_segment = NaN;
else

    % Calculate statistics for st_interval
    mean_st_segment = mean(valid_st_segment);
    std_st_segment = std(valid_st_segment);
    median_st_segment = median(valid_st_segment);
end
%% Check if any valid intervals are present for pr_segment
if isempty(valid_pr_segment)
    mean_pr_segment = NaN;
    std_pr_segment = NaN;
    median_pr_segment = NaN;
else

    % Calculate statistics for pr_segment
    mean_pr_segment = mean(valid_pr_segment);
    std_pr_segment = std(valid_pr_segment);
    median_pr_segment = median(valid_pr_segment);
end

%% Check if any valid intervals are present for valid_rt_peaks
if isempty(valid_rt_peaks)
    mean_rt_peaks = NaN;
    std_rt_peaks = NaN;
    median_rt_peaks = NaN;
else

    % Calculate statistics for valid_rt_peaks
    mean_rt_peaks = mean(valid_rt_peaks);
    std_rt_peaks = std(valid_rt_peaks);
    median_rt_peaks = median(valid_rt_peaks);
end

%% Check if any valid intervals are present for valid_pr_peaks
if isempty(valid_pr_peaks)
    mean_pr_peaks = NaN;
    std_pr_peaks = NaN;
    median_pr_peaks = NaN;
else

    % Calculate statistics for valid_pr_peaks
    mean_pr_peaks = mean(valid_pr_peaks);
    std_pr_peaks = std(valid_pr_peaks);
    median_pr_peaks = median(valid_pr_peaks);
end

%% Check if any valid intervals are present for valid_qr_peaks
if isempty(valid_qr_peaks)
    mean_qr_peaks = NaN;
    std_qr_peaks = NaN;
    median_qr_peaks = NaN;
else

    % Calculate statistics for valid_qr_peaks
    mean_qr_peaks = mean(valid_qr_peaks);
    std_qr_peaks = std(valid_qr_peaks);
    median_qr_peaks = median(valid_qr_peaks);
end

%% Check if any valid intervals are present for valid_rs_peaks
if isempty(valid_rs_peaks)
    mean_rs_peaks = NaN;
    std_rs_peaks = NaN;
    median_rs_peaks = NaN;
else

    % Calculate statistics for valid_rs_interval_peaks
    mean_rs_peaks = mean(valid_rs_peaks);
    std_rs_peaks = std(valid_rs_peaks);
    median_rs_peaks = median(valid_rs_peaks);
end


%%  Results
% % Store pr_interval results
% featureset.mean_pr_interval = mean_pr_interval;
% featureset.std_pr_interval = std_pr_interval;
% featureset.median_pr_interval = median_pr_interval;
% 
% % Store qt_interval results
% featureset.mean_qt_interval = mean_qt_interval;
% featureset.std_qt_interval = std_qt_interval;
% featureset.median_qt_interval = median_qt_interval;
% featureset.QTc_b = QTc_b;
% featureset.QTc_f = QTc_f;
% 
% % Store st_segment results
% featureset.mean_st_segment = mean_st_segment;
% featureset.std_st_segment = std_st_segment;
% featureset.median_st_segment = median_st_segment;
% 
% % Store pr_segment results
% featureset.mean_pr_segment = mean_pr_segment;
% featureset.std_pr_segment = std_pr_segment;
% featureset.median_pr_segment = median_pr_segment;
% 
% % Store pr_peaks results
% featureset.mean_pr_peaks = mean_pr_peaks;
% featureset.std_pr_peaks = std_pr_peaks;
% featureset.median_pr_peaks = median_pr_peaks;
% 
% % Store qr_peaks results
% featureset.mean_qr_peaks = mean_qr_peaks;
% featureset.std_qr_peaks = std_qr_peaks;
% featureset.median_qr_peaks = median_qr_peaks;
% 
% % Store qrs_complex results
% featureset.mean_qrs_complex = mean_qrs;
% featureset.std_qrs_complex = std_qrs;
% featureset.median_qrs_complex = median_qrs;
% 
% % Store rs_peaks results
% featureset.mean_rs_peaks = mean_rs_peaks;
% featureset.std_rs_peaks = std_rs_peaks;
% featureset.median_rs_peaks = median_rs_peaks;
% 
% % Store rt_peaks results
% featureset.mean_rt_peaks = mean_rt_peaks;
% featureset.std_rt_peaks = std_rt_peaks;
% featureset.median_rt_peaks = median_rt_peaks;

feature_vec = [mean_pr_interval, std_pr_interval, median_pr_interval, mean_qt_interval, std_qt_interval, median_qt_interval, QTc_b, QTc_f, ...
    mean_st_segment, std_st_segment, median_st_segment, mean_pr_segment, std_pr_segment, median_pr_segment, ...
    mean_pr_peaks, std_pr_peaks, median_pr_peaks, mean_qr_peaks, std_qr_peaks, median_qr_peaks, mean_qrs, std_qrs, median_qrs, ...
    mean_rs_peaks, std_rs_peaks, median_rs_peaks, mean_rt_peaks, std_rt_peaks, median_rt_peaks];

% Define feature info
feature_info.names = {'mean_pr_interval', 'std_pr_interval', 'median_pr_interval', 'mean_qt_interval', 'std_qt_interval', 'median_qt_interval', 'qtc_b','qtc_f' ...
    'mean_st_segment', 'std_st_segment', 'median_st_segment', 'mean_pr_segment', 'std_pr_segment', 'median_pr_segment', ...
    'mean_pr_peaks_interval', 'std_pr_peaks_interval', 'median_pr_peaks_interval', 'mean_qr_peaks_interval', 'std_qr_peaks_interval', 'median_qr_peaks_interval', ...
    'mean_qrs_complex_interval', 'std_qrs_complex_interval', 'median_qrs_complex_interval', ...
    'mean_rs_peaks_interval', 'std_rs_peaks_interval', 'median_rs_peaks_interval', 'mean_rt_peaks_interval', 'std_rt_peaks_interval', 'median_rt_peaks_interval'};

feature_info.units = repmat({'ms'}, 1, length(feature_info.names));

feature_info.description = {"Mean PR-interval", "Standard deviation PR-interval", "Median PR-interval", "Mean QT-interval", "Standard deviation QT-interval", "Median QT-interval", "QTc (Bazett)","QTc (Fridericia)" ...
    "Mean ST-segment (QRSoff-Ton)", "Standard deviation ST-segment (QRSoff-Ton)", "Median ST-segment (QRSoff-Ton)", "Mean PR-segment (Poff-QRSon)", "Standard deviation PR-segment (Poff-QRSon)", "Median PR-segment (Poff-QRSon)", ...
    "Mean P to R peaks interval", "Standard deviation P to R peaks interval", "Median P to R peaks interval", "Mean QRSon to R peak interval", "Standard deviation QRSon to R peak interval", "Median QRSon to R peak interval", ...
    "Mean  QRS-complex interval", "Standard deviation QRS-complex interval", "Median  QRS-complex interval", ...
    "Mean R peak to QRSoff interval", "Standard deviation R peak to QRSoff interval", "Median R peak to QRSoff interval", "Mean R peak to T peak interval", "Standard deviation R peak to T peak interval", "Median R peak to T peak interval"};

end
