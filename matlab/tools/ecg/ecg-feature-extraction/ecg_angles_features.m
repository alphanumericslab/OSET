function [feature_vec, feature_info] = ecg_angles_features(data, position, fs)
% featureset = ecg_angles_features(data, position, fs)
% Extract features related to the ECG amplitude-to-interval ratios
%
% Inputs:
%   data: ECG signal (1D array).
%   position: fiducial points of ECG signal
%   fs: Sampling frequency (Hz)
%
% Output:
%   featureset: Structure containing features related to ECG amplitude-to-interval ratios
%   Units: mv/ms (millivolts per millisecond)
%
% Author:
%   Seyedeh Somayyeh Mousavi
%   Emory University, Georgia, USA
%   Email: bmemousavi@gmail.com
%   Date: SEP 24, 2024
% Author: Sajjad Karimi
% Location: Emory University, Georgia, USA
% Email: sajjadkarimi91@gmail.com
% Date: Apr 10, 2025

%% Constant value
convert_s_ms =1000;
%% Define ECG points
p = position.P;
qrs_onset = position.QRSon;
rpeak = position.R;
qrs_offset = position.QRSoff;
t = position.T;

%% Calculate intervals
rt_interval = convert_s_ms * ( t - rpeak) / fs;
pr_interval = convert_s_ms * (rpeak - p) / fs;
qr_interval = convert_s_ms * (rpeak - qrs_onset) / fs;
rs_interval = convert_s_ms * (qrs_offset - rpeak) / fs;

% Initialize rt_amp with NaNs of the same size as rpeak
rt_amp = NaN * ones(size(rpeak));
% Loop over each index and compute the difference
for i = 1:length(rpeak)
    if ~isnan(rpeak(i)) && ~isnan(t(i))
        rt_amp(i) = data(rpeak(i)) - data(t(i));
    end
end

% Initialize pr_amp with NaNs of the same size as rpeak
pr_amp = NaN * ones(size(rpeak));
% Loop over each index and compute the difference
for i = 1:length(rpeak)
    if ~isnan(rpeak(i)) && ~isnan(p(i))
        pr_amp(i) = data(rpeak(i)) - data(p(i));
    end
end

% Initialize qr_amp with NaNs of the same size as rpeak
qr_amp = NaN * ones(size(rpeak));
% Loop over each index and compute the difference
for i = 1:length(rpeak)
    if ~isnan(rpeak(i)) && ~isnan(qrs_onset(i))
        qr_amp(i) = data(rpeak(i)) - data(qrs_onset(i));
    end
end

% Initialize rs_amp with NaNs of the same size as rpeak
rs_amp = NaN * ones(size(rpeak));
% Loop over each index and compute the difference
for i = 1:length(rpeak)
    if ~isnan(rpeak(i)) && ~isnan(qrs_offset(i))
        rs_amp(i) = data(rpeak(i)) - data(qrs_offset(i));
    end
end

% Calculate ratios
ratio_rt = rt_amp./ rt_interval;
ratio_pr = pr_amp./ pr_interval;
ratio_qr = qr_amp./ qr_interval;
ratio_rs = rs_amp./ rs_interval;

% Omit NaN values
valid_ratio_rt = ratio_rt(~isnan(ratio_rt));
valid_ratio_pr = ratio_pr(~isnan(ratio_pr));
valid_ratio_qr = ratio_qr(~isnan(ratio_qr));
valid_ratio_rs = ratio_rs(~isnan(ratio_rs));

%% Check if any valid values are present for valid_ratio_rt
if isempty(valid_ratio_rt)
    mean_rt_ratio = NaN;
    std_rt_ratio = NaN;
    median_rt_ratio = NaN;
else
    mean_rt_ratio = mean(valid_ratio_rt);
    std_rt_ratio = std(valid_ratio_rt);
    median_rt_ratio = median(valid_ratio_rt);
end

%% Check if any valid values are present for valid_ratio_pr
if isempty(valid_ratio_pr)
    mean_pr_ratio = NaN;
    std_pr_ratio = NaN;
    median_pr_ratio = NaN;
else
    mean_pr_ratio = mean(valid_ratio_pr);
    std_pr_ratio = std(valid_ratio_pr);
    median_pr_ratio = median(valid_ratio_pr);
end

%% Check if any valid values are present for valid_ratio_qr
if isempty(valid_ratio_qr)
    mean_qr_ratio = NaN;
    std_qr_ratio = NaN;
    median_qr_ratio = NaN;
else
    mean_qr_ratio = mean(valid_ratio_qr);
    std_qr_ratio = std(valid_ratio_qr);
    median_qr_ratio = median(valid_ratio_qr);
end

%% Check if any valid values are present for valid_ratio_rs
if isempty(valid_ratio_rs)
    mean_rs_ratio = NaN;
    std_rs_ratio = NaN;
    median_rs_ratio = NaN;
else
    mean_rs_ratio = mean(valid_ratio_rs);
    std_rs_ratio = std(valid_ratio_rs);
    median_rs_ratio = median(valid_ratio_rs);
end

%% Results
feature_vec = [mean_pr_ratio, std_pr_ratio, median_pr_ratio, mean_qr_ratio, std_qr_ratio, median_qr_ratio,...
    mean_rs_ratio, std_rs_ratio, median_rs_ratio, mean_rt_ratio, std_rt_ratio, median_rt_ratio];

% Define feature info
names = {    'mean_pr_ratio', 'std_pr_ratio', 'median_pr_ratio', ...
    'mean_qr_ratio', 'std_qr_ratio', 'median_qr_ratio', ...
    'mean_rs_ratio', 'std_rs_ratio', 'median_rs_ratio', ...
    'mean_rt_ratio', 'std_rt_ratio', 'median_rt_ratio'};

units = repmat({'mv/ms'}, 1, length(names));
description = {"Mean P-R amplitude-to-interval ratio", "Standard deviation P-R amplitude-to-interval ratio", "Median P-R amplitude-to-interval ratio", ...
    "Mean Q-R amplitude-to-interval ratio", "Standard deviation Q-R amplitude-to-interval ratio", "Median Q-R amplitude-to-interval ratio", ...
    "Mean R-S amplitude-to-interval ratio", "Standard deviation R-S amplitude-to-interval ratio", "Median R-S amplitude-to-interval ratio", ...
    "Mean R-T amplitude-to-interval ratio", "Standard deviation R-T amplitude-to-interval ratio", "Median R-T amplitude-to-interval ratio"};
feature_info = struct('names', {names}, 'units', {units}, 'description', {description} );

end
