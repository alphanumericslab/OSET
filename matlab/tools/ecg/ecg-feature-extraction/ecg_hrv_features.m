function [feature_vec, feature_info]  = ecg_hrv_features(rpeak_indexes, fs)
%%
% Matlab function for extracting Heart Rate Variability (HRV) features.
% based on time units (seconds).
%
% INPUT:
% R_peaks - A vector containing the R-peak indices of the ECG signal (expressed as sample points).
% fs      - Sampling frequency (Hz), used to convert sample indices to time units.
%
% OUTPUT:
% features - A structure containing the following HRV features:
%   1. RMSSD (Root Mean Square of Successive Differences) in milliseconds
%   2. SDNN (Standard Deviation of normal-to-normal (NN) intervals) in milliseconds
%   3. HR (MEDIAN, MEAN, upper_5, lower_5)
%
% Author: Seyedeh Somayyeh Mousavi
% Date: Dec 17, 2024
% Location: Emory University, Georgia, USA
% Email: bmemousavi@gmail.com
% Author: Sajjad Karimi
% Location: Emory University, Georgia, USA
% Email: sajjadkarimi91@gmail.com
% Date: Mar 14, 2025
%=======================================================================
%% Constant value
convert_s_ms =1000;

% Extract features if length(R_peaks) > 3
if length(rpeak_indexes) > 3
        
    % Compute RR intervals in samples
    RR_intervals_samples = diff(rpeak_indexes); % Time between successive R-peaks in samples

    % Convert RR intervals to time in seconds using sampling frequency (fs)
    RR_intervals_seconds = RR_intervals_samples / fs;
    
    % Filtering for non-normal beats
    RR_intervals_seconds(RR_intervals_seconds<0.2) = [];
    RR_intervals_seconds(RR_intervals_seconds>2) = [];

    % RR intervals [2.5 97.5] percentile for filtering
    IQR_RR = prctile(RR_intervals_seconds, [2.5 97.5]);
    filtered_RR = RR_intervals_seconds(RR_intervals_seconds >= IQR_RR(1) & RR_intervals_seconds <= IQR_RR(2));

    % Feature 1: RMSSD (Root Mean Square of Successive Differences)
    RMSSD = sqrt(mean(diff(filtered_RR).^2));
    % Convert to milliseconds
    RMSSD = RMSSD * convert_s_ms;

    % Feature 2: SDNN (Standard Deviation of NN intervals)
    SDNN = std(filtered_RR);
    % Convert to milliseconds
    SDNN = SDNN * convert_s_ms;

    % Calculate heart rate (HR)
    HR = 60 ./ RR_intervals_seconds;  % HR in beats per minute based on RR intervals in seconds

    % Feature 3: Calculate median heart rate
    median_HR = median(HR);

    % Feature 4: Calculate heart rate for lower 5% 
    HR_lower_5 = prctile(HR, 5);

    % Feature 5: Calculate heart rate for upper 5% 
    HR_upper_5 = prctile(HR, 95);

    % Feature 6: Calculate heart rate for the interquartile range (5%-95%)
    IQR_HR = prctile(HR, [5 95]);
    HR = HR(HR >= IQR_HR(1) & HR <= IQR_HR(2));  % Filter HR values within the range
    HR_mean = mean(HR);

    % Store the HRV features
    features.N_beats = length(rpeak_indexes);
    features.RMSSD = RMSSD;  
    features.SDNN = SDNN;
    features.HR_median = median_HR;
    features.HR_mean = HR_mean;
    features.HR_upper_5 = HR_upper_5;
    features.HR_lower_5 = HR_lower_5;
    
else
    % Return NaNs if insufficient R-peaks
    features.N_beats = length(rpeak_indexes);
    features.RMSSD = nan;
    features.SDNN = nan;
    features.HR_median = nan;
    features.HR_mean = nan;
    features.HR_upper_5 = nan;
    features.HR_lower_5 = nan;
    
end

feature_vec = [features.N_beats, features.RMSSD, features.SDNN, features.HR_median, features.HR_mean, features.HR_upper_5, features.HR_lower_5 ];
% Define feature info
feature_info.names = {'n_beats', 'rmssd', 'sdnn', 'hr_median', 'hr_mean', 'hr_upper_5', 'hr_lower_5'};
feature_info.units = {'scaler', 'ms', 'ms', 'bpm', 'bpm', 'bpm', 'bpm'};
feature_info.description = {"Number of ECG beat", "Root Mean Square of Successive Differences for HRV", ...
    "Standard Deviation of Normal-to-Normal Intervals", "Median heart-rate", "Mean heart-rate", "95% heart-rate", "5% heart-rate"};

end
