function [feature_vec, feature_info] = morphology_features(mean_beat, num_samples, norm_flag)
% feature = ecg_complexity_features(data)
% Extract features related to the complexity analysis of ECG signal
%
% Inputs:
%   data: ECG signal (1D array).
%
% Output:
%   feature: Structure containing two values (Mobility and Complexity)
%
% Author:
%   Sajjad Karimi
%   Emory University, Georgia, USA
%   Email: sajjadkarimi91@gmail.com
%   Date: Mar 11, 2025

T = length(mean_beat);
N = num_samples; % Number of samples to extract

% Compute equally spaced indices from 1 to T
sample_indices = round(linspace(1, T, N));

% Extract the corresponding values
sampled_values = mean_beat(sample_indices);

if norm_flag
    sampled_values = sampled_values/norm(sampled_values);
end

feature_vec = sampled_values;

% Define feature info
for n = 1:num_samples

    if norm_flag
        feature_info.names{n} = ['norm_morphology_sample_', num2str(n)];
        feature_info.units{n} = 'scaler';
        feature_info.description{n} = string(['Normalized average beat sample ', num2str(n)]);
    else
        feature_info.names{n} = ['morphology_sample_', num2str(n)];
        feature_info.units{n} = 'mv';
        feature_info.description{n} = string(['Average beat sample ', num2str(n)]);
    end

end

end