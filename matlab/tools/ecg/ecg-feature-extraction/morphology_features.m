function sampled_values = morphology_features(mean_beat, num_samples, norm_flag)
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
%   Email: bmemousavi@gmail.com
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

end