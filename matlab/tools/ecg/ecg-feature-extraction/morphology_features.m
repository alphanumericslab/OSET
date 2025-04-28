function [feature_vec, feature_info] = morphology_features(mean_beat, num_samples, norm_flag)
% featureset = morphology_features(mean_beat, num_samples, norm_flag)
% Extract features related to the morphology of the average ECG beat
%
% Inputs:
%   mean_beat: Average ECG beat signal (1D array)
%   num_samples: Number of samples to extract from the mean beat
%   norm_flag: Boolean flag to indicate whether to normalize the samples
%              (1: normalize, 0: do not normalize)
%
% Output:
%   feature_vec: Vector containing the sampled values from the mean beat
%   feature_info: Structure containing feature information with fields:
%       - names: Cell array of feature names
%       - units: Cell array of feature units ('mv' or 'scaler')
%       - description: Cell array of feature descriptions
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