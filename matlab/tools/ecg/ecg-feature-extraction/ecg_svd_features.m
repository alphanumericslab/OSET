function [feature_vec, feature_info] = ecg_svd_features(data, rpeak_indexes, n_svd)
% featureset = ecg_svd_features(data, R_peak_indexes, number_eigenvalues)
% Extract features from ECG using Singular Value Decomposition (SVD)
%
% Inputs:
%   data: ECG signal (1D array).
%   R_peaks_indexes - A vector containing the R-peak indices of the ECG signal (expressed as sample points).
%   number_eigenvalues: Number of eigenvalues to consider for the feature vector (scaler).
%
% Output:
%   featureset: Structure containing normalized singular values in percentage.
%
% Dependencies:
%   `events_snr` function from the OSET package
%
% Author:
%   Seyedeh Somayyeh Mousavi
%   Emory University, Georgia, USA
%   Email: bmemousavi@gmail.com
%   Date: OCT 14, 2024
% Author:
%   Sajjad Karimi
%   Emory University, Georgia, USA
%   Email: sajjadkarimi91@gmail.com
%   Date: Mar 14, 2025

% Step 1: Compute RR intervals in samples
RR_intervals_samples = diff(rpeak_indexes);

% Step 2: Determine the event bounds (median of RR intervals)
event_bounds = round(median(RR_intervals_samples));

% Ensure event_bounds is odd
if mod(event_bounds, 2) == 0
    event_bounds = event_bounds + 1;
end

% Step 3: Extract signal segments using event_stacker
[stacked_events, ~] = event_stacker(data, rpeak_indexes, event_bounds);

% Step 4: Perform Singular Value Decomposition (SVD)
singular_values = svd(stacked_events);

% Step 5: Normalize singular values
singular_values = singular_values / sum(singular_values);

% Store the normalized singular values into the output array
if length(singular_values) >= n_svd
    % Truncate if there are more eigenvalues than desired
    singularvalues_normalized = singular_values(1:n_svd);
elseif length(singular_values) < n_svd
    % Zero-pad if there are fewer eigenvalues than desired
    singularvalues_normalized = [singular_values; zeros(n_svd - length(singular_values), 1)];
end

% Step 6: Convert normalized singular values to percentages
feature_vec = singularvalues_normalized' * 100;

% Define feature info
for n = 1:n_svd
    feature_info.names{n} = ['svd', num2str(n)];
    feature_info.units{n} = 'scaler';
    feature_info.description{n} = string(['Normalized SVD value ', num2str(n)]);
end



end

