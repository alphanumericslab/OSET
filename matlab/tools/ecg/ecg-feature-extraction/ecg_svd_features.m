function featureset = ecg_svd_features(data, R_peak_indexes, number_eigenvalues)
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

    % Step 1: Compute RR intervals in samples
    RR_intervals_samples = diff(R_peak_indexes);

    % Step 2: Determine the event bounds (median of RR intervals)
    event_bounds = round(median(RR_intervals_samples));

    % Ensure event_bounds is odd
    if mod(event_bounds, 2) == 0
        event_bounds = event_bounds + 1;
    end

    % Step 3: Extract signal segments using event_stacker
    [stacked_events, ~] = event_stacker(data, R_peak_indexes, event_bounds);

    % Step 4: Perform Singular Value Decomposition (SVD)
    singular_values = svd(stacked_events);

    % Step 5: Normalize singular values
    singular_values = singular_values / sum(singular_values);

    % Store the normalized singular values into the output array
    if length(singular_values) >= number_eigenvalues
        % Truncate if there are more eigenvalues than desired
        singularvalues_normalized = singular_values(1:number_eigenvalues);
    elseif length(singular_values) < number_eigenvalues
        % Zero-pad if there are fewer eigenvalues than desired
        singularvalues_normalized = [singular_values; zeros(number_eigenvalues - length(singular_values), 1)];
    end

    % Step 6: Convert normalized singular values to percentages
    featureset = singularvalues_normalized' * 100;
end

