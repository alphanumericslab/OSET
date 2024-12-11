function [data_out, data_out_smoothed] = conditional_median_filter(data_in, wlen, threshold, fc)
% CONDITIONAL_MEDIAN_FILTER Applies a conditional median filter to the input data
%
%   [data_out, data_out_smoothed] = CONDITIONAL_MEDIAN_FILTER(data_in, wlen, threshold, fc)
%   performs a conditional median filter on the input data 'data_in'. The function replaces
%   data points that deviate from the sliding window median by more than a defined threshold.
%   After the conditional median filtering, the result is also passed through a low-pass
%   filter with zero phase.
%
% INPUTS:
%   data_in   - Input data (1D array); may contain outliers or missing values (NaNs)
%   wlen      - Length of the sliding window for median filtering. Must be odd.
%   threshold - Threshold for detecting outliers. Points deviating more than this threshold
%               from the median (in abs) within the sliding window are considered outliers.
%   fc        - Cutoff frequency for the low-pass filter applied after
%               median filtering normalized by the sampling frequency
%
% OUTPUTS:
%   data_out           - Data after applying the conditional median filter, with outliers
%                        replaced by the local median.
%   data_out_smoothed  - Data after applying both the conditional median filter and the
%                        low-pass filter.
%
% STEPS:
%   1. Ensure 'wlen' is odd. If it's even, increment it by 1 and raise a warning.
%   2. Compute the sliding window median for 'data_in' using a window of length 'wlen'.
%   3. Identify outliers based on the threshold. Outliers are defined as data points that deviate
%       from the sliding median by more than 'threshold', or points that are NaN.
%   4. Replace outliers in 'data_in' with the corresponding values from the sliding median.
%   5. Replace any remaining NaN values in the result with the global median of the filtered data.
%   6. Apply a low-pass filter with normalized cutoff frequency 'fc' to the filtered data to smooth it.
%
% EXAMPLE USAGE:
%   data_in = [1, 2, 100, 4, 5, NaN, 7, 8, 200];
%   wlen = 3;
%   threshold = 20;
%   fc = 0.2;
%   [data_out, data_out_smoothed] = conditional_median_filter(data_in, wlen, threshold, fc);
%
% The output 'data_out' will have outliers replaced by local medians, and
% 'data_out_smoothed' will be a smoothed version of 'data_out'.
%
%
% Reza Sameni, Sep 2024
% The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET

% Step 1: Ensure the window length 'wlen' is odd
if mod(wlen, 2) == 0
    wlen = wlen + 1;
    warning('wlen needs to be odd, increasing it by 1')
end

% Step 2: Compute the sliding window median
data_sliding_median = movmedian(data_in, wlen, "omitnan");

% Step 3: Identify outliers based on the threshold and NaN values
I_outliers = abs(data_in - data_sliding_median) > threshold | isnan(data_in) | isnan(data_sliding_median);

% Step 4: Replace outliers in the input data with the sliding median
data_out = data_in;
data_out(I_outliers) = data_sliding_median(I_outliers);

% Step 5: Make a copy of data_out and replace any remaining NaNs with the global median of the filtered data
data_out_nanless = data_out;
data_out_nanless(isnan(data_out_nanless)) = median(data_out_nanless, "omitnan");

% Step 6: Apply a low-pass filter to smooth the copy without nans
data_out_smoothed = lp_filter_zero_phase(data_out_nanless, fc);
end


