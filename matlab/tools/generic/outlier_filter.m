function x_filtered = outlier_filter(x_raw, method, half_wlen, percentile)
% outlier_filter - Detects outlier samples with diff values beyond a given
%   percentile and replaces them with mean or median of a sliding window
%   around the sample. Other samples remain unchanged.
%
% Syntax: x_filtered = outlier_filter(x_raw, method, half_wlen, percentile)
%
% Inputs:
%   x_raw: Raw input signal (num_channels x time).
%   method: method used for sliding window averaging 'MEAN' or 'MEDIAN' (default)
%   half_wlen: Half window length of the sliding median filter
%   percentile: percentile threshold above which samples are replaced with
%       the median of the sliding window
%
% Output:
%   x_filtered: Filtered output
%
%   Revision History:
%       2023: First release
%
% Applications: outlier detection, spike detection, QRS detection, baseline
%   wander removal, and heart rate time-series smoothing
% 
%   Reza Sameni, 2023
%   The Open-Source Electrophysiological Toolbox
%   https://github.com/alphanumericslab/OSET

df = diff(x_raw, 1, 2); % first order difference on the second dimension
diff_threshold = prctile(abs(df), percentile, 2); % find the percentile edge

T = size(x_raw, 2);

x_filtered = x_raw;
for t = 1 : T
    indexes = max(1, t - half_wlen) : min(T, t + half_wlen); % find the indexes of the surrounding window
    switch method
        case 'MEAN'
            avg = mean(x_raw(:, indexes), 2);
        case 'MEDIAN'
            avg = median(x_raw(:, indexes), 2);
    end
    er = abs(x_raw(:, t) - avg);
    replace = er >= diff_threshold;
    x_filtered(replace, t) = avg(replace);
end