function [interpeaks, ref_to_out_lag] = intermediate_peak_det(x, peaks, wlen, n_peaks, varargin)
% intermediate_peak_det - Searches for intermediate peaks between a
%   set of readily found landmark peaks. Useful for detecting P and T wave
%   peaks from the reference R-peaks or for PCG S1 and S2 detection from ECG
%   R-peaks
% 
% Syntax:
%   [interpeaks, ref_to_out_lag] = intermediate_peak_det(x, peaks, wlen, n_peaks, flag, method)
%
% Inputs:
%   x: Vector of input data.
%   peaks: Marker peaks indicating the landmarks (e.g., ECG R-peaks).
%   wlen: Running search window length.
%   n_peaks: Number of intermediate peaks to report.
%   flag (optional): - Search for positive (flag=1) or negative (flag=0) peaks. By default,
%       the maximum absolute value of the signal determines the peak sign.
%   method (optional): Model used for calculating the intermediate peaks (see code for functionality).
%       method = 1 by default.
%
% Outputs:
%   interpeaks: Vector of intermediate peaks impulse train.
%   ref_to_out_lag: Lag between input reference peaks and output peaks (in method = 2 only).
%
% Notes:
% - The signal baseline wander should be removed before the R-peak detection.
%
% Revision History:
%   2021: First release
%   2023: Renamed from deprecated version IntermediatePeakDetection
%
% Reza Sameni, 2021-2023
% The Open-Source Electrophysiological Toolbox

% Handle optional input arguments
if nargin > 4
    flag = varargin{1};
else
    flag = abs(max(x)) > abs(min(x)); % Default: Search for the maximum absolute value
end

if nargin > 5
    method = varargin{2};
else
    method = 1; % Default: Use the second method
end


N = length(x);
I_refpeaks = find(peaks);
if I_refpeaks(1) > 1
    I_refpeaks = [1, I_refpeaks];
end
if I_refpeaks(end) < N
    I_refpeaks = [I_refpeaks, N];
end

% Find the dominant 'n_peaks' local peaks between reference peaks
if method == 0 % Method 0: start from an empty set and add peaks one-by-one
    ref_to_out_lag = [];
    interpeaks = zeros(1, N);
    for j = 1 : length(I_refpeaks) - 1
        ref_ind_range = I_refpeaks(j) : I_refpeaks(j + 1) - 1;
        x_segment = x(ref_ind_range);
        interpeaks_tmp = peak_det_local_search(x_segment, 1 / wlen, flag); % Find all potential local peaks
        JJ = find(interpeaks_tmp);
        I_amps = x_segment(JJ);
        if(flag)
            [~, I] = sort(I_amps, 'descend'); % Sort amplitudes in descending order
        else
            [~, I] = sort(I_amps, 'ascend'); % Sort amplitudes in ascending order
        end
        interpeaks(ref_ind_range(JJ(I(1 : n_peaks)))) = 1; % Mark the dominant 'n_peaks' local peaks between references
    end
    
elseif method == 1 % Method 1: start from the superset of points and remove the less probable ones one-by-one
    all_peaks = peak_det_local_search(x, 1 / wlen, flag); % Find all potential local peaks
    
    interpeaks = all_peaks;
    ref_to_out_lag = zeros(1, length(I_refpeaks));
    for j = 1 : length(I_refpeaks) - 1
        ref_ind_range = I_refpeaks(j) : I_refpeaks(j + 1) - 1;
        x_segment = x(ref_ind_range);
        segment_peaks = all_peaks(ref_ind_range);
        I_segment_peaks = find(segment_peaks);
        segment_peaks_amps = x_segment(I_segment_peaks);
        if flag
            [~, I] = sort(segment_peaks_amps, 'descend'); % Sort amplitudes in descending order
        else
            [~, I] = sort(segment_peaks_amps, 'ascend'); % Sort amplitudes in ascending order
        end
        I_absolute_poistion = ref_ind_range(I_segment_peaks(I));
        if length(I_absolute_poistion) > n_peaks
            interpeaks(I_absolute_poistion(n_peaks + 1 : end)) = 0; % Remove less dominant peaks
        end
    end
    
elseif method == 2 % Method 2: Find the first 'n_peaks' local peaks after each reference in a window of length 'wlen'
    N = length(x);
    all_peaks = peak_det_local_search(x, 1 / wlen, flag); % Call peak_det_local_search function
    
    interpeaks = zeros(1, N);
    ref_to_out_lag = zeros(1, length(I_refpeaks));
    for j = 1 : length(I_refpeaks)
        ref_ind_range = I_refpeaks(j) : min(I_refpeaks(j) + wlen - 1, N); % Define the search window
        x_segment = x(ref_ind_range);
        segment_peaks = all_peaks(ref_ind_range);
        I_segment_peaks = find(segment_peaks);
        segment_peaks_amps = x_segment(I_segment_peaks);
        if flag
            [~, I] = sort(segment_peaks_amps, 'descend'); % Sort amplitudes in descending order
        else
            [~, I] = sort(segment_peaks_amps, 'ascend'); % Sort amplitudes in ascending order
        end
        I_absolute_poistion = ref_ind_range(I_segment_peaks(I));
        if(~isempty(I_absolute_poistion))
            interpeaks(I_absolute_poistion(1 : n_peaks)) = 1; % Mark the first 'n_peaks' local peaks after each reference
            ref_to_out_lag(j) = I_absolute_poistion(1) - I_refpeaks(j); % Calculate the lag between references and output peaks
        end
    end
else
    error('Unknown method');
end
