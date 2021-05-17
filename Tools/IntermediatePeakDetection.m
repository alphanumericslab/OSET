function interpeaks = IntermediatePeakDetection(x, peaks, wlen, num, varargin)
%
% interpeaks = IntermediatePeakDetection(x, peaks, wlen, num, flag),
% Searches for intermediate peaks between a set of readily found landmark
% peaks. Useful for detecting P and T wave peaks from the reference R-peaks
% or for PCG S1 and S2 detection from ECG R-peaks
%
% inputs:
% x: vector of input data
% peaks: marker peaks
% wlen: running search window length
% num: number of intermediate peaks to report
% flag: search for positive (flag=1) or negative (flag=0) peaks. By default
% the maximum absolute value of the signal, determines the peak sign.
%
% output:
% interpeaks: vector of intermediate peaks impulse train
%
% Notes:
% - The signal baseline wander should be removed before the R-peak detection
%
%
% Open Source Electrophysiological Toolbox, version 3.14, May 2021
% Released under the GNU General Public License
% Copyright (C) 2021  Reza Sameni
% Emory University
% reza.sameni@gmail.com

if nargin > 4
    flag = varargin{1};
else
    flag = abs(max(x)) > abs(min(x));
end

if nargin > 5
    method = varargin{2};
else
    method = 1; % use second method by default
end

if(method == 0) % method 0: find the dominant 'num' local peaks between refs
    N = length(x);
    interpeaks = zeros(1, N);
    
    I_refpeaks = find(peaks);
    if(I_refpeaks(1) > 1)
        I_refpeaks = [1, I_refpeaks];
    end
    if(I_refpeaks(end) < N)
        I_refpeaks = [I_refpeaks, N];
    end
    
    for j = 1 : length(I_refpeaks) - 1
        ref_ind_range = I_refpeaks(j) : I_refpeaks(j + 1) - 1;
        x_segment = x(ref_ind_range);
        interpeaks_tmp = PeakDetection(x_segment, 1 / wlen, flag);
        JJ = find(interpeaks_tmp);
        I_amps = x_segment(JJ);
        if(flag)
            [~, I] = sort(I_amps, 'descend');
        else
            [~, I] = sort(I_amps, 'ascend');
        end
        interpeaks(ref_ind_range(JJ(I(1 : num)))) = 1;
    end
    
elseif(method == 1) % method 1: find the dominant 'num' local peaks between refs (second appaorach)
    N = length(x);
    I_refpeaks = find(peaks);
    if(I_refpeaks(1) > 1)
        I_refpeaks = [1, I_refpeaks];
    end
    if(I_refpeaks(end) < N)
        I_refpeaks = [I_refpeaks, N];
    end
    
    all_peaks = PeakDetection(x, 1 / wlen, flag);
    
    interpeaks = all_peaks;
    for j = 1 : length(I_refpeaks) - 1
        ref_ind_range = I_refpeaks(j) : I_refpeaks(j + 1) - 1;
        x_segment = x(ref_ind_range);
        segment_peaks = all_peaks(ref_ind_range);
        I_segment_peaks = find(segment_peaks);
        segment_peaks_amps = x_segment(I_segment_peaks);
        if(flag)
            [~, I] = sort(segment_peaks_amps, 'descend');
        else
            [~, I] = sort(segment_peaks_amps, 'ascend');
        end
        I_absolute_poistion = ref_ind_range(I_segment_peaks(I));
        if(length(I_absolute_poistion) > num)
            interpeaks(I_absolute_poistion(num + 1 : end)) = 0;
        end
    end
elseif(method == 2) % method 2: find the first num local peak after each ref in a window of length wlen
    N = length(x);
    I_refpeaks = find(peaks);
    if(I_refpeaks(1) > 1)
        I_refpeaks = [1, I_refpeaks];
    end
    if(I_refpeaks(end) < N)
        I_refpeaks = [I_refpeaks, N];
    end
    
    all_peaks = PeakDetection(x, 1 / wlen, flag);
    
    interpeaks = zeros(1, N);
    for j = 1 : length(I_refpeaks) - 1
        ref_ind_range = I_refpeaks(j) : min(I_refpeaks(j) + wlen - 1, N);
        x_segment = x(ref_ind_range);
        segment_peaks = all_peaks(ref_ind_range);
        I_segment_peaks = find(segment_peaks);
        segment_peaks_amps = x_segment(I_segment_peaks);
        if(flag)
            [~, I] = sort(segment_peaks_amps, 'descend');
        else
            [~, I] = sort(segment_peaks_amps, 'ascend');
        end
        I_absolute_poistion = ref_ind_range(I_segment_peaks(I));
        if(~isempty(I_absolute_poistion))
            interpeaks(I_absolute_poistion(1 : num)) = 1;
        end
    end
else
    error('Unknown method');
end