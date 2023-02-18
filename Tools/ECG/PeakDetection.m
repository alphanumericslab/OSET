function [peaks, peak_indexes] = PeakDetection(x,ff,varargin)
%
% peaks = PeakDetection(x,f,flag, num_rounds),
% R-peak detector based on max search
%
% inputs:
% x: vector of input data
% f: approximate ECG beat-rate in Hertz, normalized by the sampling frequency
% flag: search for positive (flag=1) or negative (flag=0) peaks. By default
% the maximum absolute value of the signal, determines the peak sign.
% num_rounds: the number of iterations to find the R-peaks, up to 3
% (everytime updating the expected R-peak rates). Default = 1 (no iterations)
%
% output:
% peaks: vector of R-peak impulse train
% peak_indexes: vector of R-peak indexes
%
% Notes:
% - The R-peaks are found from a peak search in windows of length N; where
% N corresponds to the R-peak period calculated from the given f. R-peaks
% with periods smaller than N/2 or greater than N are not detected.
% - The signal baseline wander is recommended to be removed before the
% R-peak detection
%
% Revision History:
% 2006: First release
% 2022: Added multi-iteration feature (remains backward compatible)
%
% The Open-Source Electrophysiological Toolbox
% Reza Sameni, 2006
% https://github.com/alphanumericslab/OSET
% Sharif University of Technology, Tehran, Iran -- GIPSA-Lab, INPG, Grenoble, France

if nargin > 2 && ~isempty(varargin{1})
    flag = varargin{1};
else
    flag = abs(max(x)) > abs(min(x));
end

if nargin > 3 && ~isempty(varargin{2})
    num_rounds = varargin{2};
    if num_rounds > 3
        error('Number of R-peak detections not supported')
    end
else
    num_rounds = 1;
end

[peaks, peak_indexes] = PeakDetection_internal(x, ff, flag);
if num_rounds > 1
    for k = 1 : num_rounds - 1
        rr_intervals = diff(peak_indexes);
        ff = 1.05 / median(rr_intervals); % refined heart rate (in Hz) used for R-peak detection
        [peaks, peak_indexes] = PeakDetection_internal(x, ff, flag);
    end
end
end

% The internal R-peak detector function
function [peaks, peak_indexes] = PeakDetection_internal(x, ff, flag)
N = length(x);
peaks = zeros(1,N);

th = .5;
rng = floor(th/ff);

if flag
    for j = 1 : N
        %         index = max(j-rng,1):min(j+rng,N);
        if j > rng && j < N-rng
            index = j-rng:j+rng;
        elseif j > rng
            index = N-2*rng:N;
        else
            index = 1:2*rng;
        end

        if max(x(index))==x(j)
            peaks(j) = 1;
        end
    end
else
    for j = 1 : N
        %         index = max(j-rng,1):min(j+rng,N);
        if j > rng && j < N-rng
            index = j-rng:j+rng;
        elseif j > rng
            index = N-2*rng:N;
        else
            index = 1:2*rng;
        end

        if min(x(index))==x(j)
            peaks(j) = 1;
        end
    end
end

% remove fake peaks
I = find(peaks);
d = diff(I);
peaks(I(d < rng))=0;
peak_indexes = find(peaks);

end