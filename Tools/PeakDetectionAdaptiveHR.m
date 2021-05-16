function peaks = PeakDetectionAdaptiveHR(x, ff, fs, varargin)
%
% peaks = PeakDetectionAdaptiveHR(x, f, fs, flag),
% R-peak detector based on max search over sliding window with adaptive
% width
%
% inputs:
% x: vector of input data
% f: approximate ECG beat-rate in Hz
% fs: sampling frequency in Hz
% flag: search for positive (flag=1) or negative (flag=0) peaks. By default
% the maximum absolute value of the signal, determines the peak sign.
%
% output:
% peaks: vector of R-peak impulse train
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

N = length(x);
peaks0 = zeros(1,N);
peaks = zeros(1,N);

th = 0.5;
th2 = 0.8;
rng = floor(th/(ff/fs));

if nargin == 4
    flag = varargin{1};
else
    flag = abs(max(x)) > abs(min(x));
end

if(flag)
    % First round peak detection
    for j = 1 : N
        if j > rng && j < N - rng
            index = j-rng : j+rng;
        elseif j > rng
            index = N-2*rng : N;
        else
            index = 1 : 2*rng;
        end
        
        if max(x(index)) == x(j)
            peaks0(j) = 1;
        end
    end
    
    II = find(peaks0);
    RR_intervals = diff(II);
    RR_intervals_smoothed = BaseLine1(RR_intervals, 3, 'md');
    RR_intervals_smoothed2 = BaseLine1(RR_intervals_smoothed, 3, 'mn');
    ff = fs ./ [RR_intervals_smoothed2(1) RR_intervals_smoothed2];
    if(II(1) > 1)
        II = [1, II];
        ff = [mean(ff(1:2)), ff];
    end
    if(II(end) < N)
        II = [II, N];
        ff = [ff, mean(ff(end-1:end))];
    end
    ff_interpolated = interp1(II, ff, 1 : N);
    rng2 = floor(th2./(ff_interpolated/fs));
    
    % Second round peak detection
    for j = 1 : N
        if j > rng2(j) && j < N-rng2(j)
            index = j-rng2(j) : j+rng2(j);
        elseif j > rng2(j)
            index = N-2*rng2(j) : N;
        else
            index = 1 : 2*rng2(j);
        end
        
        if max(x(index)) == x(j)
            peaks(j) = 1;
        end
    end
else
    % First round peak detection
    for j = 1 : N
        if j > rng && j < N-rng
            index = j-rng : j+rng;
        elseif j>rng
            index = N-2*rng : N;
        else
            index = 1 : 2*rng;
        end
        
        if min(x(index))==x(j)
            peaks0(j) = 1;
        end
    end
    
    II = find(peaks0);
    RR_intervals = diff(II);
    RR_intervals_smoothed = BaseLine1(RR_intervals, 3, 'md');
    RR_intervals_smoothed2 = BaseLine1(RR_intervals_smoothed, 3, 'mn');
    ff = fs ./ [RR_intervals_smoothed2(1) RR_intervals_smoothed2];
    if(II(1) > 1)
        II = [1, II];
        ff = [mean(ff(1:2)), ff];
    end
    if(II(end) < N)
        II = [II, N];
        ff = [ff, mean(ff(end-1:end))];
    end
    ff_interpolated = interp1(II, ff, 1 : N);
    rng2 = floor(th2./(ff_interpolated/fs));
    
    % Second round peak detection
    for j = 1 : N
        if j > rng2(j) && j < N-rng2(j)
            index = j-rng2(j) : j+rng2(j);
        elseif j > rng2(j)
            index = N-2*rng2(j) : N;
        else
            index = 1 : 2*rng2(j);
        end
        
        if min(x(index)) == x(j)
            peaks(j) = 1;
        end
    end
end

% remove fake peaks (below 25% of the median peak amplitudes)
I = find(peaks);
peak_amps = median(x(I));
J = abs(x(I)) < 0.3 * abs(peak_amps);
peaks(I(J)) = 0;
% d = diff(I);
% z = find(d < rng);
% peaks(I(d < rng))=0;