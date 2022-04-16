function [HR peaks HRPeriodSD] = HRCalculation4(x, f, fs, trim, wlen, method)
%
% [HR peaks HRPeriodSD] = HRCalculation4(x,f,fs,trim,wlen,method),
% Heartrate calculator. The R-peak detector is based on max search and the
% HR is the median of the HR over the data window.
% Identical to HRCalculation2, except that the maximum point of the
% original signal is returned after the initial peak detection
%
% inputs:
% x: vector of input data
% f: approximate ECG beat-rate in Hz
% fs: sampling frequency in Hz
% trim: number of beats to trim from averaging
% wlen: window length used for positive local peak detection
% method: 'trmean' or 'median'
%
% output:
% HR: the trimmed mean/median of the HR over the data window.
% peaks: the detected R-peaks
% HRPeriodSD: the standard deviation of the RR-intervals (in samples)
%
% Notes:
% - The HR is reported as number of beats per minute (BPM)
% - The R-peaks are found from a peak search in windows of length N; where
% N corresponds to the R-peak period calculated from the given f. R-peaks
% with periods smaller than N/2 or greater than N are not detected.
% - The signal baseline wander is recommended to be removed before the
% R-peak detection
%
%
% Open Source ECG Toolbox, version 3.14, April 2019
% Released under the GNU General Public License
% Copyright (C) 2019  Reza Sameni
% Shiraz University, Shiraz, Iran
% reza.sameni@gmail.com

half_wlen = round(wlen/2);
T = length(x);

yy = filtfilt(ones(1, wlen), wlen, x.^2);
peaks0 = PeakDetection(yy, f/fs, 1);

I0 = find(peaks0);
I0_len = length(I0);
peaks = zeros(size(peaks0));
for i = 1 : I0_len
    index = max(I0(i) - half_wlen, 1) : min(I0(i) + half_wlen, T);
    [~, ii] = max(abs(x(index)));
    peaks(index(1) + ii - 1) = 1;
end
I = find(peaks);
D = diff(I);

% median approach:
if(strcmp(method,'median')) % median approach
    HR = 60*fs/median(D);
    HRPeriodSD = std(D);
elseif(strcmp(method, 'trmean')) % trimmed-mean approach
    D = sort(D);
    HR = 60*fs/mean(D(trim+1:end-trim));
    HRPeriodSD = std(D(trim+1:end-trim));
else
    error('unknown method');
end
