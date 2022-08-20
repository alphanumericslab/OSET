function [HR peaks HRPeriodSD] = HRCalculation3(x, f, fs, trim, wlen, method)
%
% [HR peaks HRPeriodSD] = HRCalculation3(x,f,fs,trim,wlen,method),
% Multichannel Heartrate calculator. The R-peak detector is based on max search 
% over the energy envelope of multiple channels and the HR is the median 
% of the HR over the data window.
%
% inputs:
% x: matrix of input data (N channels by T samples)
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
% Open Source Electrophysiological Toolbox, version 3.14, February 2019
% Released under the GNU General Public License
% Copyright (C) 2019  Reza Sameni
% Shiraz University, Shiraz, Iran
% reza.sameni@gmail.com

y = zeros(size(x));
for k = 1 : size(y, 1),
    y(k, :) = filtfilt(ones(1, wlen), wlen, x(k, :).^2);
end
yy = sqrt(sum(y, 1));
peaks = PeakDetection(yy, f/fs, 1);
% peaks = PeakDetection5(x,f/fs);
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
