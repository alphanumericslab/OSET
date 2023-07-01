function [HR peaks HRPeriodSD] = HRCalculation2(x, f, fs, trim, wlen, method)
%
% [HR peaks HRPeriodSD] = HRCalculation2(x,f,fs,trim,wlen,method),
% Heartrate calculator. The R-peak detector is based on max search and the
% HR is the median of the HR over the data window.
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
% Open Source ECG Toolbox, version 2.0, November 2009
% Released under the GNU General Public License
% Copyright (C) 2009  Reza Sameni
% Shiraz University, Shiraz, Iran
% reza.sameni@gmail.com

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.

% peaks = PeakDetection5(x,f/fs);
yy = filtfilt(ones(1, wlen), wlen, x.^2);
peaks = PeakDetection(yy, f/fs, 1);

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
