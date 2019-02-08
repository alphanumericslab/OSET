function [HR peaks] = HRCalculation1(x,f,fs,trim,varargin)
%
% HR = HRCalculation1(x,f,flag),
% Heartrate calculator. The R-peak detector is based on max search and the
% HR is the median of the HR over the data window.
%
% inputs:
% x: vector of input data
% f: approximate ECG beat-rate in Hertz, normalized by the sampling frequency
% flag: search for positive (flag=1) or negative (flag=0) peaks. By default
% the maximum absolute value of the signal, determines the peak sign.
%
% output:
% HR: the median of the HR over the data window.
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

peaks = PeakDetection5(x,f/fs);
I = find(peaks);
D = diff(I);

% median approach:
% HR = 60*fs/median(D);

% trimmed-mean approach
D = sort(D);
HR = 60*fs/mean(D(trim+1:end-trim));
