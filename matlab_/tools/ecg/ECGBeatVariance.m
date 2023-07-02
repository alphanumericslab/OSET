function [ECGmean, ECGsd] = ECGBeatVariance(x, peaks, wlen)
%
% [ECGmean, ECGsd] = ECGBeatVariance(x, peaks, wlen)
% Calculation of the mean and SD of ECG waveforms in different beats
%
% inputs:
% x: input ECG signal
% peaks: ECG peaks
% wlen: window length around the peaks (from each side of the peak)
% outputs:
% ECGmean: mean ECG beat
% ECGsd: standard deviation of ECG beats
%
% Open Source ECG Toolbox, version 2.0, March 2010
% Released under the GNU General Public License
% Copyright (C) 2010  Reza Sameni
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

x = [zeros(1,wlen) x zeros(1,wlen)];
peaks = [zeros(1,wlen) peaks zeros(1,wlen)];
I = find(peaks);

N = length(I);

ECG = zeros(N,2*wlen);
for i = 1 : N,
    ECG(i,:) = x(I(i)-wlen+1:I(i)+wlen);
end

ECGmean = mean(ECG,1);
ECGsd = std(ECG,1,1);


