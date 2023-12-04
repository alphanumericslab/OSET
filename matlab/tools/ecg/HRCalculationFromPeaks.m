function HR = HRCalculationFromPeaks(peaks, fs, trim, method)
%
% HR = HRCalculationFromPeaks(peaks, fs, trim, method),
% Heartrate calculator from the R-peaks impulse vector
%
% inputs:
% peaks: an impulse vector containing the R-peaks
% fs: sampling frequency in Hz
% trim: number of beats to trim from averaging
% method: 'trmean' or 'median'
%
% output:
% HR: the trimmed mean/median of the HR over the data window.
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

I = find(peaks);
D = diff(I);

% median approach:
if(strcmp(method,'median')) % median approach
    HR = 60*fs/median(D);
elseif(strcmp(method, 'trmean')) % trimmed-mean approach
    D = sort(D);
    HR = 60*fs/mean(D(trim+1:end-trim));
else
    error('unknown method');
end
