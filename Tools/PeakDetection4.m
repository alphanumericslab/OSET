function [peaks,r] = PeakDetection4(ref,fs,h,fmax)
%
% peaks = PeakDetection4(x,fs,h,fmax),
% R-peak detector based on a matched filter
%
% inputs:
% x: vector of input data
% fs: sampling rate
% h: template waveform
% fmax: maximum expected frequency of the R-peaks
%
% output:
% peaks: vector of R-peak impulse train
%
%
% Open Source ECG Toolbox, version 2.0, March 2008
% Released under the GNU General Public License
% Copyright (C) 2008  Reza Sameni
% Sharif University of Technology, Tehran, Iran -- GIPSA-Lab, INPG, Grenoble, France
% reza.sameni@gmail.com

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.

N = length(ref);
L = length(h);

h = h(end:-1:1);

w = floor(L/2);

r = filter(h,1,[ref zeros(1,w-1)]);

r = r(w:N+w-1);

peaks = PeakDetection(r , fmax/fs, 1);