%
% Test program for R-peak detection of ECG signals.
%
% Dependencies: The baseline wander and ECG filtering toolboxes of the Open Source ECG Toolbox
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

clc
clear all
close all;

load('SampleECG1.mat');
% load('SampleECG2.mat'); data = data(1:15000,5)';

fs = 1000;
t = (0:length(data)-1)/fs;

f = 1;                                      % approximate R-peak frequency

bsline = LPFilter(data,.7/fs);              % baseline wander removal
x = data-bsline;

peaks1 = PeakDetection(x,f/fs);     % peak detection (max detection)
peaks2 = PeakDetection2(x,fs);      % peak detection (Pan-Tompkins)
[Y I] = max(abs(x));
peaks3 = PeakDetection3(x,fs,x(I-50:I+49),.2,1.5);    % peak detection (matched filter)

figure;
plot(t,x,'b');
hold on
plot(t,peaks1.*x,'ro');
plot(t,peaks2.*x,'go');
plot(t,peaks3.*x,'mo');
xlabel('time (sec.)');
legend('ECG','ECG Peaks(max detection)','ECG Peaks(Pan-Tompkins)','ECG Peaks(matched filter)');
grid;
