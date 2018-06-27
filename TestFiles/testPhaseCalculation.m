%
% Test program for phase calculation of ECG signals.
%
% Dependencies: The baseline wander and ECG filtering toolboxes of the Open Source ECG Toolbox
%
% Open Source ECG Toolbox, version 1.0, November 2006
% Released under the GNU General Public License
% Copyright (C) 2006  Reza Sameni
% Sharif University of Technology, Tehran, Iran -- LIS-INPG, Grenoble, France
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

% load('SampleECG1.mat');
load('SampleECG2.mat'); data = data(1:end,15)';

fs = 1000;
t = (0:length(data)-1)/fs;

f = 1;                                      % approximate R-peak frequency

bsline = LPFilter(data,.7/fs);              % baseline wander removal (may be replaced by other approaches)
x = data-bsline;
peaks = PeakDetection(x,f/fs);              % peak detection

[phase phasepos] = PhaseCalculation(peaks); % phase calculation

[X,Y,Z] = pol2cart(phase,1,x);

I = find(peaks);
figure;
% plot(t,data*2*pi/max(data),'b');
plot(t,2*pi*x/max(x),'b');
hold on
plot(t(I),peaks(I)*2,'ro');
plot(t,phasepos,'c--');
plot(t,phase,'g','linewidth',2);
grid;
xlabel('time (sec.)');
legend('Scaled ECG','ECG Peaks','\phi_{positive}','\phi');

figure;
plot3(X,Y,Z);
grid;
title('phase-wrapped ECG');
xlabel('X');
ylabel('Y');
zlabel('ECG');