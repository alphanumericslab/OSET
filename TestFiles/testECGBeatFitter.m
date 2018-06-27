%
% Test program for mean ECG beat extraction and parameter optimization.
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
load('SampleECG2.mat'); data = data(1:15000,6)';
fs = 1000;

t = (0:length(data)-1)/fs;

f = 1;                                          % approximate R-peak frequency

bsline = LPFilter(data,.7/fs);                  % baseline wander removal (may be replaced by other approaches)
%bsline = BaseLineKF(data,.5/fs);                % baseline wander removal (may be replaced by other approaches)

data1 = data-bsline;

%//////////////////////////////////////////////////////////////////////////
% Making the data noisy
SNR = 0;
SignalPower = mean(data1.^2);
NoisePower = SignalPower / 10^(SNR/10);
x = data1 + sqrt(NoisePower)*randn(size(data1));
%//////////////////////////////////////////////////////////////////////////

peaks = PeakDetection(x,f/fs);                  % peak detection

[phase phasepos] = PhaseCalculation(peaks);     % phase calculation

teta = 0;                                       % desired phase shift
pphase = PhaseShifting(phase,teta);             % phase shifting

bins = round(fs/3);                                     % number of phase bins
[ECGmean,ECGsd,meanphase] = MeanECGExtraction(x,pphase,bins,1); % mean ECG extraction 

ECGBeatFitter(ECGmean,ECGsd,meanphase,'OptimalParams');               % ECG beat fitter GUI

% display the optimal parameters
L = length(OptimalParams)/3;
ai = OptimalParams(1:L)
bi = OptimalParams(L+1:2*L)
tetai = OptimalParams(2*L+1:3*L)


figure;
plot(t,data*6/max(data),'b');
hold on
plot(t,peaks*2,'ro');
plot(t,phase,'c--','linewidth',1);
plot(t,pphase,'g','linewidth',1);
grid;
xlabel('time (sec.)');
legend('Scaled ECG','ECG Peaks','phase','shifted phase');

figure;
errorbar(meanphase,ECGmean,ECGsd/2);
hold on;
plot(meanphase,ECGmean,'r');
legend('SD Bar','Mean ECG');
xlabel('phase (rads.)');
grid

