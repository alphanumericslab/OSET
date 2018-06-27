%
% Test program for generating single-channel ECGs with additive colored
% noise
%
% Dependencies: The synthetic ECG generator and noise generator package of
%   the Open Source ECG Toolbox
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

%//////////////////////////////////////////////////////////////////////////
clc
clear
close all;

load('SampleECG1.mat');
% load('SampleECG2.mat'); data = data(1:15000,2)';

fs = 1000;
t = (0:length(data)-1)/fs;

f = 1;                                          % approximate R-peak frequency

bsline = LPFilter(data,.7/fs);                  % baseline wander removal (may be replaced by other approaches)
%bsline = BaseLineKF(data,.5/fs);                % baseline wander removal (may be replaced by other approaches)

x = data - bsline;

peaks = PeakDetection(x,f/fs);                  % peak detection

[phase, phasepos] = PhaseCalculation(peaks);     % phase calculation

teta = 0;                                       % desired phase shift
pphase = PhaseShifting(phase,teta);             % phase shifting

bins = 250;                                     % number of phase bins
[ECGmean,ECGsd,meanphase] = MeanECGExtraction(x,pphase,bins,1); % mean ECG extraction 

ECGBeatFitter(ECGmean,ECGsd,meanphase);           % ECG beat fitter GUI

L = length(OptimumParams)/3;% number of Gaussian kernels
alphai = OptimumParams(1:L);
bi = OptimumParams(L+1:2*L);
tetai = OptimumParams(2*L+1:3*L);
% teta0 = pi/2;
teta0 = 0;
[ECG, teta]= SingleChannelECGGenerator(pphase',teta0,alphai,bi,tetai);

%//////////////////////////////////////////////////////////////////////////
% data plotting
figure;
subplot(311);
plot(t,1000*data);
grid
xlabel('time(s)');
ylabel('Amplitude(mV)');
title('Original ECG');
subplot(312);
plot(t,1000*x,'r');
grid
xlabel('time(s)');
ylabel('Amplitude(mV)');
title('Original ECG');
subplot(313);
plot(t,1000*ECG,'m');
grid
xlabel('time(s)');
ylabel('Amplitude(mV)');
title('Synthetic ECG');
