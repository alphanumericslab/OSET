%
% Test program for mean ECG extraction.
%
% By:
% Reza Sameni, September 2006
% LIS-INPG, Grenoble, France - Sharif University of Technology, Tehran, Iran
% reza.sameni@gmail.com

clc
clear
close all;

load('SampleECG2.mat'); data = data(1:15000,2)';

fs = 1000;
t = (0 : length(data)-1)/fs;

f = 1;                                          % approximate R-peak frequency

bsline = LPFilter(data,.7/fs);                  % baseline wander removal (may be replaced by other approaches)
x = data-bsline;
peaks = PeakDetection(x,f/fs);                  % peak detection

[phase, phasepos] = PhaseCalculation(peaks);     % phase calculation

teta = 0;                                       % desired phase shift
pphase = PhaseShifting(phase,teta);             % phase shifting

bins = 500;                                     % number of phase bins
[ECGmean,ECGsd,meanphase] = MeanECGExtraction(x,pphase,bins,1);

[prms,mdl,error,approach] = ECGBeatFitterAuto(ECGmean,meanphase);

plot(meanphase,ECGmean);
hold on;
plot(meanphase,mdl,'r');
grid; 

