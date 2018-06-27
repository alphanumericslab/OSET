%
% Test program for automatic ECG model fitting
%
% By:
% Reza Sameni, September 2010
% Shiraz University, Shiraz, Iran
% reza.sameni@gmail.com

clc
clear all
close all;

load('SampleECG2.mat'); data = data(1:15000,2)';

fs = 1000;
t = [0:length(data)-1]/fs;

f = 1;                                          % approximate R-peak frequency

bsline = LPFilter(data,.7/fs);                  % baseline wander removal (may be replaced by other approaches)
x = data - bsline +.5*randn(size(data));
peaks = PeakDetection(x,f/fs);                  % peak detection

[phase phasepos] = PhaseCalculation(peaks);     % phase calculation

teta = 0;                                       % desired phase shift
pphase = PhaseShifting(phase,teta);             % phase shifting

bins = fs/2;                                     % number of phase bins
[ECGmean,ECGsd,meanphase] = MeanECGExtraction(x,pphase,bins,1);

stopth = .1; % percentage of error
num = 10;
maxitr = 10;
approach = 'min error'; % fixed number/min error
energyth = 0.001;
wlen = 5;
[params model er indexes method] = ECGBeatFitterAuto(ECGmean,meanphase,approach,stopth,num,maxitr,energyth,wlen);
method

figure;
hold on;
plot(meanphase,ECGmean);
plot(meanphase,model,'r');
grid;
xlabel('phase(rad.)');
legend('real beat','fitted beat');

