% Test program for fetal ECG ploting
%
% Dependencies: The baseline wander and ECG filtering toolboxes of the Open Source ECG Toolbox
%
% Sample code for testing naive maternal template subtraction technique
%
% Reza Sameni (C)
% Email: rsameni@shirazu.ac.ir
% Web: www.sameni.info
%
% Crated 2006
% Modified June 2018


clc;
clear;
close all;

fs = 1600;

load CC20060830_ch1_fECG fecg;
% load fECG1_CC_100sec fecg

N = length(fecg);
t = (0:N-1)/fs;

% fecg2 = fecg - LPFilter(Median(fecg,N,50,100)',50/fs);
fecg2 = fecg - LPFilter(fecg,50/fs);
fecg2 = LPFilter(fecg2,150/fs);

h = fecg2(110:150);
% % % h = fecg2(31090:31145);
peaks1 = PeakDetection3(fecg2,fs,h,.001,2.5);

peaks2 = PeakDetection(fecg2,2/fs,1);

I1 = find(peaks1);
I2 = find(peaks2);

figure;
hold on;
plot(t,fecg2);
plot(t(I1),fecg2(I1),'ro');
plot(t(I2),fecg2(I2),'go');
grid
legend('fECG', 'Matched filter peak detector', 'Simple peak detector');

figure;
hold on
plot(t(I1(1:end-1)),fs*60./diff(I1),'b');
plot(t(I2(1:end-1)),fs*60./diff(I2),'r');
grid;
xlabel('time(s)','FontSize',16);
ylabel('Heart Rate (BPM)','FontSize',16);
set(gca,'Box','On','FontSize',16);
legend('Matched filter peak detector', 'Simple peak detector');

