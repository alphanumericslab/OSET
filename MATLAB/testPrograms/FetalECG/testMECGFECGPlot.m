clear
clc
close all;

fs = 1000;
load('SampleMECG-FECG1.mat');
plot(((6001:12000)-6000)/fs,testdata(6001:12000))
grid
xlabel('time(s)');
ylabel('Amplitude(mV)');
title('Sample maternal-fetal ECG segment (after incomplete mECG cancellation)');