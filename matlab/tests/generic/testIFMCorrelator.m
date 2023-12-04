clc;
close all;
clear

T = 100000;
fs = 500;
f0 = 30;
BW = 5;
f = 30.5;
n = (0:T-1);
x = sin(2*pi*f/fs*n) + 0.5*randn(1, T);
tau = 3;

ff = IFMCorrelator(x, f0/fs, BW/fs, tau);

figure;
plot(f0 - ff*fs);
grid