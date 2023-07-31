%
% Test program for RelativeSpectralPower
%
% Dependencies: RelativeSpectralPower() from the Open Source ECG Toolbox
%
%
% Reza Sameni, June 2022
% reza.sameni@gmail.com
% Open-Source Electrophysiological Toolbox (OSET)
% https://github.com/alphanumericslab/OSET.git

clc
clear
close all;

load('SampleECG2.mat'); data = data(1:15000, :)';

fs = 1000;
f0 = 50;
data = data - LPFilter(data, 1/fs);

n = 0 : length(data)-1;
harmonic1 = 0.09 * sin(2*pi*n*f0/fs + pi/9)./(std(data, [], 2)*ones(1, size(data, 2)));
harmonic2 = 0.09 * sin(2*pi*n*3*f0/fs + pi/6)./(std(data, [], 2)*ones(1, size(data, 2)));

x = data + harmonic1 + harmonic2;

freqs = [f0, 2*f0, 3*f0];
Q_factor = 35;
harmonics_ratios_percentage = 100 * RelativeSpectralPower(x, fs, freqs, Q_factor)
