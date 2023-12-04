close all;
clear;
clc;

% parameters
path = 'E:\Sameni\Projects\Seizure\Dog_1\';
subject = 'Dog_1';
number = 4;
% mode = 'interictal';
mode = 'preictal';
ch = 1:16; % desired channels
time = [5 10]*60; % in seconds
nfft = 10000 - 1;
SpectralWinLen = 10; % in seconds
SpectralOverlapPercentage = 50; % [0 99]

% load signals
[x0 fs sequence] = LoadSeizureEEG(path, subject, mode, number, ch, time);

% fs = 100;
% x = resample(x0', fs, fs0)';

% pre-process
x = LPFilter(x0, 10/fs);
% x = x0 - LPFilter(x0, 0.1/fs);
% x = x0;

% estimate spectrum
h = spectrum.welch('Hamming', round(SpectralWinLen*fs), SpectralOverlapPercentage); % Create a Welch spectral estimator.
% h = spectrum.periodogram('Hamming'); % Create a Periodogram spectral estimator.
% h = spectrum.periodogram('Rectangular');
% h = spectrum.yulear;
H = MultiChannelSpectrum(x, fs, nfft, h); % Estimate the PSD 

% [y, W, A] = SCA(x, .5/fs, 7/fs);

figure
hold on
for i = 1:size(x, 1),
%     plot(H(i).Frequencies, 10*log10(H(i).data/mean(H(1).data.^2)));
    plot(H(i).Frequencies, H(i).data/mean(H(i).data.^2));
end
grid;


PlotECG(x, 4, 'b', fs);

% PlotECG(y, 4, 'r', fs);

