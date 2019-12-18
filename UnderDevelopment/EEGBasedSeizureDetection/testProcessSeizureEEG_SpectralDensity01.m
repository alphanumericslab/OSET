close all;
clear;
clc;

% parameters
subject1 = 'Patient_1';
subject2 = 'Patient_1';
number1 = 9;
number2 = 9;
path1 = ['E:\Sameni\Projects\Seizure\' subject1 '\'];
path2 = ['E:\Sameni\Projects\Seizure\' subject2 '\'];
mode1 = 'interictal';%'preictal';
mode2 = 'preictal';%'interictal';
% epochlen = 30; % epoch length in seconds
% epochoverlap = 20; % overlap between epochs in seconds
% sourcenum = 3;

nfft = 10000 - 1;
SpectralWinLen = 5; % in seconds
SpectralOverlapPercentage = 50; % [0 99]


% load signals
% [x fsx sequencex] = LoadSeizureEEG(path, subject, mode1, number1, ch, time);
% [y fsy sequencey] = LoadSeizureEEG(path, subject, mode2, number2, ch, time);
[x0 fsx] = LoadSeizureEEGFull(path1, subject1, mode1, number1);
[y0 fsy] = LoadSeizureEEGFull(path2, subject2, mode2, number2);

fs = 160;
x = resample(x0', fs, round(fsx))';
y = resample(y0', fs, round(fsy))';

% pre-process
x = x - LPFilter(x, 1/fs);
y = y - LPFilter(y, 1/fs);

% IIR notch
f0 = 60;
Q = 30;
Wo = f0/(fs/2);  BW = Wo/Q;
[b,a] = iirnotch(Wo,BW);
x = filter(b, a, x, [], 2);
y = filter(b, a, y, [], 2);


% estimate spectrum
h = spectrum.welch('Hamming', round(SpectralWinLen*fs), SpectralOverlapPercentage); % Create a Welch spectral estimator.
% h = spectrum.periodogram('Hamming'); % Create a Periodogram spectral estimator.
% h = spectrum.periodogram('Rectangular');
% h = spectrum.yulear;
Hx = MultiChannelSpectrum(x, fs, nfft, h); % Estimate the PSD
Hy = MultiChannelSpectrum(y, fs, nfft, h); % Estimate the PSD

% [y, W, A] = SCA(x, .5/fs, 7/fs);

figure
hold on
for i = 1:size(x, 1),
    %     plot(H(i).Frequencies, 10*log10(H(i).data/mean(H(1).data.^2)));
    plot(Hx(i).Frequencies, log(Hx(i).data/sqrt(sum(Hx(i).data.^2))),'b');
    plot(Hy(i).Frequencies, log(Hy(i).data/sqrt(sum(Hy(i).data.^2))),'r');
end
grid;


PlotECG(x, 4, 'b', fs);
PlotECG(y, 4, 'r', fs);

% PlotECG(y, 4, 'r', fs);
