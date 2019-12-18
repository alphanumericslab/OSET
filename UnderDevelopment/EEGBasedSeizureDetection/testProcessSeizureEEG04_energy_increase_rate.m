close all;
clear;
clc;

% parameters
subject = 'Dog_3';
number = 10;
path = ['E:\Sameni\Projects\Seizure\' subject '\'];
mode1 = 'interictal';
mode2 = 'preictal';
ch = 1:16; % desired channels
time = [1 10]*60; % in seconds
nfft = 10000;
SpectralWinLen = 10; % in seconds
SpectralOverlapPercentage = 50; % [0 99]

% load signals
% [x fsx sequencex] = LoadSeizureEEG(path, subject, mode1, number, ch, time);
% [y fsy sequencey] = LoadSeizureEEG(path, subject, mode2, number, ch, time);
[x fsx] = LoadSeizureEEGFull(path, subject, mode1, number);
[y fsy] = LoadSeizureEEGFull(path, subject, mode2, number);

% pre-process
% x = LPFilter(x, 150.0/fsx); % lowpass filter
% x = x - LPFilter(x, 0.1/fsx); % highpass filter
%
% y = LPFilter(y, 150.0/fsy); % lowpass filter
% y = y - LPFilter(y, 0.1/fsy); % highpass filter

% resample data
ffs = 150;
if(fsx > ffs) % for the patients who have high sampling rates
    x = resample(x', ffs, round(fsx))';
    fsx = ffs;
end
if(fsy > ffs) % for the patients who have high sampling rates
    y = resample(y', ffs, round(fsy))';
    fsy = ffs;
end


f0 = 60; % powerline frequency Hz
Q = 30; % notch filter Q-factor
% IIR notch
Wo = f0/(fsx/2);  BW = Wo/Q;
[b,a] = iirnotch(Wo,BW);
x = filter(b, a, x, [], 2);
y = filter(b, a, y, [], 2);


X = fft(x, nfft, 2);
Y = fft(y, nfft, 2);

X = X./(sqrt(sum(abs(X(:,2:end/2)).^2, 2))*ones(1,size(X,2)));
Y = Y./(sqrt(sum(abs(Y(:,2:end/2)).^2, 2))*ones(1,size(Y,2)));

CumSumX = cumsum(abs(X(:,2:end/2)).^2, 2);
CumSumY = cumsum(abs(Y(:,2:end/2)).^2, 2);

meanX = mean(CumSumX, 1);
meanY = mean(CumSumY, 1);

dfX = diff(meanX);
dfY = diff(meanY);

csX = mean(meanX)
csY = mean(meanY)

fx = fsx*(0:nfft-1)/nfft;
fy = fsy*(0:nfft-1)/nfft;

figure;
hold on;
plot(fx, abs(X)', 'b');
plot(fy, abs(Y)', 'r');
grid

figure;
hold on;
plot(CumSumX', 'b');
plot(meanX, 'b', 'linewidth', 3);
plot(CumSumY', 'r');
plot(meanY, 'r', 'linewidth', 3);
grid

figure
hold on;
plot(dfX)
plot(dfY, 'r')
grid
% PlotECG(x, 4, 'b', fsx);
% PlotECG(y, 4, 'r', fsy);

