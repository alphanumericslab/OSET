close all;
clear;
clc;

% parameters
subject = 'Dog_1';
number = 2;
path = ['E:\Sameni\Projects\Seizure\' subject '\'];
mode1 = 'interictal';
mode2 = 'preictal';
ch = 1:16; % desired channels
time = [8 10]*60; % in seconds
nfft = 5000 - 1;
SpectralWinLen = 10; % in seconds
SpectralOverlapPercentage = 50; % [0 99]

% load signals
[x fsx sequencex] = LoadSeizureEEG(path, subject, mode1, number, ch, time);
[y fsy sequencey] = LoadSeizureEEG(path, subject, mode2, number, ch, time);


devfraction = 3.0;
fc = 10.0;
bw = 4;

x = LPFilter(x, (fc + bw/2)/fsx); % lowpass filter
x = x - LPFilter(x, (fc - bw/2)/fsx); % highpass filter
y = LPFilter(y, (fc + bw/2)/fsy); % lowpass filter
y = y - LPFilter(y, (fc - bw/2)/fsy); % highpass filter

xPM = zeros(size(x));
yPM = zeros(size(y));

% x(1,:) = fmmod(sin(2*pi*(1:length(x))*fc/fsx), fc, fsx, bw/2);
% x(1,:) = pmmod(sin(2*pi*(1:length(x))*fc/fsx), fc, fsx, devfraction/2);


for i = 1:size(x,1),
    xPM(i, :) = pmdemod(x(i, :), fc, fsx, pi/devfraction, 0.0);
    yPM(i, :) = pmdemod(y(i, :), fc, fsy, pi/devfraction, 0.0);
%     xPM(i, :) = fmdemod(x(i, :), fc, fsx, bw/2, fc);
%     yPM(i, :) = fmdemod(y(i, :), fc, fsy, bw/2, fc);
end

% xPM = xPM - mean(xPM, 2)*ones(1, size(xPM, 2));
% yPM = yPM - mean(yPM, 2)*ones(1, size(yPM, 2));

xPM = unwrap(xPM, devfraction, 2);
yPM = unwrap(yPM, devfraction, 2);
% xPM = cumsum(unwrap(xPM, bw/2, 2)/fsx, 2);
% yPM = cumsum(unwrap(yPM, bw/2, 2)/fsy, 2);

% xPM = xPM - LPFilter(xPM, 0.1/fsx);
% yPM = yPM - LPFilter(yPM, 0.1/fsy);

% % Cx = cov(x');
% % Cy = cov(y');
Cx = cov(xPM');
Cy = cov(yPM');

% estimate spectrum
h = spectrum.welch('Hamming', round(SpectralWinLen*fsx), SpectralOverlapPercentage); % Create a Welch spectral estimator.
% h = spectrum.periodogram('Hamming'); % Create a Periodogram spectral estimator.
% h = spectrum.periodogram('Rectangular');
% h = spectrum.yulear;
Hx = MultiChannelSpectrum(x, fsx, nfft, h); % Estimate the PSD 
Hy = MultiChannelSpectrum(y, fsy, nfft, h); % Estimate the PSD 

% [y, W, A] = SCA(x, .5/fs, 7/fs);


[Vx Dx] = eig(Cx);
[Vy Dy] = eig(Cy);

Dx = diag(Dx);
Dy = diag(Dy);

Dx = Dx(end:-1:1);
Dy = Dy(end:-1:1);

Dx = Dx/Dx(1);
Dy = Dy/Dy(1);

figure
hold on
bar(log(Dx),'b');
bar(log(Dy),'r');
h = findobj(gca,'Type','patch');
set(h(2), 'FaceColor', 'r', 'FaceAlpha', 0.5, 'EdgeColor', 'r', 'EdgeAlpha', 0.5);
set(h(1), 'FaceColor', 'b', 'FaceAlpha', 0.5, 'EdgeColor', 'b', 'EdgeAlpha', 0.5);
grid

figure
hold on
for i = 1:size(x, 1),
    plot(Hx(i).Frequencies, Hx(i).data/sqrt(mean(Hx(i).data.^2)),'b');
    plot(Hy(i).Frequencies, Hy(i).data/sqrt(mean(Hy(i).data.^2)),'r');
end
grid;

PlotECG(xPM, 4, 'b', fsx);
PlotECG(yPM, 4, 'r', fsy);

% PlotECG(x, 4, 'b', fsx);
% PlotECG(y, 4, 'r', fsy);

c = 2;
figure;
hold on;
plot(x(c, :));
plot(xPM(c, :), 'r');
grid

