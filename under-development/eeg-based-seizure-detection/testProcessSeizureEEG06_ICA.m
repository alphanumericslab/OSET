close all;
clear;
clc;

% parameters
subject = 'Dog_3';
number = 15;
path = ['J:\Seizure\' subject '\'];
mode1 = 'interictal';
mode2 = 'preictal';
ch = 1:16; % desired channels
time = [7 10]*60; % in seconds
nfft = 5000 - 1;
SpectralWinLen = 10; % in seconds
SpectralOverlapPercentage = 50; % [0 99]

% load signals
[x fsx sequencex] = LoadSeizureEEG(path, subject, mode1, number, ch, time);
[y fsy sequencey] = LoadSeizureEEG(path, subject, mode2, number, ch, time);

% pre-process
x = LPFilter(x, 50.0/fsx); % lowpass filter
x = x - LPFilter(x, 1.0/fsx); % highpass filter
% 
y = LPFilter(y, 50.0/fsy); % lowpass filter
y = y - LPFilter(y, 1.0/fsy); % highpass filter

% xxLP = LPFilter(x, 5.0/fsx);
% xxBP = BPFilter(x, 7.0/fsx, 15.0/fsx);
% yyLP = LPFilter(y, 5.0/fsy);
% yyBP = BPFilter(y, 7.0/fsy, 15.0/fsy);
% 
% rx = sum(xxLP.^2, 2)./sum(xxBP.^2, 2);
% ry = sum(yyLP.^2, 2)./sum(yyBP.^2, 2);
% 
% medpowerradio_x = median(rx)
% medpowerradio_y = median(ry)

% wx = jadeR(x);
% wy = jadeR(y);
% sx = wx*x;
% sy = wy*y;

[wxinv, sx] = sobi(x);
[wyinv, sy] = sobi(y);

% fl = 0.5;
% fu = 5.0;
% sx = SCA2(x,fl/fsx,fu/fsx);
% sy = SCA2(y,fl/fsy,fu/fsy);

% % estimate spectrum
% h = spectrum.welch('Hamming', round(SpectralWinLen*fsx), SpectralOverlapPercentage); % Create a Welch spectral estimator.
% % h = spectrum.periodogram('Hamming'); % Create a Periodogram spectral estimator.
% % h = spectrum.periodogram('Rectangular');
% % h = spectrum.yulear;
% Hx = MultiChannelSpectrum(sx, fsx, nfft, h); % Estimate the PSD 
% Hy = MultiChannelSpectrum(sy, fsy, nfft, h); % Estimate the PSD 

Cx = cov(x');
Cy = cov(y');

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
bar(Dx,'b');
bar(Dy,'r');
h = findobj(gca,'Type','patch');
set(h(2), 'FaceColor', 'r', 'FaceAlpha', 0.5, 'EdgeColor', 'r', 'EdgeAlpha', 0.5);
set(h(1), 'FaceColor', 'b', 'FaceAlpha', 0.5, 'EdgeColor', 'b', 'EdgeAlpha', 0.5);
grid

% figure
% hold on
% for i = 1:size(x, 1),
%     plot(Hx(i).Frequencies, Hx(i).data/sqrt(mean(Hx(i).data.^2)),'b');
%     plot(Hy(i).Frequencies, Hy(i).data/sqrt(mean(Hy(i).data.^2)),'r');
% end
% grid;

PlotECG(x, 4, 'k', fsx);
PlotECG(y, 4, 'm', fsy);

PlotECG(sx, 4, 'b', fsx);
PlotECG(sy, 4, 'r', fsy);

% for i = 1:size(x, 1),
%     figure;
%     hold on;
%     plot(x(i, :));
%     plot(y(i, :),'r');
%     grid
%     title(num2str(i));
% end

% for i = 1:size(x, 1),
%     figure;
%     hold on;
%     plot(cy(i, :),'r');
%     plot(cx(i, :));
%     grid
%     title(num2str(i));
% end

% c = 1;
% nwindow = round(fsx*3.0);
% noverlap = round(fsx*2.5);
% nfft = 1000;
% 
% figure
% subplot(211);
% spectrogram(x(c, :), nwindow, noverlap, nfft, fsx, 'yaxis');
% subplot(212);
% spectrogram(y(c, :), nwindow, noverlap, nfft, fsy, 'yaxis');