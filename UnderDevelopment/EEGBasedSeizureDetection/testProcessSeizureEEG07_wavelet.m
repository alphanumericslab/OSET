close all;
clear;
clc;

% parameters
subject = 'Dog_2';
number = 19;
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
% x = LPFilter(x, 50.0/fsx); % lowpass filter
x = x - LPFilter(x, 1.0/fsx); % highpass filter
% % 
% y = LPFilter(y, 50.0/fsy); % lowpass filter
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

vx = zeros(size(x, 1), 1);
vy = zeros(size(y, 1), 1);

level = 3;
for i = 1:size(x, 1),
    [cx(i, :), lx(i, :)] = wavedec(x(i, :)/std(x(i, :)), level, 'db4');
    [cy(i, :), ly(i, :)] = wavedec(y(i, :)/std(y(i, :)), level, 'db4');
    
    indx = [0 cumsum(lx(i,1:end-1),2)];
    indy = [0 cumsum(ly(i,1:end-1),2)];
    
    for j = 1 : level+1,
        px(i, j) = mean(cx(i, indx(j)+1 : indx(j + 1)).^2);
        py(i, j) = mean(cy(i, indy(j)+1 : indy(j + 1)).^2);
    end
%     cx(i, :) = cumsum(cx(i, :).^2)/length(cx(i, :));
%     cy(i, :) = cumsum(cy(i, :).^2)/length(cy(i, :));
%     
%     vx(i) = sum((1:length(cx(i, :))).*cx(i, :)/length(cx(i, :)));
%     vy(i) = sum((1:length(cy(i, :))).*cy(i, :)/length(cy(i, :)));
end

% disp(vx-vy);

% disp(vy);

% estimate spectrum
% h = spectrum.welch('Hamming', round(SpectralWinLen*fsx), SpectralOverlapPercentage); % Create a Welch spectral estimator.
% h = spectrum.periodogram('Hamming'); % Create a Periodogram spectral estimator.
% h = spectrum.periodogram('Rectangular');
% h = spectrum.yulear;
% Hx = MultiChannelSpectrum(x, fsx, nfft, h); % Estimate the PSD 
% Hy = MultiChannelSpectrum(y, fsy, nfft, h); % Estimate the PSD 

% [Vx Dx] = eig(Cx);
% [Vy Dy] = eig(Cy);
% 
% Dx = diag(Dx);
% Dy = diag(Dy);
% 
% Dx = Dx(end:-1:1);
% Dy = Dy(end:-1:1);
% 
% Dx = Dx/Dx(1);
% Dy = Dy/Dy(1);
% 
% figure
% hold on
% bar(Dx,'b');
% bar(Dy,'r');
% h = findobj(gca,'Type','patch');
% set(h(2), 'FaceColor', 'r', 'FaceAlpha', 0.5, 'EdgeColor', 'r', 'EdgeAlpha', 0.5);
% set(h(1), 'FaceColor', 'b', 'FaceAlpha', 0.5, 'EdgeColor', 'b', 'EdgeAlpha', 0.5);
% grid

% figure
% hold on
% for i = 1:size(x, 1),
%     plot(Hx(i).Frequencies, Hx(i).data/sqrt(mean(Hx(i).data.^2)),'b');
%     plot(Hy(i).Frequencies, Hy(i).data/sqrt(mean(Hy(i).data.^2)),'r');
% end
% grid;
% 
% PlotECG(x, 4, 'b', fsx);
% PlotECG(y, 4, 'r', fsy);

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

figure
hold on
plot(cumsum(px,2)'/(1),'b');
plot(cumsum(py,2)'/(1),'r');
plot(cumsum(px,2)'/(1),'bo');
plot(cumsum(py,2)'/(1),'ro');
grid

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