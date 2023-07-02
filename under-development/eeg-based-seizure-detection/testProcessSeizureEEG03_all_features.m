close all;
clear;
clc;

% parameters
subject = 'Dog_2';
path = ['J:\Seizure\' subject '\'];
ch = 1:16; % desired channels
time = [7 10]*60; % in seconds

d = dir([path '*.mat']);
NumRecords = length(d);
mode = zeros(1, NumRecords);
powerradio = zeros(1, NumRecords);
A = zeros(1, NumRecords);
skew = zeros(1, NumRecords);
kurt = zeros(1, NumRecords);
f0 = zeros(1, NumRecords);
bw = zeros(1, NumRecords);
fskew = zeros(1, NumRecords);
for i = 1 : NumRecords,
    load([path d(i).name]);
    reg = regexp(d(i).name, '_');
    varname = [d(i).name(reg(2)+1:reg(4)-1) '_' num2str(str2double(d(i).name(end-7:end-4)))];
    md = d(i).name(reg(2)+1:reg(3)-1);
    if(isequal(md, 'interictal'))
        mode(i) = 1;
    elseif(isequal(md, 'preictal'))
        mode(i) = 2;
    elseif(isequal(md, 'test'))
        mode(i) = 3;
    end
    
    eval(['data = ' varname '; clear ' varname ';']);
    
    disp([num2str(i) ': ' md]);
    
    x = data.data;
    fs = data.sampling_frequency;
    %     sequence = data.sequence;
    samples = round(fs*time);
    samples(samples < 1) = 1;
    samples(samples > size(x,2)) = size(x,2);
    x = x(ch, samples(1):samples(2));
    
    
    % pre-process
    x = LPFilter(x, 40/fs); % lowpass filter
    x = x - LPFilter(x, 0.1/fs); % highpass filter
    
    % power ratio in bands feature
    xxLP = LPFilter(x, 5.0/fs);
    %     xxBP = BPFilter(x, 7.0/fs, 15.0/fs);
    %     rx = sum(xxLP.^2, 2)./sum(xxBP.^2, 2);
    rx = sum(xxLP.^2, 2)./sum(x.^2, 2);
    powerradio(i) = median(rx);
    
%     % frequency features
%     NFFT = 1000;
%     halfNFFT = round(NFFT/2);
%     
%     X = fft(x,NFFT);
%     XX = abs(X(1:halfNFFT)).^2;
%     f = (fs/2)*(0:halfNFFT-1)/halfNFFT;
%     A(i) = sqrt(mean(abs(x).^2));
%     skew(i) = skewness(x);
%     kurt(i) = kurtosis(x);
%     f0(i) = mean(f .* XX) / mean(XX);
%     bw(i) = sqrt(mean( (f - f0(i)).^2 .* XX ) / mean(XX));
%     fskew(i) = nthroot(mean( (f - f0(i)).^3 .* XX ) / mean(XX) , 3);
end

bins = 0.1:.01:.8;
figure
hist(powerradio(mode == 1), bins);
hold on
hist(powerradio(mode == 2), bins);
hist(powerradio(mode == 3), bins);
h = findobj(gca,'Type','patch');
set(h(2), 'FaceColor', 'r', 'FaceAlpha', 0.5, 'EdgeColor', 'r', 'EdgeAlpha', 0.5);
set(h(1), 'FaceColor', 'b', 'FaceAlpha', 0.5, 'EdgeColor', 'b', 'EdgeAlpha', 0.5);
set(h(3), 'FaceColor', 'g', 'FaceAlpha', 0.5, 'EdgeColor', 'g', 'EdgeAlpha', 0.5);
title('powerradio');
grid

% hist(powerradio(mode == 3),0:.1:5);
% h = findobj(gca,'Type','patch');
% set(h,'FaceColor','g','EdgeColor','w');

%
%
% mode1 = 'interictal';
% mode2 = 'preictal';
% nfft = 1000 - 1;
% SpectralWinLen = 10; % in seconds
% SpectralOverlapPercentage = 50; % [0 99]
%
% % load signals
% [x fsx sequencex] = LoadSeizureEEG(path, subject, mode1, number, ch, time);
% [y fsy sequencey] = LoadSeizureEEG(path, subject, mode2, number, ch, time);
%
% % pre-process
% x = LPFilter(x, 40/fsx); % lowpass filter
% x = x - LPFilter(x, 0.1/fsx); % highpass filter
%
% y = LPFilter(y, 40/fsy); % lowpass filter
% y = y - LPFilter(y, 0.1/fsy); % highpass filter
%
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
%
% % estimate spectrum
% h = spectrum.welch('Hamming', round(SpectralWinLen*fsx), SpectralOverlapPercentage); % Create a Welch spectral estimator.
% % h = spectrum.periodogram('Hamming'); % Create a Periodogram spectral estimator.
% % h = spectrum.periodogram('Rectangular');
% % h = spectrum.yulear;
% Hx = MultiChannelSpectrum(x, fsx, nfft, h); % Estimate the PSD
% Hy = MultiChannelSpectrum(y, fsy, nfft, h); % Estimate the PSD
%
% % [y, W, A] = SCA(x, .5/fs, 7/fs);
%
% % % Cx = cov(x');
% % % Cy = cov(y');
% % %
% % % [Vx Dx] = eig(Cx);
% % % [Vy Dy] = eig(Cy);
% %
% % figure
% % hold on
% % bar(diag(Dx),'b');
% % bar(diag(Dy),'r');
% % grid
%
% figure
% hold on
% for i = 1:size(x, 1),
%     plot(Hx(i).Frequencies, Hx(i).data/sqrt(mean(Hx(i).data.^2)),'b');
%     plot(Hy(i).Frequencies, Hy(i).data/sqrt(mean(Hy(i).data.^2)),'r');
% end
% grid;
%
% % PlotECG(x, 4, 'b', fsx);
% % PlotECG(y, 4, 'r', fsy);
%
% % PlotECG(y, 4, 'r', fs);
%
