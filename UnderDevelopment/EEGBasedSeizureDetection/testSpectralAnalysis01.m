close all;
clear;
clc;
resultfname = 'testSpectralAnalysis01.txt';

% parameters
subject = {'Dog_1', 'Dog_2', 'Dog_3', 'Dog_4', 'Dog_5', 'Patient_1', 'Patient_2'};
f0 = 60; % powerline frequency Hz
Q = 30; % notch filter Q-factor
wlen = [0.5 1.0 3.0 5.0 10.0]; % energy window lengths in seconds
drop = 200; % remove data head and tail due to filtering transitions
ffs = 160; % resample rate

nfft = 1000;
SpectralWinLen = 3.0; % in seconds
SpectralOverlapPercentage = 50; % [0 99]
SegmentLen = 180; % segment analyzed in seconds

% bins = 0:.01:3;
% epochlen = 15; % epoch length in seconds
% epochoverlap = 12; % overlap between epochs in seconds
% smoothing_iterations = 100;
x = [];
% estimate spectrum
h = spectrum.welch('Hamming', round(SpectralWinLen*ffs), SpectralOverlapPercentage); % Create a Welch spectral estimator.
% h = spectrum.periodogram('Hamming'); % Create a Periodogram spectral estimator.
% h = spectrum.periodogram('Rectangular');
% h = spectrum.yulear;
for m = 1 : length(subject),
    % search path and subjects
    path = ['E:\Sameni\Projects\Seizure\' subject{m} '\'];
    d = dir([path '*.mat']);
    NumRecords = length(d);
    for k =  1: NumRecords, % 475: 485%
        % clear previous variable to save memory
        if(~isempty(x))
            clear x;
        end
        
        % load data
        [x fs mode] = LoadSeizureEEGFull2(path, d(k).name);
        
        % resample data
        if(fs > ffs) % for the patients who have high sampling rates
            x = resample(x', ffs, round(fs))';
            fs = ffs;
        end
        
        % pre-process
        % x = LPFilter(x, 80.0/fs); % lowpass filter
        x = x - LPFilter(x, 0.25/fs); % highpass filter
        % IIR notch
        Wo = f0/(fs/2);  BW = Wo/Q;
        [b,a] = iirnotch(Wo,BW);
        x = filter(b, a, x, [], 2);
        
        % normalize
        x = (x - mean(x, 2)*ones(1, size(x,2)))./(std(x, [], 2)*ones(1, size(x,2)));
        
        % remove data head and tail due to transitions
        x1 = x(:, drop+1:drop+round(SegmentLen*fs));
        x2 = x(:, end-drop-round(SegmentLen*fs)+1:end-drop);
        
        H1 = MultiChannelSpectrum(x1, fs, nfft, h); % Estimate the PSD
        H2 = MultiChannelSpectrum(x2, fs, nfft, h); % Estimate the PSD
        %         H12 = MultiChannelSpectrum(x1.*x2, fs, nfft, h); % Estimate the PSD
        
        
        %         figure
        %         hold on
        S1 = zeros(size(x, 1), nfft/2);
        S2 = zeros(size(x, 1), nfft/2);
        %             S12 = zeros(size(x, 1), nfft/2);
        %         S12 = zeros(size(x, 1), nfft-1);
        for i = 1:size(x, 1),
            S1(i, :) = ( H1(i).data/sqrt(mean(H1(i).data.^2)) )';
            S2(i, :) = ( H2(i).data/sqrt(mean(H2(i).data.^2)) )';
            %             S12(i, :) = ( H12(i).data/sqrt(mean(H12(i).data.^2)) )';
            %             S12(i, :) = ( H12(i).data/sqrt(mean(H12(i).data.^2)) )';
            %                 H12 = MultiChannelSpectrum(xcorr(x1(i, :), x2(i, :)), fs, nfft, h); % Estimate the PSD
            %                 S12(i, :) = ( H12.data/sqrt(mean(H12.data.^2)) )';
            %             dev = (S1(i, :) - S2(i, :))./S1(i, :);
            %             plot(H12(i).Frequencies, dev);
            %             plot(H12(i).Frequencies, S12(i, :), 'g');
            %             S1(i, :) = ( H1(i).data )';
            %             S2(i, :) = ( H2(i).data )';
            %             plot(H1(i).Frequencies, S1(i, :), 'b');
            %             plot(H2(i).Frequencies, S2(i, :), 'r');
            %                 plot(H12.Frequencies, S12(i, :), 'g');
        end
        %         grid;
        
        absolute_spectral_variation = (S1 - S2)./S1;
        relative_spectral_variation = (S1 - S2)./S1;
        
        f1 = mean(relative_spectral_variation(:));
        f2 = median(relative_spectral_variation(:));
        f3 = std(relative_spectral_variation(:));
        
        f4 = mean(absolute_spectral_variation(:));
        f5 = median(absolute_spectral_variation(:));
        f6 = std(absolute_spectral_variation(:));
        
        % write results
        fid = fopen(resultfname,'a');
        fprintf(fid, '%s\t%d\t%d\t%d\t', d(k).name, m, k, mode);
        fprintf(fid, '%8.6f\t%8.6f\t%8.6f\t%8.6f\t%8.6f\t%8.6f\t', f1, f2, f3, f4, f5, f6);
        fprintf(fid, '\n');
        fclose(fid);
        disp(['subject: ', d(k).name]);
    end
end

% for i = 1 : size(x, 1)
%     figure;
%     t = (0:length(x)-1)/fs;
%     plot(t, x(i, :),'k');
%     hold on
%     t = (0:length(E)-1)/fs;
%     plot(t, E(i, :));
%     grid
% end
%
% figure
% imshow(Cx);
% grid
%
clock
