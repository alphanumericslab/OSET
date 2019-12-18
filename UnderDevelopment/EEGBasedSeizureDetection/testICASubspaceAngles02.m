close all;
clear;
clc;
resultfname = 'testICASubspaceAngles02.txt';

% parameters
subject = {'Dog_1', 'Dog_2', 'Dog_3', 'Dog_4', 'Dog_5', 'Patient_1', 'Patient_2'};
f0 = 60; % powerline frequency Hz
Q = 30; % notch filter Q-factor
% wlen = [0.5 1.0 3.0 5.0 10.0]; % energy window lengths in seconds
drop = 200; % remove data head and tail due to filtering transitions
ffs = 400;%160; % resample rate

nfft = 20;
SpectralWinLen = 3.0; % in seconds
SpectralOverlapPercentage = 50; % [0 99]
SegmentLen = 30; % segment analyzed in seconds
sourcenum = 8;

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
for m = 1 : 1,%length(subject),
    % search path and subjects
    path = ['E:\Sameni\Projects\Seizure\' subject{m} '\'];
    d = dir([path '*.mat']);
    NumRecords = length(d);
    for k = 489,%475: 485,%1: 2,NumRecords, %
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
        
%         % pre-process
%         % x = LPFilter(x, 80.0/fs); % lowpass filter
        x = x - LPFilter(x, 0.25/fs); % highpass filter
%         % IIR notch
%         Wo = f0/(fs/2);  BW = Wo/Q;
%         [b,a] = iirnotch(Wo,BW);
%         x = filter(b, a, x, [], 2);
        
        % normalize
        %         x = (x - mean(x, 2)*ones(1, size(x,2)))./(std(x, [], 2)*ones(1, size(x,2)));
        
        % remove data head and tail due to transitions
        x1 = x(:, drop+1:drop+round(SegmentLen*fs));
        x2 = x(:, end-drop-round(SegmentLen*fs)+1:end-drop);
        
        % PCA
        Cx1 = cov(x1');
        Cx2 = cov(x2');
        [V1 D1] = eig(Cx1);
        [V2 D2] = eig(Cx2);
        d1 = diag(D1);
        d2 = diag(D2);
        [Y1, I1] = sort(d1, 1, 'descend');
        [Y2, I2] = sort(d2, 1, 'descend');
        W1PCA = V1(:, I1(1:sourcenum))';
        W2PCA = V2(:, I2(1:sourcenum))';
        
        % JADE
        W1JADE = jadeR(x1, sourcenum);
        W2JADE = jadeR(x2, sourcenum);
        W1JADE = normr(W1JADE);
        W2JADE = normr(W2JADE);
        s1 = W1JADE*x1;
        s2 = W2JADE*x2;
        
%         % SOBI
%                 [~, s1] = sobi(x1, [], 500);
%                 [~, s2] = sobi(x2, [], 500);
        
        % SOBI
        %             w = pinv(sobi(xx));
        %             wSOBI(i,:,:) = w(1:sourcenum, :);
        
        % SCA
        %             [~, w] = SCA2(xx, 10.0/fs, 10.0/fs, 5);
        %             wSCA(i,:,:) = w(1:sourcenum, :);
        
        PlotECG(x1, 4, 'b', fs);
        PlotECG(s1, 4, 'm', fs);
        PlotECG(x2, 4, 'k', fs);
        PlotECG(s2, 4, 'r', fs);
        % % %         % write results
        % % %         fid = fopen(resultfname,'a');
        % % %         fprintf(fid, '%s\t%d\t%d\t%d\t', d(k).name, m, k, mode);
        % % %         fprintf(fid, '%8.6f\t%8.6f\t%8.6f\t%8.6f\t%8.6f\t%8.6f\t', f1, f2, f3, f4, f5, f6);
        % % %         fprintf(fid, '\n');
        % % %         fclose(fid);
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
