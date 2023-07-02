clear;
close all;
clc;
addpath 'G:\Sameni\Documents\Papers\Journal 40 (Spline ECG Filtering)\CODES'

% % % dbstop if warning

type = 'Arrhythmia';

switch(type)
    case 'Normal'
        pth = 'J:\ECGData\Normal';
        fs = 128;
        segstart_list = 0 : 60 : 30*60; % segstart_list = 0 : 600 : 12*3600;
        diff_h = 2;%load DifferentiatorImpulseResponse01 diff_h
        seg_len = 10;
    case 'Arrhythmia'
        pth = 'J:\ECGData\Arrhythmia';
        fs = 360;
        segstart_list = 0 : 60 : 20*60;
        diff_h = 2;
        seg_len = 10;
    case 'PTB'
        pth = 'J:\ECGData\PTB';
        fs = 1000;
        segstart_list = [0 20];
        diff_h = 2;
        seg_len = 10;
end

noisetypes_list = {'WHITE'}; % 'EM' 'MA'
snr_list = 30 : -3: -12;
% denoiserparam_list = logspace(-4, 2, 25);
% denoiserparam_list = logspace(-4, 1, 50);
% denoiserparam_list = logspace(-8, 1, 50);
% % % denoiserparam_list = linspace(.01*(1/fs)^2, 100*(1/fs)^2, 50);
wlen = 100e-3; % window length (s)
trial_num = 1;
mode1 = 3;
mode2 = 7;
LCurveSweepLength = 100;
ForgettingFactor = 1;%%%0.9; % forgetting rate of the smoothness factor calculated in each block
NoisePowerOverEstimationFactor = 1.0;
fresults = ['AllDataSmootherResults' type '02.txt'];
fout = fresults;
foutid = fopen(fout, 'w');
% strng = ['snr0 = ' num2str(snr0) '(dB), snrimp1 = ' num2str(snr1- snr0) '(dB), snrimp2 = ' num2str(snr2 - snr0) '(dB), denparam = ' num2str(denparam)];
% fprintf(foutid, strng);
fclose(foutid);

randn('seed', 0);
rand('seed', 0);
wd = cd;
cd(pth);
lst = dir('*.dat');
for kk = 1 : length(lst)
    for nn = 1 : length(segstart_list),
        fname = lst(kk).name;
        fname = fname(1:end-4);
        segstart = segstart_list(nn);
        segstop = segstart_list(nn) + seg_len;
        system(['rdsamp -r ' fname ' -f ' num2str(segstart) ' -t ' num2str(segstop) ' > sampledata.txt']);
        data_all_channels = load([pth '\sampledata.txt'])';
        data_all_channels = data_all_channels(2:end,:);
        %         figure
        %         PlotECG(data_all_channels, 2, 'b', fs);
        for dd = 1 : size(data_all_channels, 1),
            data = data_all_channels(dd, :);
            %             data = data - mean(data); % make the data zero mean
            data = data - LPFilter(data, 1.0/fs); % we're not targeting baseline wanders at this stage
            %             peaks = find(PeakDetection(data, 1.0/fs) == 1);
            SignalPower = var(data);%mean(data.^2);
            for nt = 1 : length(noisetypes_list),
                for ss = 1 : length(snr_list),
                    snr = snr_list(ss);
                    trial_snr0 = zeros(1, trial_num);
                    trial_snr1 = zeros(1, trial_num);
                    trial_snr2 = zeros(1, trial_num);
                    trial_snr3 = zeros(1, trial_num);
                    trial_snr4 = zeros(1, trial_num);
                    %                     trial_snr4 = zeros(length(denoiserparam_list), trial_num);
                    for trial = 1 : trial_num, % trial is the most internal loop
                        %                             noise = NoiseGenerator(noisetypes_list(nt), SignalPower, snr, length(data), fs)';
                        NoisePower = SignalPower / 10^(snr/10);
                        noise =  NoiseGeneratorSimple(noisetypes_list{nt}, NoisePower, fs, length(data));
                        x = data + noise;
                        smoothnessparam = NoisePowerOverEstimationFactor*NoisePower;
                        [~, x_filtered, ~, optim_gammas2] = ECGSmoothnessPriorsDenoiserBW(x, smoothnessparam, mode1, diff_h, round(wlen*fs), 1e-8, 250, ForgettingFactor);
                        x_filtered_LTI1 = ECGSmoothnessPriorsDenoiserLTI(x, max(optim_gammas2), [], diff_h);
                        x_filtered_LTI2 = ECGSmoothnessPriorsDenoiserLTI(x, median(optim_gammas2), [], diff_h);
                        [~, x_filteredL, ~, optim_gammas2L] = ECGSmoothnessPriorsDenoiserBW(x, smoothnessparam, mode2, diff_h, round(wlen*fs), 1e-8, 250, ForgettingFactor, LCurveSweepLength);
                        %                             [~, x_smoothedKF, Pbar, Phat, PSmoothed, Kgain, alpha] = ECGSmoothnessPriorsDenoiserKF(x, diff_h, NoisePower, (fs*wlen)^2/denoiserparam_list(ll));
                        trial_snr0(trial) = sum(data.^2)/sum((data - x).^2);
                        trial_snr1(trial) = sum(data.^2)/sum((data - x_filtered).^2);
                        trial_snr2(trial) = sum(data.^2)/sum((data - x_filtered_LTI1).^2);
                        trial_snr3(trial) = sum(data.^2)/sum((data - x_filtered_LTI2).^2);
                        trial_snr4(trial) = sum(data.^2)/sum((data - x_filteredL).^2);
                        %                             trial_snr4(ll, trial) = sum(data.^2)/sum((data - x_smoothedKF).^2);
                    end
                    snr0 = 10*log10(mean(trial_snr0));
                    snr1 = 10*log10(mean(trial_snr1));
                    snr2 = 10*log10(mean(trial_snr2));
                    snr3 = 10*log10(mean(trial_snr3));
                    snr4 = 10*log10(mean(trial_snr4));
                    %                         snr4 = 10*log10(mean(trial_snr4(ll, :)));
                    %                         strng = [num2str(snr, '%4.4f') '\t' num2str(snr0, '%4.4f') '\t' num2str(snr1 - snr0, '%4.4f') '\t' num2str(snr2 - snr0, '%4.4f') '\t' num2str(snr3 - snr0, '%4.4f') '\t' num2str(snr4 - snr0, '%4.4f') '\t' num2str(denparam, '%2.6g') '\n'];
                    strng = [num2str(kk) '\t' fname '\t' num2str(nn) '\t' num2str(dd)  '\t' num2str(nt)  '\t' num2str(ss) '\t' num2str(snr, '%4.4f') '\t' num2str(snr0, '%4.4f') '\t' num2str(snr1 - snr0, '%4.4f') '\t' num2str(snr2 - snr0, '%4.4f') '\t' num2str(snr3 - snr0, '%4.4f') '\t' num2str(snr4 - snr0, '%4.4f') '\t' num2str(mean(optim_gammas2), '%2.6g') '\t' num2str(median(optim_gammas2), '%2.6g') '\t' num2str(min(optim_gammas2), '%2.6g') '\t' num2str(max(optim_gammas2), '%2.6g') '\t' num2str(mean(optim_gammas2L), '%2.6g') '\n'];
                    foutid = fopen(fout, 'a');
                    fprintf(foutid, strng);
                    fclose(foutid);
                    
                    disp(['SNR = ' num2str(snr)]);
                    % [~, x_smoothedKF_best, Pbar, Phat, PSmoothed, Kgain, alpha] = ECGSmoothnessPriorsDenoiserKF(x, diff_h, NoisePower, (fs*wlen)^2/denoiserparam_list(I4));
% % %                     t = (0:length(data)-1)/fs;
% % %                     figure
% % %                     hold on
% % %                     plot(t, x, 'c', 'linewidth', 1);
% % %                     plot(t, x_filtered, 'b', 'linewidth', 2);
% % %                     plot(t, x_filtered_LTI1, 'r', 'linewidth', 2);
% % %                     plot(t, x_filtered_LTI2, 'm', 'linewidth', 2);
% % %                     %                     plot(t, x_smoothedKF_best, 'c', 'linewidth', 2);
% % %                     plot(t, data, 'k', 'linewidth', 1);
% % %                     %                     title(['SNR0 = ' num2str(snr) '(dB), SNR1imp = ' num2str(SNR1_max - snr) '(dB), SNR2imp = ' num2str(SNR2_max - snr) '(dB), SNR3imp = ' num2str(SNR3_max - snr) '(dB), SNR4imp = ' num2str(SNR4_max - snr) '(dB), Params = [' num2str(denoiserparam_list(I1)) ', ' num2str(denoiserparam_list(I2)) ', ' num2str(denoiserparam_list(I3)) ', ' num2str(denoiserparam_list(I4)) ']']);
% % %                     title(['SNR0 = ' num2str(snr) '(dB), SNR1imp = ' num2str(snr1 - snr) '(dB), SNR2imp = ' num2str(snr2 - snr) '(dB), SNR3imp = ' num2str(snr3 - snr) '(dB), Params = [' num2str(max(optim_gammas2)) ', ' num2str(median(optim_gammas2)) ']']);
% % %                     grid
                end
            end
            %                 delete([pth '\sampledata.txt']);
        end
    end
end
cd(wd);

% system(['rdsamp -r ' fname ' -f ' num2str(segstart) ' -t ' num2str(segstop) ' -s ' channels ' > data.txt']);
