clear;
close all;
clc;
addpath 'G:\Sameni\Documents\Papers\Journal 40 (Spline ECG Filtering)\CODES'

type = 'PTB';

switch(type)
    case 'Normal'
        pth = 'J:\ECGData\Normal';
        fs = 128;
        segstart_list = 0 : 60 : 30*60; % segstart_list = 0 : 600 : 12*3600;
        denoiserparam_list = logspace(-8, 3, 100);
        diff_h = 2;%load DifferentiatorImpulseResponse01 diff_h
        seg_len = 10;
    case 'Arrhythmia'
        pth = 'J:\ECGData\Arrhythmia';
        fs = 360;
        segstart_list = 0 : 60 : 20*60;
        diff_h = 2;
        denoiserparam_list = logspace(-3, 3, 100);
        seg_len = 10;
    case 'PTB'
        pth = 'J:\ECGData\PTB';
        fs = 1000;
        segstart_list = [0 20];
        denoiserparam_list = logspace(-4, 6, 100);
        diff_h = 2;
        seg_len = 10;
end

noisetypes_list = {'WHITE'}; % 'EM' 'MA'
snr_list = 25 : -5: -5;
% denoiserparam_list = logspace(-4, 2, 25);
% denoiserparam_list = logspace(-4, 1, 50);
% denoiserparam_list = logspace(-8, 1, 50);
% % % denoiserparam_list = linspace(.01*(1/fs)^2, 100*(1/fs)^2, 50);
wlen = 150e-3; % window length (s)
trial_num = 1;
mode = 5;
fresults = ['AllDataSmootherResults' type '01.txt'];
fout = [fresults];
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
        %         system(['rdsamp -r ' fname ' > sampledata.txt']);
        data_all_channels = load([pth '\sampledata.txt'])';
        data_all_channels = data_all_channels(2:end,:);
        %         figure
        %         PlotECG(data_all_channels, 2, 'b', fs);
        for dd = 1 : size(data_all_channels, 1),
            data = data_all_channels(dd, :);
            %             data = data - mean(data); % make the data zero mean
            data = data - LPFilter(data, 1.0/fs); % we're not targeting baseline wanders at this stage
            SignalPower = var(data);%mean(data.^2);
            for nt = 1 : length(noisetypes_list),
                for ss = 1 : length(snr_list),
                    snr = snr_list(ss);
                    trial_snr0 = zeros(length(denoiserparam_list), trial_num);
                    trial_snr1 = zeros(length(denoiserparam_list), trial_num);
                    trial_snr2 = zeros(length(denoiserparam_list), trial_num);
                    trial_snr3 = zeros(length(denoiserparam_list), trial_num);
                    %                     trial_snr4 = zeros(length(denoiserparam_list), trial_num);
                    for ll = 1 : length(denoiserparam_list), % the algorithmic parameters are in the second loop
                        denparam = denoiserparam_list(ll);
                        for trial = 1 : trial_num, % trial is the most internal loop
                            %                             noise = NoiseGenerator(noisetypes_list(nt), SignalPower, snr, length(data), fs)';
                            NoisePower = SignalPower / 10^(snr/10);
                            noise =  NoiseGeneratorSimple(noisetypes_list{nt}, NoisePower, fs, length(data));
                            x = data + noise;
                            [~, x_filtered] = ECGSmoothnessPriorsDenoiserBW(x, denparam, 0, diff_h, round(wlen*fs));
                            [~, x_filtered_adaptive] = ECGSmoothnessPriorsDenoiserBW(x, denparam, 5, diff_h, round(wlen*fs));
                            x_filtered_LTI = ECGSmoothnessPriorsDenoiserLTI(x, denparam, [], diff_h);
                            %                             [~, x_smoothedKF, Pbar, Phat, PSmoothed, Kgain, alpha] = ECGSmoothnessPriorsDenoiserKF(x, diff_h, NoisePower, (fs*wlen)^2/denoiserparam_list(ll));
                            trial_snr0(ll, trial) = sum(data.^2)/sum((data - x).^2);
                            trial_snr1(ll, trial) = sum(data.^2)/sum((data - x_filtered).^2);
                            trial_snr2(ll, trial) = sum(data.^2)/sum((data - x_filtered_adaptive).^2);
                            trial_snr3(ll, trial) = sum(data.^2)/sum((data - x_filtered_LTI).^2);
                            %                             trial_snr4(ll, trial) = sum(data.^2)/sum((data - x_smoothedKF).^2);
                        end
                        snr0 = 10*log10(mean(trial_snr0(ll, :)));
                        snr1 = 10*log10(mean(trial_snr1(ll, :)));
                        snr2 = 10*log10(mean(trial_snr2(ll, :)));
                        snr3 = 10*log10(mean(trial_snr3(ll, :)));
                        %                         snr4 = 10*log10(mean(trial_snr4(ll, :)));
                        %                         strng = [num2str(snr, '%4.4f') '\t' num2str(snr0, '%4.4f') '\t' num2str(snr1 - snr0, '%4.4f') '\t' num2str(snr2 - snr0, '%4.4f') '\t' num2str(snr3 - snr0, '%4.4f') '\t' num2str(snr4 - snr0, '%4.4f') '\t' num2str(denparam, '%2.6g') '\n'];
                        strng = [num2str(snr, '%4.4f') '\t' num2str(snr0, '%4.4f') '\t' num2str(snr1 - snr0, '%4.4f') '\t' num2str(snr2 - snr0, '%4.4f') '\t' num2str(snr3 - snr0, '%4.4f') '\t' num2str(denparam, '%2.6g') '\n'];
                        foutid = fopen(fout, 'a');
                        fprintf(foutid, strng);
                        fclose(foutid);
                        disp(num2str(ll));
                    end
                    all_params_snr1 = 10*log10(mean(trial_snr1, 2));
                    all_params_snr2 = 10*log10(mean(trial_snr2, 2));
                    all_params_snr3 = 10*log10(mean(trial_snr3, 2));
                    %                     all_params_snr4 = 10*log10(mean(trial_snr4, 2));
                    [SNR1_max, I1] = max(all_params_snr1);
                    [SNR2_max, I2] = max(all_params_snr2);
                    [SNR3_max, I3] = max(all_params_snr3);
                    %                     [SNR4_max, I4] = max(all_params_snr4);
                    %                     noise = NoiseGenerator(noisetypes_list(nt), SignalPower, snr, length(data), fs)';
                    noise =  NoiseGeneratorSimple(noisetypes_list{nt}, NoisePower, fs, length(data));
                    x = data + noise;
                    [~, x_filtered_best] = ECGSmoothnessPriorsDenoiserBW(x, denoiserparam_list(I1), 0, diff_h, round(wlen*fs));
                    [~, x_filtered_adaptive_best] = ECGSmoothnessPriorsDenoiserBW(x, denoiserparam_list(I2), 5, diff_h, round(wlen*fs));
                    x_filtered_LTI_best = ECGSmoothnessPriorsDenoiserLTI(x, denoiserparam_list(I3), [], diff_h);
                    %                     [~, x_smoothedKF_best, Pbar, Phat, PSmoothed, Kgain, alpha] = ECGSmoothnessPriorsDenoiserKF(x, diff_h, NoisePower, (fs*wlen)^2/denoiserparam_list(I4));
                    t = (0:length(data)-1)/fs;
                    figure
                    hold on
                    plot(t, x, 'c', 'linewidth', 1);
                    plot(t, x_filtered_best, 'b', 'linewidth', 2);
                    plot(t, x_filtered_adaptive_best, 'r', 'linewidth', 2);
                    plot(t, x_filtered_LTI_best, 'm', 'linewidth', 2);
                    %                     plot(t, x_smoothedKF_best, 'c', 'linewidth', 2);
                    plot(t, data, 'k', 'linewidth', 1);
                    %                     title(['SNR0 = ' num2str(snr) '(dB), SNR1imp = ' num2str(SNR1_max - snr) '(dB), SNR2imp = ' num2str(SNR2_max - snr) '(dB), SNR3imp = ' num2str(SNR3_max - snr) '(dB), SNR4imp = ' num2str(SNR4_max - snr) '(dB), Params = [' num2str(denoiserparam_list(I1)) ', ' num2str(denoiserparam_list(I2)) ', ' num2str(denoiserparam_list(I3)) ', ' num2str(denoiserparam_list(I4)) ']']);
                    title(['SNR0 = ' num2str(snr) '(dB), SNR1imp = ' num2str(SNR1_max - snr) '(dB), SNR2imp = ' num2str(SNR2_max - snr) '(dB), SNR3imp = ' num2str(SNR3_max - snr) '(dB), Params = [' num2str(denoiserparam_list(I1)) ', ' num2str(denoiserparam_list(I2)) ', ' num2str(denoiserparam_list(I3)) ']']);
                    grid
                end
                %                 delete([pth '\sampledata.txt']);
            end
        end
    end
end
cd(wd);

% system(['rdsamp -r ' fname ' -f ' num2str(segstart) ' -t ' num2str(segstop) ' -s ' channels ' > data.txt']);
