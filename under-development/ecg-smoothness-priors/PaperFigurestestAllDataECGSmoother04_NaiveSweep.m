clear;
close all;
clc;
addpath 'G:\Sameni\Documents\Papers\Journal 40 (Spline ECG Filtering)\CODES'

randn('seed', 0);
rand('seed', 0);

% % % dbstop if warning

type = {'Normal', 'Arrhythmia', 'PTB'};

for ttp = 1 : length(type),
    switch(type{ttp})
        case 'Normal'
            pth = 'J:\ECGData\Normal';
            fs = 128;
            segstart_list = 0 : 60 : 30*60; % segstart_list = 0 : 600 : 12*3600;
            diff_h = [2 4 6];%load DifferentiatorImpulseResponse01 diff_h
            seg_len = 10;
            smoothnessparam = logspace(-2.2, 4.2, 250);
            %             smoothnessparam = 0 : 0.1 : 15.0;
        case 'Arrhythmia'
            pth = 'J:\ECGData\Arrhythmia';
            fs = 360;
            segstart_list = 0 : 60 : 20*60;
            diff_h = [2 4 6];
            seg_len = 10;
            smoothnessparam = logspace(-2.2, 3.5, 250);
            %             smoothnessparam = 0 : 0.1 : 15.0;
        case 'PTB'
            pth = 'J:\ECGData\PTB';
            fs = 1000;
            segstart_list = [0 20];
            diff_h = [2 4 6];
            seg_len = 10;
            smoothnessparam = logspace(-2.2, 3.5, 250);
            %             smoothnessparam = 0 : 0.1 : 15.0;
    end
    
    noisetypes_list = {'WHITE'}; % 'EM' 'MA'
    snr_list = 30 : -3: -15;
    % denoiserparam_list = logspace(-4, 2, 25);
    % denoiserparam_list = logspace(-4, 1, 50);
    % denoiserparam_list = logspace(-8, 1, 50);
    % % % denoiserparam_list = linspace(.01*(1/fs)^2, 100*(1/fs)^2, 50);
    %     wlen = 100e-3; % window length (s)
    trial_num = 1;
    mode = 1;
    LCurveSweepLength = 1;%100;
    ForgettingFactor = 1;%%%0.9; % forgetting rate of the smoothness factor calculated in each block
    NoisePowerOverEstimationFactor = 1.0;
    fresults = ['AllDataSmootherResults' type{ttp} '04_NaiveSweepGammaBased.txt'];
    wd = cd;
    cd(pth);
    fout = fresults;
    foutid = fopen(fout, 'w');
    % strng = ['snr0 = ' num2str(snr0) '(dB), snrimp1 = ' num2str(snr1- snr0) '(dB), snrimp2 = ' num2str(snr2 - snr0) '(dB), denparam = ' num2str(denparam)];
    % fprintf(foutid, strng);
    fclose(foutid);
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
                for oo = 1 : length(diff_h),
                    data = data_all_channels(dd, :);
                    %             data = data - mean(data); % make the data zero mean
                    data = data - LPFilter(data, 1.0/fs); % we're not targeting baseline wanders at this stage
                    %             peaks = find(PeakDetection(data, 1.0/fs) == 1);
                    SignalPower = var(data);%mean(data.^2);
                    for nt = 1 : length(noisetypes_list),
                        for ss = 1 : length(snr_list),
                            snr = snr_list(ss);
                            NoisePower = SignalPower / 10^(snr/10);
                            snr0 = zeros(1, length(smoothnessparam));
                            snr1 = zeros(1, length(smoothnessparam));
                            %                         snr2 = zeros(1, length(smoothnessparam));
                            for pm = 1 : length(smoothnessparam),
                                trial_snr0 = zeros(1, trial_num);
                                trial_snr1 = zeros(1, trial_num);
                                %                             trial_snr2 = zeros(1, trial_num);
                                %                             trial_snr2 = zeros(1, trial_num);
                                %                             trial_snr3 = zeros(1, trial_num);
                                %                             trial_snr4 = zeros(1, trial_num);
                                %                     trial_snr4 = zeros(length(denoiserparam_list), trial_num);
                                for trial = 1 : trial_num, % trial is the most internal loop
                                    %                             noise = NoiseGenerator(noisetypes_list(nt), SignalPower, snr, length(data), fs)';
                                    noise =  NoiseGeneratorSimple(noisetypes_list{nt}, NoisePower, fs, length(data));
                                    x = data + noise;
                                    %smoothnessparam = NoisePowerOverEstimationFactor*NoisePower;
                                    %                             [~, x_filtered, ] = ECGSmoothnessPriorsDenoiserBW(x, smoothnessparam(pm), mode1, diff_h(oo), round(wlen*fs), 1e-8, 250, ForgettingFactor);
                                    %                             [~, x_smoothedKF, Pbar, Phat, PSmoothed, Kgain, alpha] = ECGSmoothnessPriorsDenoiserKF(x, diff_h(oo), NoisePower, (fs*wlen)^2/denoiserparam_list(ll));
                                    %                                 [~, x_filtered] = ECGSmoothnessPriorsDenoiserBW(x, smoothnessparam(pm), mode, diff_h(oo), round(wlen*fs), 1e-8, 250, ForgettingFactor);
                                    x_filtered_LTI = ECGSmoothnessPriorsDenoiserLTI(x, smoothnessparam(pm), mode, diff_h(oo));
                                    trial_snr0(trial) = sum(data.^2)/sum((data - x).^2);
                                    trial_snr1(trial) = sum(data.^2)/sum((data - x_filtered_LTI).^2);
                                    %                                 trial_snr2(trial) = sum(data.^2)/sum((data - x_filtered).^2);
                                    %                             trial_snr4(ll, trial) = sum(data.^2)/sum((data - x_smoothedKF).^2);
                                end
                                snr0(pm) = 10*log10(mean(trial_snr0));
                                snr1(pm) = 10*log10(mean(trial_snr1));
                                %                             snr2(pm) = 10*log10(mean(trial_snr2));
                                %                             snr2(pm) = 10*log10(mean(trial_snr2));
                                %                         snr4 = 10*log10(mean(trial_snr4(ll, :)));
                                %                         strng = [num2str(snr, '%4.4f') '\t' num2str(snr0, '%4.4f') '\t' num2str(snr1 - snr0, '%4.4f') '\t' num2str(snr2 - snr0, '%4.4f') '\t' num2str(snr3 - snr0, '%4.4f') '\t' num2str(snr4 - snr0, '%4.4f') '\t' num2str(denparam, '%2.6g') '\n'];
                                
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
                            [Y I] = max(snr1 - snr0);
                            strng = [num2str(kk) '\t' fname '\t' num2str(nn) '\t' num2str(dd) '\t' num2str(diff_h(oo)) '\t' num2str(nt)  '\t' num2str(ss) '\t' num2str(snr, '%4.4f') '\t' num2str(snr0(I), '%4.4f') '\t' num2str(Y, '%4.4f') '\t' num2str(smoothnessparam(I), '%2.6g') '\t' num2str(I, '%d') '\n'];
                            %                         strng = [num2str(kk) '\t' fname '\t' num2str(nn) '\t' num2str(dd)  '\t' num2str(nt)  '\t' num2str(ss) '\t' num2str(snr, '%4.4f') '\t' num2str(snr0(I), '%4.4f') '\t' num2str(snr1(I) - snr0(I), '%4.4f') '\t' num2str(snr2(I) - snr0(I), '%4.4f') '\t' num2str(smoothnessparam(I), '%2.6g') '\t' num2str(I, '%d') '\n'];
                            foutid = fopen(fout, 'a');
                            fprintf(foutid, strng);
                            fclose(foutid);
                            disp(['SNR = ' num2str(snr)]);
                        end
                    end
                end
                %                 delete([pth '\sampledata.txt']);
            end
        end
    end
end
cd(wd);

% system(['rdsamp -r ' fname ' -f ' num2str(segstart) ' -t ' num2str(segstop) ' -s ' channels ' > data.txt']);
