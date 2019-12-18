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
            segstart_list = 20;%0 : 60 : 30*60; % segstart_list = 0 : 600 : 12*3600;
            diff_h = [2 4 6];%load DifferentiatorImpulseResponse01 diff_h
            seg_len = 10;
        case 'Arrhythmia'
            pth = 'J:\ECGData\Arrhythmia';
            fs = 360;
            segstart_list = 20;%0 : 60 : 20*60;
            diff_h = [2 4 6];
            seg_len = 10;
        case 'PTB'
            pth = 'J:\ECGData\PTB';
            fs = 1000;
            segstart_list = 20;%[0 20];
            diff_h = [2 4 6];
            seg_len = 10;
    end
    
    noisetypes_list = {'WHITE'}; % 'EM' 'MA'
    snr_list = 30 : -3: -15;
    wlen = [100e-3 200e-3];%%%100e-3; % window length (s)
    trial_num = 1;
    NoisePowerOverEstimationFactor = [0.8 1 1.2];%%%sqrt(fs*wlen)];
    ForgettingFactor = 1;%%%0.9 : 0.1 : 1; % forgetting rate of the smoothness factor calculated in each block
    mode = 3;
    MAX_ITERATION = 100;
    ACCURACY = 1e-8;
    %     LCurveSweepLength = 1;%100;
    fresults = ['AllDataSmootherResults' type{ttp} '08_KnownVarianceOptimized.txt'];
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
            for dd = 1 : size(data_all_channels, 1),
                data = data_all_channels(dd, :);
                %             data = data - mean(data); % make the data zero mean
                data = data - LPFilter(data, 1.0/fs); % we're not targeting baseline wanders at this stage
                SignalPower = mean(data.^2);%var(data);
                for nt = 1 : length(noisetypes_list),
                    for ss = 1 : length(snr_list),
                        snr = snr_list(ss);
                        NoisePower = SignalPower / 10^(snr/10);
                        %                             noise = NoiseGenerator(noisetypes_list(nt), SignalPower, snr, length(data), fs)';
                        noise =  NoiseGeneratorSimple(noisetypes_list{nt}, NoisePower, fs, length(data));
                        x = data + noise;
                        snr0 = 10*log10(sum(data.^2)/sum((data - x).^2));
                        
                        for ww = 1 : length(wlen),
                            for oo = 1 : length(diff_h),
                                for ff = 1 : length(ForgettingFactor),
                                    for ovr = 1 : length(NoisePowerOverEstimationFactor),
                                        [x_filtered1, x_filtered2, optim_params1, optim_params2] = ECGSmoothnessPriorsDenoiserBW(x, NoisePowerOverEstimationFactor(ovr)*NoisePower, mode, diff_h(oo), round(wlen(ww)*fs), ACCURACY, MAX_ITERATION, ForgettingFactor(ff));
                                        
                                        snr1 = 10*log10(sum(data.^2)/sum((data - x_filtered1).^2));
                                        snr2 = 10*log10(sum(data.^2)/sum((data - x_filtered2).^2));
                                        
                                        strng = [num2str(kk) '\t' fname '\t' num2str(nn) '\t' num2str(dd) '\t' num2str(nt)  '\t' num2str(ss) '\t' num2str(ww) '\t' num2str(diff_h(oo)) '\t' num2str(ff) '\t' num2str(ovr) '\t' num2str(snr, '%4.4f') '\t' num2str(snr0, '%4.4f') '\t' num2str(snr1 - snr0, '%2.6g') '\t' num2str(snr2 - snr0, '%4.4f') '\t' num2str(nanmean(1./optim_params1), '%d') '\t' num2str(nanmean(1./optim_params2), '%d') '\n'];
                                        %                         strng = [num2str(kk) '\t' fname '\t' num2str(nn) '\t' num2str(dd)  '\t' num2str(nt)  '\t' num2str(ss) '\t' num2str(snr, '%4.4f') '\t' num2str(snr0(I), '%4.4f') '\t' num2str(snr1(I) - snr0(I), '%4.4f') '\t' num2str(snr2(I) - snr0(I), '%4.4f') '\t' num2str(smoothnessparam(I), '%2.6g') '\t' num2str(I, '%d') '\n'];
                                        foutid = fopen(fout, 'a');
                                        fprintf(foutid, strng);
                                        fclose(foutid);
                                        
                                        % % %                                     figure
                                        % % %                                     hold on
                                        % % %                                     plot(x, 'c');
                                        % % %                                     plot(x_filtered1, 'b', 'linewidth', 2);
                                        % % %                                     plot(x_filtered2, 'r', 'linewidth', 2);
                                        % % %                                     plot(data, 'k');
                                        % % %                                     grid
                                        % % %                                     title(['SNR = ' num2str(snr) ', FF = '  num2str(ForgettingFactor(ff)) ', NPowerOVR = '  num2str(NoisePowerOverEstimationFactor(ovr))])
                                    end
                                end
                                %                             disp(['SNR = ' num2str(snr)]);
                            end
                        end
                    end
                end
                %                 delete([pth '\sampledata.txt']);
            end
        end
        disp(fname);
    end
end
cd(wd);

% system(['rdsamp -r ' fname ' -f ' num2str(segstart) ' -t ' num2str(segstop) ' -s ' channels ' > data.txt']);
