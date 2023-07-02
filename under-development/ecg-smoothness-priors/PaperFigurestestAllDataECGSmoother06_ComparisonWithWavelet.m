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
            %             segstart_list = 0 : 60 : 30*60; % segstart_list = 0 : 600 : 12*3600;
            segstart_list = 20; % for test
            seg_len = 10;
        case 'Arrhythmia'
            pth = 'J:\ECGData\Arrhythmia';
            fs = 360;
            %             segstart_list = 0 : 60 : 20*60;
            segstart_list = 20; % for test
            seg_len = 10;
        case 'PTB'
            pth = 'J:\ECGData\PTB';
            fs = 1000;
            %             segstart_list = [0 20];
            segstart_list = 20; % for test
            seg_len = 10;
    end
    
    noisetypes_list = {'WHITE'}; % 'EM' 'MA'
    snr_list = 30 : -3: -15;
    TPTR = {'rigrsure'};%%%{'rigrsure', 'heursure', 'sqtwolog', 'minimaxi'};
    SORH = {'s'};%%%{'s', 'h'};
    SCAL = {'sln'};%%%{'one', 'sln', 'mln'};
    WLEVELS = 1 : 10;
    WNAME = {'haar', 'db2', 'db3' ,'db4', 'db5', 'db6', 'db8', 'db12', 'db16', 'coif2', 'coif3', 'coif4', 'coif5', 'sym4', 'sym5', 'sym6', 'sym7', 'sym8', 'bior1.5', 'bior2.6', 'bior2.8', 'bior5.5', 'bior6.8'};
    %%% 'dmey', 'meyr', 'morl', 'mexh', 'gaus???', 'cgau???', 'shan???', 'fbsp???', 'cmor???', 
    
    % denoiserparam_list = logspace(-4, 2, 25);
    % denoiserparam_list = logspace(-4, 1, 50);
    % denoiserparam_list = logspace(-8, 1, 50);
    % % % denoiserparam_list = linspace(.01*(1/fs)^2, 100*(1/fs)^2, 50);
    %     wlen = 100e-3; % window length (s)
    trial_num = 1;
    fresults = ['AllDataSmootherResults' type{ttp} '06_ComparisonWithWavelet.txt'];
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
                data = data_all_channels(dd, :);
                data = data - LPFilter(data, 1.0/fs); % we're not targeting baseline wanders at this stage
                SignalPower = var(data);
                for nt = 1 : length(noisetypes_list),
                    for ss = 1 : length(snr_list),
                        snr = snr_list(ss);
                        NoisePower = SignalPower / 10^(snr/10);
                        
                        snr0 = zeros(length(TPTR), length(SORH), length(SCAL), length(WNAME), length(WLEVELS));
                        snr1 = snr0;
                        SYZ = size(snr0);
                        for tptr = 1 : length(TPTR),
                            for sorh = 1 : length(SORH),
                                for scal = 1 : length(SCAL),
                                    for wname = 1 : length(WNAME),
                                        for wlevel = 1 : length(WLEVELS),
                                            trial_snr0 = zeros(1, trial_num);
                                            trial_snr1 = zeros(1, trial_num);
                                            for trial = 1 : trial_num, % trial is the most internal loop
                                                %                             noise = NoiseGenerator(noisetypes_list(nt), SignalPower, snr, length(data), fs)';
                                                noise =  NoiseGeneratorSimple(noisetypes_list{nt}, NoisePower, fs, length(data));
                                                x = data + noise;
                                                x_filtered_Wavelet = wden(x, TPTR{tptr}, SORH{sorh}, SCAL{scal}, WLEVELS(wlevel), WNAME{wname});
                                                
                                                trial_snr0(trial) = sum(data.^2)/sum((data - x).^2);
                                                trial_snr1(trial) = sum(data.^2)/sum((data - x_filtered_Wavelet).^2);
                                            end
                                            snr0(tptr, sorh, scal, wname, wlevel) = 10*log10(mean(trial_snr0));
                                            snr1(tptr, sorh, scal, wname, wlevel) = 10*log10(mean(trial_snr1));
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
                                end
                            end
                        end
                        [YMAX, IMAX] = max(snr1(:) - snr0(:));
                        [I1, I2, I3, I4, I5] = ind2sub(SYZ, IMAX);
                        %                         [~, I1] = max(snr1 - snr0, [], 1);
                        %                         [~, I2] = max(snr1 - snr0, [], 2);
                        %                         [~, I3] = max(snr1 - snr0, [], 3);
                        %                         [~, I4] = max(snr1 - snr0, [], 4);
                        %                         [~, I5] = max(snr1 - snr0, [], 5);
                        strng = [num2str(kk) '\t' fname '\t' num2str(nn) '\t' num2str(dd) '\t' TPTR{I1} '\t' SORH{I2} '\t' SCAL{I3} '\t' WNAME{I4} '\t' num2str(WLEVELS(I5)) '\t' num2str(nt)  '\t' num2str(snr, '%4.4f') '\t' num2str(snr0(IMAX), '%4.4f') '\t' num2str(YMAX, '%4.4f') '\n'];
                        %                         strng = [num2str(kk) '\t' fname '\t' num2str(nn) '\t' num2str(dd)  '\t' num2str(nt)  '\t' num2str(ss) '\t' num2str(snr, '%4.4f') '\t' num2str(snr0(I), '%4.4f') '\t' num2str(snr1(I) - snr0(I), '%4.4f') '\t' num2str(snr2(I) - snr0(I), '%4.4f') '\t' num2str(NormalizedCutOff(I), '%2.6g') '\t' num2str(I, '%d') '\n'];
                        foutid = fopen(fout, 'a');
                        fprintf(foutid, strng);
                        fclose(foutid);
                        disp(['SNR = ' num2str(snr)]);
                    end
                end
            end
        end
        %                 delete([pth '\sampledata.txt']);
    end
end
cd(wd);

% system(['rdsamp -r ' fname ' -f ' num2str(segstart) ' -t ' num2str(segstop) ' -s ' channels ' > data.txt']);
