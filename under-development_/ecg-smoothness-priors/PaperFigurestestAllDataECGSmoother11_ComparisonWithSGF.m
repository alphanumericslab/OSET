% compare with this paper:
% S. Krishnan and C. Seelamantula, “On the selection of optimum savitzky-golay filters,” Signal Processing, IEEE Transactions on, vol. 61, no. 2, pp. 380–391, Jan 2013.

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
    WLEN = 5:2:75;
    
    % denoiserparam_list = logspace(-4, 2, 25);
    % denoiserparam_list = logspace(-4, 1, 50);
    % denoiserparam_list = logspace(-8, 1, 50);
    % % % denoiserparam_list = linspace(.01*(1/fs)^2, 100*(1/fs)^2, 50);
    %     wlen = 100e-3; % window length (s)
    trial_num = 1;
    fresults = ['AllDataSmootherResults' type{ttp} '11_ComparisonWithGSF.txt'];
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
                        
                        for winlen = 1 : length(WLEN),
                            trial_snr0 = zeros(1, trial_num);
                            trial_snr1 = zeros(1, trial_num);
                            trial_snr2 = zeros(1, trial_num);
                            trial_snr3 = zeros(1, trial_num);
                            for trial = 1 : trial_num, % trial is the most internal loop
                                %                             noise = NoiseGenerator(noisetypes_list(nt), SignalPower, snr, length(data), fs)';
                                noise =  NoiseGeneratorSimple(noisetypes_list{nt}, NoisePower, fs, length(data));
                                x = data + noise;
                                
                                [x_filtered1, ~] = SGBandwidth(x, WLEN(winlen));
                                [x_filtered2, ~] = SGBandwidthReg(x, WLEN(winlen));
                                x_filtered3 = SGOrder(x, WLEN(winlen));
                                
                                trial_snr0(trial) = sum(data.^2)/sum((data - x).^2);
                                trial_snr1(trial) = sum(data.^2)/sum((data - x_filtered1).^2);
                                trial_snr2(trial) = sum(data.^2)/sum((data - x_filtered2).^2);
                                trial_snr3(trial) = sum(data.^2)/sum((data - x_filtered3).^2);
                            end
                            snr0 = 10*log10(mean(trial_snr0));
                            snr1 = 10*log10(mean(trial_snr1));
                            snr2 = 10*log10(mean(trial_snr2));
                            snr3 = 10*log10(mean(trial_snr3));
                            
                            % % %                             figure
                            % % %                             hold on
                            % % %                             plot(x);
                            % % %                             plot(x_filtered1, 'r');
                            % % %                             plot(x_filtered2, 'g');
                            % % %                             plot(x_filtered3, 'm');
                            % % %                             grid
                            
                            strng = [num2str(kk) '\t' fname '\t' num2str(nn) '\t' num2str(dd) '\t' num2str(nt)  '\t' num2str(WLEN(winlen)) '\t' num2str(snr, '%4.4f') '\t' num2str(snr0, '%4.4f') '\t' num2str(snr1 - snr0, '%4.4f') '\t' num2str(snr2 - snr0, '%4.4f') '\t' num2str(snr3 - snr0, '%4.4f') '\n']; % '\t' num2str(optwin1) '\t' num2str(optwin2)
                            foutid = fopen(fout, 'a');
                            fprintf(foutid, strng);
                            fclose(foutid);
                        end
                    end
                end
            end
        end
        disp(['Subject: ' fname]);
        %         pack
    end
end
cd(wd);

% system(['rdsamp -r ' fname ' -f ' num2str(segstart) ' -t ' num2str(segstop) ' -s ' channels ' > data.txt']);
