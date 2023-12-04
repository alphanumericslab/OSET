% report all wavelet level results
clear;
close all;
clc;
addpath 'G:\Sameni\Documents\Papers\Journal 40 (Spline ECG Filtering)\CODES'

randn('seed', 0);
rand('seed', 0);

% % % dbstop if warning

type = {'PTB'};%%%{'Normal', 'Arrhythmia', 'PTB'};

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
    WNAME = {'haar', 'db2', 'db3' ,'db4', 'db5', 'db6', 'db7', 'db8', 'db9', 'db10', 'db12', 'db16', 'coif1', 'coif2', 'coif3', 'coif4', 'coif5', 'sym2', 'sym3', 'sym4', 'sym5', 'sym6', 'sym7', 'sym8', 'bior1.5', 'bior2.6', 'bior2.8', 'bior5.5', 'bior6.8'};
    
    % denoiserparam_list = logspace(-4, 2, 25);
    % denoiserparam_list = logspace(-4, 1, 50);
    % denoiserparam_list = logspace(-8, 1, 50);
    % % % denoiserparam_list = linspace(.01*(1/fs)^2, 100*(1/fs)^2, 50);
    %     wlen = 100e-3; % window length (s)
    trial_num = 1;
    fresults = ['AllDataSmootherResults' type{ttp} '10_ComparisonWithWavelet3.txt'];
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
                                            snr0 = 10*log10(mean(trial_snr0));
                                            snr1 = 10*log10(mean(trial_snr1));
                                            
%                                             strng = [num2str(kk) '\t' fname '\t' num2str(nn) '\t' num2str(dd) '\t' TPTR{tptr} '\t' SORH{sorh} '\t' SCAL{scal} '\t' WNAME{wname} '\t' num2str(WLEVELS(wlevel)) '\t' num2str(nt)  '\t' num2str(snr, '%4.4f') '\t' num2str(snr0, '%4.4f') '\t' num2str(snr1 - snr0, '%4.4f') '\n'];
                                            strng = [num2str(kk) '\t' fname '\t' num2str(dd) '\t' WNAME{wname} '\t' num2str(WLEVELS(wlevel)) '\t' num2str(snr, '%4.4f') '\t' num2str(snr1 - snr0, '%4.4f') '\n'];
                                            foutid = fopen(fout, 'a');
                                            fprintf(foutid, strng);
                                            fclose(foutid);
                                        end
                                    end
                                end
                            end
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
