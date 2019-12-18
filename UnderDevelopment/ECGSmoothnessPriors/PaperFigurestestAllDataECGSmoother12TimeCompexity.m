% compare with this paper:
% S. Krishnan and C. Seelamantula, “On the selection of optimum savitzky-golay filters,” Signal Processing, IEEE Transactions on, vol. 61, no. 2, pp. 380–391, Jan 2013.

clear;
close all;
clc;
addpath 'G:\Sameni\Documents\Papers\Journal 40 (Spline ECG Filtering)\CODES'

randn('seed', 0);
rand('seed', 0);

% % % dbstop if warning

type = {'Normal'};%, 'Arrhythmia', 'PTB'};

% profile on
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
    snr_list = 15;%%%30 : -3: -15;
    
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
        for nn = 1 : 1;%%%length(segstart_list),
            fname = lst(kk).name;
            fname = fname(1:end-4);
            segstart = segstart_list(nn);
            segstop = segstart_list(nn) + seg_len;
            system(['rdsamp -r ' fname ' -f ' num2str(segstart) ' -t ' num2str(segstop) ' > sampledata.txt']);
            data_all_channels = load([pth '\sampledata.txt'])';
            data_all_channels = data_all_channels(2:end,:);
            %         figure
            %         PlotECG(data_all_channels, 2, 'b', fs);
            for dd = 1 : 1;%%%size(data_all_channels, 1),
                data = data_all_channels(dd, :);
                data = data - LPFilter(data, 1.0/fs); % we're not targeting baseline wanders at this stage
                SignalPower = var(data);
                for nt = 1 : 1;%%%length(noisetypes_list),
                    for ss = 1 : 1;%%%length(snr_list),
                        snr = snr_list(ss);
                        NoisePower = SignalPower / 10^(snr/10);
                        
                        noise =  NoiseGeneratorSimple(noisetypes_list{nt}, NoisePower, fs, length(data));
                        x = data + noise;
                        
                        [x_filtered1, ~] = SGBandwidth(x, 49);
                        [x_filtered2, ~] = SGBandwidthReg(x, 49);
                        x_filtered3 = SGOrder(x, 49);
                        [x_filteredBW1, x_filteredBW2] = ECGSmoothnessPriorsDenoiserBW(x, NoisePower, 3, 2, round(100e-3*fs), 1e-8, 100, 1.0);
                        [x_filteredBW3, x_filteredBW4] = ECGSmoothnessPriorsDenoiserBW(x, NoisePower, 1, 2, round(100e-3*fs), 1e-8, 100, 1.0);
                        x_LTI = ECGSmoothnessPriorsDenoiserLTI(x, NoisePower, 1, 2);
                        x_filtered_Wavelet = wden(x, 'rigrsure', 's', 'sln', 10, 'coif5');
                        
                    end
                end
            end
        end
        disp(['Subject: ' fname]);
    end
    %         pack
end
% profile off
% profile viewer
cd(wd);

% system(['rdsamp -r ' fname ' -f ' num2str(segstart) ' -t ' num2str(segstop) ' -s ' channels ' > data.txt']);
