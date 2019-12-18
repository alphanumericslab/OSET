clear;
close all;
clc;
addpath 'G:\Sameni\Documents\Papers\Journal 40 (Spline ECG Filtering)\CODES'

randn('seed', 0);
rand('seed', 0);

% % % dbstop if warning

type = {'Normal', 'Arrhythmia', 'PTB'};

INITIAL_PARAM_LIST_NUM = 1000;
MAX_ITR = 100;
for ttp = 1 : length(type),
    switch(type{ttp})
        case 'Normal'
            pth = 'J:\ECGData\Normal';
            fs = 128;
            segstart_list = 20;%0 : 60 : 30*60; % segstart_list = 0 : 600 : 12*3600;
            diff_h = [2 4 6 8];%load DifferentiatorImpulseResponse01 diff_h
            seg_len = 10;
            smoothnessparam_min = -12;
            smoothnessparam_max = 6;
            %             smoothnessparam = logspace(-2.2, 4.2, 250);
            %             smoothnessparam = 0 : 0.1 : 15.0;
        case 'Arrhythmia'
            pth = 'J:\ECGData\Arrhythmia';
            fs = 360;
            segstart_list = 20;%0 : 60 : 20*60;
            diff_h = [2 4 6 8];
            seg_len = 10;
            smoothnessparam_min = -12;
            smoothnessparam_max = 6;
            %             smoothnessparam = logspace(-2.2, 3.5, 250);
            %             smoothnessparam = 0 : 0.1 : 15.0;
        case 'PTB'
            pth = 'J:\ECGData\PTB';
            fs = 1000;
            segstart_list = 20;%[0 20];
            diff_h = [2 4 6 8];
            seg_len = 10;
            smoothnessparam_min = -12;
            smoothnessparam_max = 6;
            %             smoothnessparam = logspace(-2.2, 3.5, 250);
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
    fresults = ['AllDataSmootherResults' type{ttp} '07_NaiveSweepGammaBasedOptimized.txt'];
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
                            %                             noise = NoiseGenerator(noisetypes_list(nt), SignalPower, snr, length(data), fs)';
                            noise =  NoiseGeneratorSimple(noisetypes_list{nt}, NoisePower, fs, length(data));
                            x = data + noise;
                            snr0 = 10*log10(sum(data.^2)/sum((data - x).^2));
                            
                            parameters = [0 logspace(smoothnessparam_min, smoothnessparam_max, INITIAL_PARAM_LIST_NUM)];
                            snrimp = zeros(1, INITIAL_PARAM_LIST_NUM);
                            for pm = 1 : length(parameters),
                                x_filtered_LTI = ECGSmoothnessPriorsDenoiserLTI(x, parameters(pm), mode, diff_h(oo));
                                snr1 = 10*log10(sum(data.^2)/sum((data - x_filtered_LTI).^2));
                                snrimp(pm) = snr1 - snr0;
                            end
                            [YY,II] = sort(snrimp, 2, 'descend');
                            
                            param_down = parameters(II(2));
                            param_up = parameters(II(1));
                            param_mid = (param_up + param_down)/2;
                            
                            x_filtered_LTI_down = ECGSmoothnessPriorsDenoiserLTI(x, param_down, mode, diff_h(oo));
                            x_filtered_LTI_up = ECGSmoothnessPriorsDenoiserLTI(x, param_up, mode, diff_h(oo));
                            
                            snr_down = 10*log10(sum(data.^2)/sum((data - x_filtered_LTI_down).^2));
                            snr_up = 10*log10(sum(data.^2)/sum((data - x_filtered_LTI_up).^2));
                            % % %                             if(snr_down > snr_up)
                            % % %                                 tmp = snr_down;
                            % % %                                 snr_down = snr_up;
                            % % %                                 snr_up = tmp;
                            % % %
                            % % %                                 tmp = param_down;
                            % % %                                 param_down = param_up;
                            % % %                                 param_up = tmp;
                            % % %                             end
                            for pm = 1 : MAX_ITR,
                                x_filtered_LTI_mid = ECGSmoothnessPriorsDenoiserLTI(x, param_mid, mode, diff_h(oo));
                                snr_mid = 10*log10(sum(data.^2)/sum((data - x_filtered_LTI_mid).^2));
                                if(snr_up >= snr_mid && snr_mid >= snr_down)
                                    param_down = param_mid;
                                    snr_down = snr_mid;
                                elseif(snr_down >= snr_mid && snr_mid >= snr_up)
                                    param_up = param_mid;
                                    snr_up = snr_mid;
                                else
                                    break;
                                    % % %                                     tmp = snr_down;
                                    % % %                                     snr_down = snr_up;
                                    % % %                                     snr_up = tmp;
                                    % % %
                                    % % %                                     tmp = param_down;
                                    % % %                                     param_down = param_up;
                                    % % %                                     param_up = tmp;
                                end
                                param_mid = (param_up + param_down)/2;
                            end
                            %                             snr_mid = sum(data.^2)/sum((data - x_filtered_LTI).^2);
                            
                            % % %                             [Y I] = max(snr_mid - snr);
                            strng = [num2str(kk) '\t' fname '\t' num2str(nn) '\t' num2str(dd) '\t' num2str(diff_h(oo)) '\t' num2str(nt)  '\t' num2str(ss) '\t' num2str(snr, '%4.4f') '\t' num2str(snr0, '%4.4f') '\t' num2str(snr_mid, '%2.6g') '\t' num2str(snr_mid - snr0, '%4.4f') '\t' num2str(param_mid, '%d') '\n'];
                            %                         strng = [num2str(kk) '\t' fname '\t' num2str(nn) '\t' num2str(dd)  '\t' num2str(nt)  '\t' num2str(ss) '\t' num2str(snr, '%4.4f') '\t' num2str(snr0(I), '%4.4f') '\t' num2str(snr1(I) - snr0(I), '%4.4f') '\t' num2str(snr2(I) - snr0(I), '%4.4f') '\t' num2str(smoothnessparam(I), '%2.6g') '\t' num2str(I, '%d') '\n'];
                            foutid = fopen(fout, 'a');
                            fprintf(foutid, strng);
                            fclose(foutid);
%                             disp(['SNR = ' num2str(snr)]);
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
