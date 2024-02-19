clear;
close all;
clc;

type = 'Normal';

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
wlen = 150e-3; % window length (s)
trial_num = 1;
mode = 3;
ForgettingFactor = 0.9; % forgetting rate of the smoothness factor calculated in each block
NoisePowerOverEstimationFactor = 1.0;
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
    for nn = 1 : length(segstart_list)
        fname = lst(kk).name;
        fname = fname(1:end-4);
        segstart = segstart_list(nn);
        segstop = segstart_list(nn) + seg_len;
        system(['rdsamp -r ' fname ' -f ' num2str(segstart) ' -t ' num2str(segstop) ' > sampledata.txt']);
        %         system(['rdsamp -r ' fname ' > sampledata.txt']);
        data_all_channels = load([pth '\sampledata.txt'])';
        data_all_channels = data_all_channels(2:end,:);
        %%%%%% DO WHAT YOU WANT HERE
    end
end
cd(wd);

% system(['rdsamp -r ' fname ' -f ' num2str(segstart) ' -t ' num2str(segstop) ' -s ' channels ' > data.txt']);
