clc
clear all
close all;
addpath('E:\Programs\Baseline Wander Comparison\ECG PTB Database');


f = 1.1;
fs = 1000;
bins = fs/4;

N = 30*fs;
M = 12;

start = 0;
stop = 30;

DataFileNames;

for kk = 1:length(fname),
    InputFileName = fname{kk};
    h = fopen('ECGExtracter.bat','w');
    fprintf(h,'rdsamp -r %s -f %f -t %f >  ECG.txt',InputFileName,start,stop);
    fclose(h);
    !ECGExtracter.bat

    data = load('ECG.txt');
    data = data(:,2:end)';

    bsline = LPFilter(data,.3/fs);                  % baseline wander removal (may be replaced by other approaches)
    x = data - bsline;

    ref = x(1,:);
    peaks = PeakDetection(ref,f/fs);                  % peak detection
    [phase phasepos] = PhaseCalculation(peaks);     % phase calculation

    InputParamFileName = ['dat',num2str(kk,'%3.3d')];
    load(InputParamFileName);

    for i = 1:M;
        [ECGmean,ECGsd,meanphase] = MeanECGExtraction(x(i,:),phase,bins,1); % mean ECG extraction
        params{i}.mean = ECGmean;
        params{i}.std = ECGsd;

        prs = ECG_params{i};
        params{i}.a = prs(1:end/3);
        params{i}.b = prs(end/3+1:2*end/3);
        params{i}.theta = prs(2*end/3+1:end);
        params{i}.label = InputFileName;
    end

    OutputFileName = ['C:\Documents and Settings\a\Desktop\ECG Parameter Database\PTB Database\','params_',num2str(kk+999,'%5.5d')];
    save(OutputFileName,'params');
end

% % % figure;
% % % plot(x');
% % % grid;
% % %
% % % figure;
% % % plot(ECG');
% % % grid;

