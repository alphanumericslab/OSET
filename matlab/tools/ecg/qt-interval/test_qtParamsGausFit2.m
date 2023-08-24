
%% clear
clear
close all
clc

%% load ecg
load ecgptb.mat

%% prepare data
ecg=ecg(10*fs+1:20*fs,:); % 10 sec multichannel ecg
L=size(ecg,1); % ecg length
tst=(1:L)/fs; % time stamp
w1=.75; % window length of the first median filter used in base line removing
w2=.9; % window length of the second median filter used in base line removing
ecg=ecg-(baseline_sliding_window(baseline_sliding_window(ecg', round(w1*fs), 'md'), round(w2*fs), 'mn'))'; % baseline removal
chnl=2;
snr=400;
varNoise=var(ecg(:,chnl))/snr;


%% set the input values
soi.q = [-0.045; -0.015];
soi.t = [.1; .5];
plot_results = true;

%% Bys framework
rng(1); % fix the seed for noise generating
ecgn = ecg + randn(size(ecg)).*sqrt(varNoise); % adding noise to the signal;

[soi_ML, soi_BYS] = QT_analysis_GaussianFit(ecgn, fs, soi, varNoise, plot_results);