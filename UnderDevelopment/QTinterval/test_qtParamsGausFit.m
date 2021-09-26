
%% clear
clear all
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
ecg=ecg-(BaseLine1(BaseLine1(ecg', round(w1*fs), 'md'), round(w2*fs), 'mn'))'; % baseline removal
chnl=2;


%% set the input values
soi.q=[-0.045; -0.015];
soi.t=[.1; .5];


%% ml framework
[mlGaussParams, rPeaks, soi,  waveParams, qtInt]=qtParamsGausFit(ecg, fs, [],[],soi);

% polt the evaluated gaussians on the signal
figure; p1=plot(tst,ecg(:,chnl)); % plot the channel
hold on
tq=soi.q(1):1/fs:soi.q(2);
tt=soi.t(1):1/fs:soi.t(2);
for j=1:length(rPeaks)
    p2=plot(tq+rPeaks(j)/fs,GausVal(tq,mlGaussParams.q(:,j,chnl)),'r-');
    p3=plot(tt+rPeaks(j)/fs,GausVal(tt,mlGaussParams.t(:,j,chnl)),'r-');
    p4=plot([waveParams.q(1,chnl) waveParams.t(2,chnl)]+rPeaks(j)/fs, ...
        ecg(floor(fs*[waveParams.q(1,chnl) waveParams.t(2,chnl)])+rPeaks(j),chnl),'c*');
end
legend([p1 p2 p4],'ecg', 'Gaussians', 'q/t onset/offset')
title 'ML framework'


%% Bys framework
snr=400;
varNoise=var(ecg(:,chnl))/snr;
rng(1); % fix the seed for noise generating
ecgn=ecg + randn(size(ecg)).*sqrt(varNoise); % adding noise to the signal;

PrMu.q=mean(mlGaussParams.q(:,:,chnl),2);
PrCov.q=cov(mlGaussParams.q(:,:,chnl)');

PrMu.t=mean(mlGaussParams.t(:,:,chnl),2);
PrCov.t=cov(mlGaussParams.t(:,:,chnl)');

[bysGaussParams, rPeaks, soi, waveParams, qtInt]=qtParamsGausFit(ecgn(:,chnl), fs, [], [], soi, [], [], [], PrMu, PrCov, varNoise);

% polt the evaluated gaussians on the signal
figure; p1=plot(tst,ecgn(:,chnl)); % plot the channel
hold on
tq=soi.q(1):1/fs:soi.q(2);
tt=soi.t(1):1/fs:soi.t(2);
for j=1:length(rPeaks)
    p2=plot(tq+rPeaks(j)/fs,GausVal(tq,bysGaussParams.q(:,j)),'r-');
    p3=plot(tt+rPeaks(j)/fs,GausVal(tt,bysGaussParams.t(:,j)),'r-');
    p4=plot([waveParams.q(1,chnl) waveParams.t(2,chnl)]+rPeaks(j)/fs, ...
        ecgn(floor(fs*[waveParams.q(1,chnl) waveParams.t(2,chnl)])+rPeaks(j),chnl),'c*');
end
legend([p1 p2 p4],'ecg', 'Gaussians', 'q/t onset/offset')
title 'BYS framework'

