
%% clear
clear 
% close all
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


[~, rPeaks] = peak_det_amp_threshold(ecg(:,chnl),70/60/fs,.6, 2, 'MEDIAN');

%% set the input values
soi1.q=[-0.080; -0.020];
soi1.t=[.1; .5];

[~, Q,~, ~, T] = ecgWavesSoI_RR(ecg(:,chnl), rPeaks, 1);
soi2.q = Q(:,[1,3])/fs;
soi2.t = T(:,[1,3])/fs;

[~, Q,~, ~, T] = ecgWavesSoI_k(ecg(:,chnl), rPeaks, 1);
soi3.q = Q(:,[1,3])/fs;
soi3.t = T(:,[1,3])/fs;
%% initial params
% p0.q=[0; .02;  -.05];
% p0.t=[0; .05; .25];
p0=[];

%% bounds
lb=[]; ub=[];

%% beta
beta = 3;

%% optimization setting
options = struct('SpecifyObjectiveGradient',true);


%% qt interval estimation
qtInt1=qtIntGausFit_v1(ecg(:,chnl), fs, rPeaks, 3, 3);
qtInt2=qtIntGausFit_v2(ecg(:,chnl), fs, rPeaks, 3, 3);
qtInt3=qtIntGausFit_v3(ecg(:,chnl), fs, rPeaks, 3, 3);

%% ml framework
soi=soi1;
[mlGaussParams, rPeaks, soi,  waveParams, qtInt]=qtParamsGausFit_cl(ecg(:,chnl), fs, rPeaks, p0, soi, lb, ub, beta, options);

% polt the evaluated gaussians on the signal
figure; p1=plot(tst,ecg(:,chnl)); % plot the channel
hold on

for j=1:length(rPeaks)
    tq=soi.q(j,1):1/fs:soi.q(j,2);
    tt=soi.t(j,1):1/fs:soi.t(j,2);
    p2=plot(tq+rPeaks(j)/fs,GausVal(tq,mlGaussParams.q(:,j)),'r-');
    p3=plot(tt+rPeaks(j)/fs,GausVal(tt,mlGaussParams.t(:,j)),'r-');
    p4=plot([waveParams.q(1,j) waveParams.t(2,j)]+rPeaks(j)/fs, ...
        ecg(floor(fs*[waveParams.q(1,j) waveParams.t(2,j)])+rPeaks(j),chnl),'c*');
end
legend([p1 p2 p4],'ecg', 'Gaussians', 'q/t onset/offset')
title 'ML framework'


%% Bys framework
snr=400;
varNoise=var(ecg(:,chnl))/snr;
rng(1); % fix the seed for noise generating
ecgn=ecg + randn(size(ecg)).*sqrt(varNoise); % adding noise to the signal;

PrMu.q=mean(mlGaussParams.q(:,:),2);
PrCov.q=cov(mlGaussParams.q(:,:)');

PrMu.t=mean(mlGaussParams.t(:,:),2);
PrCov.t=cov(mlGaussParams.t(:,:)');

[bysGaussParams, rPeaks, soi, waveParams, qtInt]=qtParamsGausFit_cl(ecgn(:,chnl), fs, rPeaks, p0, soi, lb, ub, beta, PrMu, PrCov, varNoise);

% polt the evaluated gaussians on the signal
figure; p1=plot(tst,ecgn(:,chnl)); % plot the channel
hold on

for j=1:length(rPeaks)
    tq=soi.q(j,1):1/fs:soi.q(j,2);
    tt=soi.t(j,1):1/fs:soi.t(j,2);
    p2=plot(tq+rPeaks(j)/fs,GausVal(tq,bysGaussParams.q(:,j)),'r-');
    p3=plot(tt+rPeaks(j)/fs,GausVal(tt,bysGaussParams.t(:,j)),'r-');
    p4=plot([waveParams.q(1,j) waveParams.t(2,j)]+rPeaks(j)/fs, ...
    ecg(floor(fs*[waveParams.q(1,j) waveParams.t(2,j)])+rPeaks(j),chnl),'c*');
end
legend([p1 p2 p4],'ecg', 'Gaussians', 'q/t onset/offset')
title 'BYS framework'

