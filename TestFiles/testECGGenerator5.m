%
% Test program for generating synthetic multichannel maternal/fetal ECGs
% mixtures with additive colored noise
%
% Dependencies: The synthetic ECG generator and noise generator package of
%   the Open Source ECG Toolbox
%
% Open Source ECG Toolbox, version 2.0, December 2008
% Released under the GNU General Public License
% Copyright (C) 2008  Reza Sameni
% Shiraz University, Shiraz, Iran
% reza.sameni@gmail.com

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.

clc
clear;
close all;

%//////////////////////////////////////////////////////////////////////////
% Signal/Noise parameters

fs = 1000;   % sampling rate
dim = 32;    % number of channels

SIR = -15;  % fetal signal power to maternal interference power ratio in dB
SNR = -10;  % fetal signal power to noise power ratio in dB
beta = 2;   % noise color (for colored noise only)


w_bw = 1;       % weight of baseline wander noise in the generated noise (for noisetype = 5)
w_em = 1;       % weight of electrode movement noise in the generated noise (for noisetype = 5)
w_ma = 1;       % weight of muscle artifact noise in the generated noise (for noisetype = 5)

%//////////////////////////////////////////////////////////////////////////
% Load maternal dipole parameters from ECG parameter database
BPM_m = [84*ones(1,50), 84:120, 120*ones(1,20) 120:-1:72, 72*ones(1,20)]; % fetal heart rate vector
BPM_m = BPM_m + 1*randn(size(BPM_m));
teta0_m = pi/3.1;
load('params_02000.mat');
params_m = params;

alphai_m.x = params_m{1}.a;
alphai_m.y = params_m{2}.a;
alphai_m.z = params_m{3}.a;

bi_m.x = params_m{1}.b;
bi_m.y = params_m{2}.b;
bi_m.z = params_m{3}.b;

tetai_m.x = params_m{1}.theta;
tetai_m.y = params_m{2}.theta;
tetai_m.z = params_m{3}.theta;

[DIP_m teta_m] = DipoleGenerator4(BPM_m,fs,alphai_m,bi_m,tetai_m,teta0_m);


%//////////////////////////////////////////////////////////////////////////
% Load fetal dipole parameters from ECG parameter database
BPM_f = [120*ones(1,50), 120:140, 140:-1:100, 100*ones(1,100)]; % fetal heart rate vector
BPM_f = BPM_f + 1*randn(size(BPM_f));
teta0_f = pi/3.1;
load('params_02001.mat');
params_f = params;

alphai_f.x = params_f{1}.a;
alphai_f.y = params_f{2}.a;
alphai_f.z = params_f{3}.a;

bi_f.x = params_f{1}.b;
bi_f.y = params_f{2}.b;
bi_f.z = params_f{3}.b;

tetai_f.x = params_f{1}.theta;
tetai_f.y = params_f{2}.theta;
tetai_f.z = params_f{3}.theta;

[DIP_f teta_f] = DipoleGenerator4(BPM_f,fs,alphai_f,bi_f,tetai_f,teta0_f);


%//////////////////////////////////////////////////////////////////////////
% making the raw data
Nm = size(DIP_m.x,2);
Nf = size(DIP_f.x,2);

%//////////////////////////////////////////////////////////////////////////
% Noise generation

N = min([Nm,Nf]);
ns = zeros(dim,N);
for i = 1:dim,
    ns(i,:) = NoiseGenerator(1,1,0,N,fs,beta,mod(i,2));                     % colored noise with 2 independent dimensions
%     ns(i,:) = NoiseGenerator(0,1,0,N,mod(i,3));                           % white noise with 3 independent dimensions
%      ns(i,:) =  NoiseGenerator(5,1,0,N,fs,[w_bw,w_em,w_ma],mod(i,8)*1000); % real noise with 8 independent dimensions
end
noise = diag(rand(dim,1))*ns;   % randomizing the noise power in each channel


%//////////////////////////////////////////////////////////////////////////
% Generating a random mixing matrix with a desired degree between its subspaces

% [A,B,ang] = RandomMatrices(5,dim,0);    % maximum subspace angles smaller or equal to 5 degrees 
[A,B,ang] = RandomMatrices(40,dim,1); % minimum subspace angles greater or equal to 40 degrees

if isempty(A)
    error('Random matrix not found....');
end

%//////////////////////////////////////////////////////////////////////////
% mixing up the generated signals and noises with the desired SIRs and SNRs

xm = A*[DIP_m.x ; DIP_m.y ; DIP_m.z];
xf = B*[DIP_f.x ; DIP_f.y ; DIP_f.z];
ref = [DIP_f.x ; DIP_f.y ; DIP_f.z];    % can be used to rank the signals extracted by signal processing methods

xmpower = sum(sum(xm.^2));
xfpower = sum(sum(xf.^2));
npower = sum(sum(noise.^2));

alpha = sqrt(xfpower/xmpower/10^(SIR/10));
beta = sqrt(xfpower/npower/10^(SNR/10));

x = alpha*xm(:,1:N) + xf(:,1:N) + beta*noise(:,1:N);

%//////////////////////////////////////////////////////////////////////////
% a sample processing of the data

% baseline wander (BW) removal
y = x - LPFilter(x,.7/fs);

peaks_m = PeakDetection(xm(1,1:N),mean(BPM_m)/60/fs);
peaks_f = PeakDetection(xf(1,1:N),mean(BPM_m)/60/fs);
Im = find(peaks_m);
If = find(peaks_f);
mat_hr = fs*60./diff(Im);
fet_hr = fs*60./diff(If);

t = (0:N-1)/fs;

% % % % jade before BW removal
% % % W = jadeR(x);
% % % s1 = W*x;
% % % 
% jade after BW removal
% % % W = jadeR(y);
% % % s2 = W*y;

PlotECG(x,4,'b',fs);

figure
plot(mat_hr);
hold on
plot(fet_hr,'r');
grid
figure
% subplot(2,1,1);
plot(t,x(1,:));
hold on;
plot(t(Im),x(1,Im),'ro');
plot(t(If),x(1,If),'go');
grid
axis tight
% % % PlotECG(s1,8,'r',fs);
% % % PlotECG(s2,8,'m',fs);


figure
subplot(2,1,1);
plot(BPM_m,'linewidth',2);
xlabel('beat','fontsize',12);
ylabel('BPM','fontsize',12);
grid
axis tight
box on
legend('maternal');
set(gca,'fontsize',12)

subplot(2,1,2);
plot(BPM_f,'r','linewidth',2);
xlabel('beat','fontsize',12);
ylabel('BPM','fontsize',12);
grid
axis tight
box on
legend('fetal');
set(gca,'fontsize',12)
