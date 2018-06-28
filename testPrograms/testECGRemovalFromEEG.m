% This program is for filtering the ECG components from EEG signals
% recorded during an fMRI experiment. Here we focus on the AF8 channel
% which was more noisy. The steps of the program may be summarized as
% follows:
% 1- find a rough estimate of the baseline wander and remove it. Different
% baseline wander removal techniques may be used for this step.
% 2- find the R-peaks from the ECG channel and calculate the required 'ECG
% phase' for the KF from these R-peaks.
% 3- find the average of the AF8 channel by synchronous averaging, using
% the R-peaks extracted from the ECG signal. This gives the average ECG
% artifact which has corrupted the AF8 EEG channel.
% 4- fit a Gaussian mixture with a desired number of Gaussian kernels on
% the mean ECG artifact.
% 5- use the Gaussian template for extracting the parameters of the KF
% 6- filter the noisy EEG and estimate the ECG noise
% 7- subtract the ECG estimated by EKF and EKS from the original data.
% 8- the performance of the EKF and EKS may be compared with each other.
% The EKS generally performs better than the EKF.
% 
% You can change or simplify the program if you are interested in filtering
% other signals or you are only interested in the AF8 channel. But you
% always require some ECG signal for synchronizing the R-peaks.
%
%
% Open Source ECG Toolbox, version 1.0, November 2006
% Released under the GNU General Public License
% Copyright (C) 2006  Reza Sameni
% Modified June 2018
%
% Sharif University of Technology, Tehran, Iran -- LIS-INPG, Grenoble, France
% reza.sameni@gmail.com
%
%
% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.

clear;
close all;
load('SampleEEG1.mat');

dat = dat(:,1:25000); % select the desired segment of the data

fs = 250;
t = (0:length(dat)-1)/fs;

z1 = BaseLine2(dat(1,:),fs*.3,fs*.6,'md');      % Cz
dat1 = dat(1,:)-z1;

z2 = BaseLine2(dat(2,:),fs*.3,fs*.6,'md');      % AF8
dat2 = dat(2,:)-z2;

z3 = BaseLine2(dat(3,:),fs*.3,fs*.6,'md');      % ECG
dat3 = dat(3,:)-z3;

peaks = PeakDetection(dat3,1.2/fs);
[phase, phasepos] = PhaseCalculation(peaks);

% a slight phase-shift which helps in bringing all the critical peaks at
% the center of the 0 and 2pi range.
phase1 = PhaseShifting(phase,-pi/2.5);

% [mn1,sd1,Phasemn1] = MeanECGExtraction(dat1,phase1,fs,1);
[mn2,sd2,Phasemn2] = MeanECGExtraction(dat2,phase1,fs,1);
%[mn3,sd3,Phasemn3] = MeanECGExtraction(dat3,phase1,fs,1);

% h1 = ECGBeatFitter(mn1,Phasemn1,'ParamsChannel1','Cz');           % ECG beat fitter GUI
% h2 = ECGBeatFitter(mn2,Phasemn2,'ParamsChannel2','AF8');           % ECG beat fitter GUI
% h3 = ECGBeatFitter(mn3,Phasemn3,'ParamsChannel3','ECG');           % ECG beat fitter GUI
ECGBeatFitter(mn2,sd2,Phasemn2,'ParamsChannel2');               % ECG beat fitter GUI

%//////////////////////////////////////////////////////////////////////////
% KF required parameters
N = length(ParamsChannel2)/3;% number of Gaussian kernels
JJ = find(peaks);
fm = fs./diff(JJ);                  % heart-rate
w = mean(2*pi*fm);                  % average heart-rate in rads.
wsd = std(2*pi*fm,1);               % heart-rate standard deviation in rads.

y = [phase1 ; dat2];                 % observation vector
X0 = [-pi 0]';                                              % initial state vector
P0 = [(2*pi)^2 0 ;0 (10*max(abs(dat2))).^2];                % covariance matrix of the initial state vector
Q = diag( [ (.1*ParamsChannel2(1:N)).^2 (.05*ones(1,N)).^2 (.05*ones(1,N)).^2 (wsd)^2 , (.01*mean(sd2(1:round(end/10))))^2] );  % process noise covariance matrix
R = [(w/fs).^2/12 0 ;0 (mean(sd2(1:round(end/10)))).^2];    % observation noise covariance matrix
Wmean = [ParamsChannel2 w 0]';  % mean of process noise vector
Vmean = [0 0]';                 % mean of observation noise vector
Inits = [ParamsChannel2 w fs];  % parameters for initializing the KF

InovWlen = ceil(.5*fs);         % innovations monitoring window length
tau = [];                       % Kalman filter forgetting time. tau=[] for no forgetting factor
gamma = 1;                      % observation covariance adaptation-rate. 0<gamma<1 and gamma=1 for no adaptation
RadaptWlen = ceil(fs/2);        % window length for observation covariance adaptation

%//////////////////////////////////////////////////////////////////////////
% Kalman filtering and smoothing
[Xekf,Phat,Xeks,PSmoothed,ak] = EKSmoother(y,X0,P0,Q,R,Wmean,Vmean,Inits,InovWlen,tau,gamma,RadaptWlen,1);

%//////////////////////////////////////////////////////////////////////////
% ploting the results
[X,Y,Z] = pol2cart(phase1,1,dat1);
figure;
plot3(X,Y,Z);
grid;
title('the phase representation of Cz');

[X,Y,Z] = pol2cart(phase1,1,dat2);
figure;
plot3(X,Y,Z);
grid;
title('the phase representation of AF8');

[X,Y,Z] = pol2cart(phase1,1,dat3);
figure;
plot3(X,Y,Z);
grid;
title('the phase representation of the ECG');

figure
plot(t,dat(2,:));
hold on;
plot(t,dat(2,:)-Xekf(2,:),'g');
plot(t,dat(2,:)-Xeks(2,:),'r');
grid;
legend('Noisy','Denoised with EKF','Denoised with EKS');
title('EKF & EKS results without baseline wander removal');

figure
plot(t,dat2);
hold on;
plot(t,dat2-Xekf(2,:),'g');
plot(t,dat2-Xeks(2,:),'r');
grid;
legend('Noisy','Denoised with EKF','Denoised with EKS');
title('EKF & EKS results with baseline wander removal');
