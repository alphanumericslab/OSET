%
% Test program for synthetic noise generation.
%
% Dependencies: The baseline wander toolbox of the Open Source ECG Toolbox
% and the sample noises
%
% Open Source ECG Toolbox, version 1.0, November 2006
% Released under the GNU General Public License
% Copyright (C) 2006  Reza Sameni
% Sharif University of Technology, Tehran, Iran -- LIS-INPG, Grenoble, France
% reza.sameni@gmail.com

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details. You should have received a copy of the
% GNU General Public License along with this program; if not, write to the
% Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
% MA  02110-1301, USA.

clc
clear
close all;

N = 10000;      % number of samples
NumCh = 8;      % number of channels      
fs = 500;       % sampling rate

% the ratio between the different types of noise
w_bw = 5;       % weight of baseline wander noise in the generated noise (for noisetype = 5)
w_em = 1;       % weight of electrode movement noise in the generated noise (for noisetype = 5)
w_ma = 1;       % weight of muscle artifact noise in the generated noise (for noisetype = 5)

M = 10800;
% original noise template
template =  NoiseGenerator(5,1,0,M,360,[w_bw,w_em,w_ma],1000);

% parameters required for estimating the AR coefficients using a Kalman Filter (KF)
order = 12;                         % AR model order for modeling the ECG noise
[a0,e] = aryule(template,order);    % a global AR model
q = (.05*max(abs(a0(2:end))))^2;    % AR coefficients covariance
R = 1;                              % unit variance noise
p0 = 1e6*q;                         % covariance of the KF initial state
alpha = 1;                          % KF forgetting factor

% time-variant AR parameter estimation using KF and Kalman Smoother (KS)
[Ahat,Asmoothed] = TimeVariantAR(template,order,a0(2:end)',q,R,p0,alpha);

% generating different instances of ECG noise using the time-variant AR parameters
noise1 =  zeros(NumCh,N);
noise2 =  zeros(NumCh,N);
for j = 1:NumCh
    x = randn(1,M);
    y1 =  zeros(M,1);
    y2 =  zeros(M,1);
    for i = order+1:M-1
        y1(i) = (sqrt(1)*x(i)-Ahat(:,i)'*y1(i-1:-1:i-order))/1;         % KF
        y2(i) = (sqrt(1)*x(i)-Asmoothed(:,i)'*y2(i-1:-1:i-order))/1;    % KS
    end
    % resampling the noise matrix to the desired sampling rate 
    n1 = resample((y1-mean(y1))/std(y1),fs,360);
    n2 = resample((y2-mean(y2))/std(y2),fs,360);
    noise1(j,:) = n1(101:N+100);
    noise2(j,:) = n2(101:N+100);
end
template = template(101:N+100);

% plotting the results
time = (0:N-1)/fs;
figure;
hold on;
plot(time,noise1');
plot(time,template,'k','Linewidth',2);
title('template noise vs. synthetic noise using time-variant AR model estimated by Kalman Filter');
xlabel('time(s)');
grid;

figure;
hold on;
plot(time,noise2');
plot(time,template,'k','Linewidth',2);
title('template noise vs. synthetic noise using time-variant AR model estimated by Kalman Smoother');
xlabel('time(s)');
grid;
