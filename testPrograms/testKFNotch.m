%
% Test program for Kalman notch filter
%
% Dependencies: The baseline wander toolbox of the Open Source ECG Toolbox
%
% Open Source ECG Toolbox, version 1.0, October 2007
% Released under the GNU General Public License
% Copyright (C) 2007  Reza Sameni
% Sharif University of Technology, Tehran, Iran -- GIPSA-LAB, INPG, Grenoble, France
% reza.sameni@gmail.com

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.

clc
clear all
close all;

load('SampleECG1.mat'); data = data;
% load('SampleECG2.mat'); data = data(1:15000,6)';
% load('sig'); data = sig(1:15000,1)'; clear sig;

fs = 1000;
f0 = 50;

% data = .1*randn(1,20000);
n = (0:length(data)-1);
x = 100*data + (.055+.02*sin(2*pi*n*.1/fs)).*sin(2*pi*n*f0/fs)/std(data);
% x = data + sin(2*pi*n*f0/fs+pi/9)/std(data);
% % % % % % x = data;% + .2*randn(size(data));


% [y1,y2,Pbar,Phat,PSmoothed,Kgain] = KFNotch(x,f0,fs);
% [y1,y2,Pbar,Phat,PSmoothed,Kgain] = KFNotch(x,f0,fs,1e-3,.1*var(x),.9);
[y1,y2,Pbar,Phat,PSmoothed,Kgain] = KFNotch(x,f0,fs,1e-3,.1*var(x),1);


t = n/fs;

figure;
hold on;
plot(t,data,'b');
plot(t,x,'r');
plot(t,x-y1,'m');
plot(t,x-y2,'g');
plot(t,Kgain(:,:)','k');
grid;
xlabel('time (sec.)');
legend('original ECG','noisy ECG','Kalman filter','Kalman smoother');
% legend('noisy ECG','Kalman filter','Kalman smoother');

figure;
psd(x,length(x)/2,fs);
title('noisy spectrum');

figure;
psd(y1,length(y1)/2,fs);
title('Kalman filter output spectrum');

figure;
psd(y2,length(y2)/2,fs);
title('Kalman smoother output spectrum');

figure;
plot(t,Kgain(1,:)-Kgain(2,:));
grid
