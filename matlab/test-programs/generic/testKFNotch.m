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
%
% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.

clc
clear all
close all;

% load('SampleECG1.mat'); data = data;
load('SampleECG2.mat'); data = data(1:15000,6)';
% load('sig'); data = sig(1:15000,1)'; clear sig;

fs = 1000;
f0 = 50;
data = data - LPFilter(data, 1.0/fs);

% data = .1*randn(1,20000);
n = (0:length(data)-1);
% x = 1.0*data + (.055+.02*sin(2*pi*n*.1/fs)).*sin(2*pi*n*f0/fs)/std(data);
x = data + 0.1 * sin(2*pi*n*f0/fs+pi/9)/std(data);
% % % % % % x = data;% + .2*randn(size(data));

% [y1,y2,Pbar,Phat,PSmoothed,Kgain] = KFNotch(x,f0,fs);
% [y1,y2,Pbar,Phat,PSmoothed,Kgain] = KFNotch(x,f0,fs,1e-3,.1*var(x),.9);
[y1,y2,Pbar,Phat,PSmoothed,Kgain] = KFNotch(x, f0, fs, 1e-3, 1*var(x), 1.0);


t = n/fs;

figure;
hold on;
lbl = {};
plot(t,data,'b'); lbl = cat(2, lbl, 'Original ECG');
plot(t,x,'r'); lbl = cat(2, lbl, 'Noisy ECG');
% plot(t,x-y1,'m'); lbl = cat(2, lbl, 'KF powerline estimate');
% plot(t,x-y2,'g'); lbl = cat(2, lbl, 'KS powerline estimate');
plot(t,y1,'m'); lbl = cat(2, lbl, 'Kalman filter');
plot(t,y2,'g'); lbl = cat(2, lbl, 'Kalman smoother');
plot(t,Kgain(:,:)','k'); lbl = cat(2, lbl, 'Kalman gain');
grid;
xlabel('time (s)');
legend(lbl);
title('Powerline noise cancellation with a Kalman filter');

figure;
pwelch(x, length(x)/2, length(x)/8, 1024, fs);
title('Noisy spectrum');

figure;
pwelch(y1,length(y1)/2, length(y1)/8, 1024, fs);
title('Kalman filter output spectrum');

figure;
pwelch(y2,length(y2)/2, length(y2)/8, 1024, fs);
title('Kalman smoother output spectrum');

figure;
plot(t,Kgain(1,:));
plot(t,Kgain(2,:));
grid
title('The Kalman filter gain');
