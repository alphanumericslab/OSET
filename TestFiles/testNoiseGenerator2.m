%
% Test program for noise generation.
%
% Dependencies: The baseline wander toolbox of the Open Source ECG Toolbox
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
clear all
close all;

load('SampleECG1.mat');
%load('SampleECG2.mat'); data = data(1:15000,6)';

fs = 1000;
N = length(data);

% baseline wander removal
bsline = LPFilter(data,.7/fs);          % baseline wander removal (may be replaced by other approaches)
data = data-bsline;
data = data(:);

t = [0:N-1]/fs;

% noise variance calculation
SNR = 5;
SignalPower = mean(data.^2);

beta = 1.5;     % noise color (for noisetype = 1)
w_bw = 1;       % weight of baseline wander noise in the generated noise (for noisetype = 5)
w_em = 1;       % weight of electrode movement noise in the generated noise (for noisetype = 5)
w_ma = 1;       % weight of muscle artifact noise in the generated noise (for noisetype = 5)

noise0 =  NoiseGenerator(0,SignalPower,SNR,N,1);
noise1 =  NoiseGenerator(1,SignalPower,SNR,N,fs,beta,10);
noise2 =  NoiseGenerator(2,SignalPower,SNR,N,fs,0);
noise3 =  NoiseGenerator(3,SignalPower,SNR,N,fs,100);
noise4 =  NoiseGenerator(4,SignalPower,SNR,N,fs,100);
noise5 =  NoiseGenerator(5,SignalPower,SNR,N,fs,[w_bw,w_em,w_ma],1000);

x0 = data + noise0;
x1 = data + noise1;
x2 = data + noise2;
x3 = data + noise3;
x4 = data + noise4;
x5 = data + noise5;

% Plot results
figure;
hold on
plot(t,noise0,'g');
plot(t,x0,'r');
plot(t,data,'b');
grid;
xlabel('time (sec.)');
legend('Noise','Noisy ECG','Original ECG');
title('White Gaussian Noise');

figure;
hold on
plot(t,noise1,'g');
plot(t,x1,'r');
plot(t,data,'b');
grid;
xlabel('time (sec.)');
legend('Noise','Noisy ECG','Original ECG');
title('Colored Gaussian Noise');

figure;
hold on
plot(t,noise2,'g');
plot(t,x2,'r');
plot(t,data,'b');
grid;
xlabel('time (sec.)');
legend('Noise','Noisy ECG','Original ECG');
title('Real Muscle Artifact Noise');

figure;
hold on
plot(t,noise3,'g');
plot(t,x3,'r');
plot(t,data,'b');
grid;
xlabel('time (sec.)');
legend('Noise','Noisy ECG','Original ECG');
title('Real Electrode Movement Noise');

figure;
hold on
plot(t,noise4,'g');
plot(t,x4,'r');
plot(t,data,'b');
grid;
xlabel('time (sec.)');
legend('Noise','Noisy ECG','Original ECG');
title('Real Baseline Wander Noise');

figure;
hold on
plot(t,noise5,'g');
plot(t,x5,'r');
plot(t,data,'b');
grid;
xlabel('time (sec.)');
legend('Noise','Noisy ECG','Original ECG');
title('Mixture or BW, EM, and MA Noise');
