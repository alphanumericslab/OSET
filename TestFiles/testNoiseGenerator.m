%
% Test program for colored noise generation.
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

% load('SampleECG1.mat');
load('SampleECG2.mat'); data = data(1:15000,6)';

fs = 1000;
N = length(data);

% baseline wander removal
bsline = LPFilter(data,.7/fs);          % baseline wander removal (may be replaced by other approaches)
%bsline = BaseLineKF(data,.5/fs);       % baseline wander removal (may be replaced by other approaches)
data = data-bsline;

t = [0:N-1]/fs;

% noise variance calculation
SNR = 5;
SignalPower = mean(data.^2);
NoisePower = SignalPower / 10^(SNR/10);

noisetype = 5;

beta = 1.5;     % noise color (for noisetype = 1)
w_bw = 1;       % weight of baseline wander noise in the generated noise (for noisetype = 5)
w_em = 1;       % weight of electrode movement noise in the generated noise (for noisetype = 5)
w_ma = 1;       % weight of muscle artifact noise in the generated noise (for noisetype = 5)

switch noisetype
    case 0     % white noise
        noise = sqrt(NoisePower)*randn(size(data));
        
    case 1     % colored noise
        noise = ColoredNoise(sqrt(NoisePower),N,fs,beta);
        
    case 2     % real muscle artifacts
        load('MA.mat');artifact = MA(:,2);
        artifact = resample(artifact,fs,360);
        artifact = artifact(1:N)';
        noise = sqrt(NoisePower)*(artifact - mean(artifact))/std(artifact,1);
        
    case 3     % real electrode movements
        load('EM.mat');artifact = EM(:,3);
        artifact = resample(artifact,fs,360);
        artifact = artifact(1:N)';
        noise = sqrt(NoisePower)*(artifact - mean(artifact))/std(artifact,1);
        
    case 4     % real baseline wander
        load('BW.mat');artifact = BW(:,3);
        artifact = resample(artifact,fs,360);
        artifact = artifact(1:N)';
        noise = sqrt(NoisePower)*(artifact - mean(artifact))/std(artifact,1);
    
    case 5     % mixture of real baseline wander, electrode movements, muscle artifacts
        load('BW.mat'); bw = BW(:,3);    bw = (bw-mean(bw))/std(bw);
        load('EM.mat'); em = EM(:,3);    em = (em-mean(em))/std(em);
        load('MA.mat'); ma = MA(:,3);    ma = (ma-mean(ma))/std(ma);
        artifact = (w_bw*bw + w_em*em + w_ma*ma)/(w_bw + w_em + w_ma);
        artifact = resample(artifact,fs,360);
        artifact = artifact(1:N)';
        noise = sqrt(NoisePower)*(artifact - mean(artifact))/std(artifact,1);
end

x = data + noise;

figure;
plot(t,noise);
grid;

figure;
plot(t,x,'r');
hold on
plot(t,data,'b');
grid;
xlabel('time (sec.)');
legend('Noisy ECG','Original ECG');
