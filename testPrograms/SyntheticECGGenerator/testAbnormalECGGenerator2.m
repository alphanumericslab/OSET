%
% Test program for generating synthetic abnormal multichannel ECGs with additive
% colored noise
% The current example is for the TWA abnormality
%
% Dependencies: The synthetic ECG generator and noise generator package of
%   the Open Source ECG Toolbox
%
% Open Source ECG Toolbox, version 2.0, April 2008
% Released under the GNU General Public License
% Copyright (C) 2008  Reza Sameni
% Sharif University of Technology, Tehran, Iran -- GIPSA-Lab, Grenoble, France
% reza.sameni@gmail.com
%
% Reference:
%   Clifford, Gari D., Shamim Nemati, and Reza Sameni. 
%   "An artificial vector model for generating abnormal electrocardiographic rhythms." Physiological measurement 31.5 (2010): 595.
%

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.

%//////////////////////////////////////////////////////////////////////////
clc
close all;
clear;
% randn('state',2); % For fixing a specific result

%//////////////////////////////////////////////////////////////////////////
% General parameters
N = 10000;       % # of signal samples
fs = 500;       % desired sampling rate

% Noise parameters
snr = 50;       % signal to noise ratio
beta  = 1.5;    % noise color

% Twave parameters
twafx = 0.9; % factor my which to modify T-wave in x-dirn
twafy = 0.9; % factor my which to modify T-wave in y-dirn
twafz = 0.9; % factor my which to modify T-wave in z-dirn

% Heart location
heartlocation = [-25 7 20];     % heart location with respect to the navel (coordinate reference)

% Electrode pair locations
ElecPos = [-10 5 15 ; -10 11 24 ; -10 0 23 ; -5 7 -7 ; -5 -1 -5 ; -10 10 18; -10 0 15; -10 10 15];
ElecNeg = [-10 5 24 ; -35 10 25 ; -10 10 24 ; 0 0 0 ; 0 0 0 ; -35 10 18; -10 10 15; -10 10 24];
NumCh = size(ElecPos,1);

% Dipole parameters
F = .9;                     % heart rate
k = 1;                      % dipole attenuation parameter
R0 = Rotate3D(0,0,0);       % dipole rotation matrices (tetax,tetay,tetaz)
Lambda = eye(3);

teta0 = -pi/3;              % initial phase of the ECG

%//////////////////////////////////////////////////////////////////////////
% Normal Beat model (STATE = 1)
tetai(1).x  = [-1.09  -0.83   -0.19     -.07  0 .06        0.22    1.2 1.42 1.68 2.9];
alphai(1).x = [0.03   .08    -0.13    .85 1.11 .75     0.06   0.1  0.17 0.39 .03];
bi(1).x     = [0.0906    0.1057    0.0453    0.0378    0.0332    0.0302    0.0378    0.6040 0.3020  0.1812 .5];

tetai(1).y  = [-1.1  -0.9 -0.76       -0.11   -.01       0.065  0.8 1.58 2.9];
alphai(1).y = [0.035 0.015 -0.019     0.32    .51     -0.32    0.04   0.08 .014];
bi(1).y     = [0.07  .07  0.04        0.055    0.037    0.0604  0.450  0.3 .5];

tetai(1).z  = [-1.1  -0.93 -0.7      -.4     -0.15    .095    1.05 1.25 1.55 2.8];
alphai(1).z = [-0.03 -0.14 -0.035    .045     -0.4    .46    -.12 -.2 -.35 -.035];
bi(1).z     = [.03  .12  .04         .4    .045       .05    .8 .4 .2 .4];

% Abnormal Beat model (STATE = 2)
tetai(2).x  = [-1.09  -0.83   -0.19     -.07  0 .06        0.22    1.2 1.42 1.68 2.9];
alphai(2).x = [0.03   .08    -0.13    .85 1.11 .75     0.06   0.1*twafx^2  0.17*twafx^2 0.39*twafx .03];
bi(2).x     = [0.0906    0.1057    0.0453    0.0378    0.0332    0.0302    0.0378    0.6040  0.3020*abs(1+twafx/4)   0.1812 .5];

tetai(2).y  = [-1.1  -0.9 -0.76       -0.11   -.01       0.065  0.8 1.58 2.9];
alphai(2).y = [0.035 0.015 -0.019     0.32    .51     -0.32*twafy^2    0.04*twafy^2   0.08*twafy .014];
bi(2).y     = [0.07  .07  0.04        0.055    0.037    0.0604  0.450*abs(1+twafy/4)   0.3 .5];

tetai(2).z  = [-1.1  -0.93 -0.7      -.4     -0.15    .095    1.05 1.25 1.55 2.8];
alphai(2).z = [-0.03 -0.14 -0.035    .045     -0.4    .46    -.12*twafz^2 -.2*twafz^2 -.35*twafz -.035];
bi(2).z     = [.03  .12  .04         .4    .045       .05    .8 .4*abs(1+twafz/4)  .2 .4];

% Note: STM is the State Transition Matrix
%     - For each entry of STM, Sij represents the probability of
%       going from state i to state j in the next beat
%     - Each row of STM should sum up to 1
%     - STM is usually asymmetric

% STM = [0 1 ; 1 0]         % exact alternation: the T-wave will alternate in each beat
STM = [.2 .8 ; .8 .2]       % probabilistic alternation with high probability of state transition
% STM = [.9 .1 ; .99 .01]     % abnormalities with low probability of occurrence

S0 = 1; % Start with a normal beat

%//////////////////////////////////////////////////////////////////////////
% Noise generation
noise =  cumsum(randn(NumCh,N),2);
% % % noise =  zeros(NumCh,N);
% % % for j = 1:NumCh,
% % %     noise(j,:) = NoiseGenerator(1,1,snr,N,fs,beta);
% % % end

%//////////////////////////////////////////////////////////////////////////
% ECG calculation
H = zeros(NumCh, 3);
for i = 1:NumCh
    for j = 1:3
        H(i,j) = k* ((ElecPos(i,j)-heartlocation(j))/sqrt(sum((ElecPos(i,:)-heartlocation).^2))^3 - (ElecNeg(i,j)-heartlocation(j))/sqrt(sum((ElecNeg(i,:)-heartlocation).^2))^3);
    end
end

[DIP, teta] = DipoleGeneratorAbnormal(N,fs,F,alphai,bi,tetai,teta0,STM,S0);

VCG = R0*Lambda*[DIP.x ; DIP.y ; DIP.z];
s0 = H*VCG;

s = s0 + (sqrt(sum(s0.^2,2))./sqrt(sum(noise.^2,2))/sqrt(10^(snr/10))*ones(1,size(s0,2))).*noise;

%//////////////////////////////////////////////////////////////////////////
% data plotting
t = (0 : N-1)/fs;
figure;
plot(t,1000*s');
grid
xlabel('time(s)');
ylabel('Amplitude(mV)');
title('Synthetic multi-channel ECG with additive colored noise');