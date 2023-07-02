%
% Test program for generating synthetic multichannel ECGs with beat-wise deviations and additive noise
%
% Open Source ECG Toolbox, version 2.0, April 2008
% Released under the GNU General Public License
% Copyright (C) 2008  Reza Sameni
% Sharif University of Technology, Tehran, Iran -- GIPSA-Lab, Grenoble, France
% reza.sameni@gmail.com

%
% The Open Source Electrophysiological Toolbox (OSET), version 3.14, April 2022
% URL: https://github.com/alphanumericslab/OSET
% Copyright (C) 2022  Reza Sameni
% reza.sameni@gmail.com


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
snr = 500;       % signal to noise ratio
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
HR = 45;                    % heart rate in BPM
F = HR/60;                  % heart rate in Hz
k = 1;                      % dipole attenuation parameter
R0 = Rotate3D(0,0,0);       % dipole rotation matrices (tetax,tetay,tetaz)
Lambda = eye(3);

teta0 = -pi/3;              % initial phase of the ECG

%//////////////////////////////////////////////////////////////////////////
% Normal Beat model
tetai.x  = [-1.09  -0.83   -0.19     -.07  0 .06        0.22    1.2 1.42 1.68 2.9];
alphai.x = [0.03   .08    -0.13    .85 1.11 .75     0.06   0.1  0.17 0.39 .03];
bi.x     = [0.0906    0.1057    0.0453    0.0378    0.0332    0.0302    0.0378    0.6040 0.3020  0.1812 .5];

tetai.y  = [-1.1  -0.9 -0.76       -0.11   -.01       0.065  0.8 1.58 2.9];
alphai.y = [0.035 0.015 -0.019     0.32    .51     -0.32    0.04   0.08 .014];
bi.y     = [0.07  .07  0.04        0.055    0.037    0.0604  0.450  0.3 .5];

tetai.z  = [-1.1  -0.93 -0.7      -.4     -0.15    .095    1.05 1.25 1.55 2.8];
alphai.z = [-0.03 -0.14 -0.035    .045     -0.4    .46    -.12 -.2 -.35 -.035];
bi.z     = [.03  .12  .04         .4    .045       .05    .8 .4 .2 .4];

% Abnormal Beat model
% tetai.x  = [-1.09  -0.83   -0.19     -.07  0 .06        0.22    1.2 1.42 1.68 2.9];
% alphai.x = [0.03   .08    -0.13    .85 1.11 .75     0.06   0.1*twafx^2  0.17*twafx^2 0.39*twafx .03];
% bi.x     = [0.0906    0.1057    0.0453    0.0378    0.0332    0.0302    0.0378    0.6040  0.3020*abs(1+twafx/4)   0.1812 .5];
% 
% tetai.y  = [-1.1  -0.9 -0.76       -0.11   -.01       0.065  0.8 1.58 2.9];
% alphai.y = [0.035 0.015 -0.019     0.32    .51     -0.32*twafy^2    0.04*twafy^2   0.08*twafy .014];
% bi.y     = [0.07  .07  0.04        0.055    0.037    0.0604  0.450*abs(1+twafy/4)   0.3 .5];
% 
% tetai.z  = [-1.1  -0.93 -0.7      -.4     -0.15    .095    1.05 1.25 1.55 2.8];
% alphai.z = [-0.03 -0.14 -0.035    .045     -0.4    .46    -.12*twafz^2 -.2*twafz^2 -.35*twafz -.035];
% bi.z     = [.03  .12  .04         .4    .045       .05    .8 .4*abs(1+twafz/4)  .2 .4];

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

amp_deviations = 0.1; % percentage of beat-wise amplitude deviations
center_deviations = 0.1; % percentage of beat-wise wave center deviations
width_deviations = 0.1; % percentage of beat-wise wave width deviations
F_deviations = 0.5;  % percentage of beat-wise HR deviations
[DIP, teta] = DipoleGeneratorStochastic(N, fs, F, F_deviations, alphai, amp_deviations, bi, width_deviations, tetai,center_deviations, teta0);

VCG = R0*Lambda*[DIP.x ; DIP.y ; DIP.z];
s0 = H*VCG;

s = s0 + (sqrt(sum(s0.^2,2))./sqrt(sum(noise.^2,2))/sqrt(10^(snr/10))*ones(1,size(s0,2))).*noise;

ref_ch = 2;
[~, peak_indexes] = PeakDetection(s(ref_ch, :), F/fs);

%//////////////////////////////////////////////////////////////////////////
% data plotting
t = (0 : N-1)/fs;
figure;
hold on
plot(t,1000*s');
plot(t(peak_indexes), 1000*s(ref_ch, peak_indexes), 'bo', 'markersize', 14);
grid
xlabel('time(s)');
ylabel('Amplitude(mV)');
title('Synthetic multi-channel ECG with additive colored noise');

rr_interbvals = diff(peak_indexes);
disp(['HR range = ', num2str(mean(rr_interbvals)), '+-' num2str(std(diff(peak_indexes)))]);