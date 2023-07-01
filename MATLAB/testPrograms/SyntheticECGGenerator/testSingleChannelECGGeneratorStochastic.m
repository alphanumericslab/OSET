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
% Dipole parameters
HR = 60;                    % heart rate in BPM
F = HR/60;                  % heart rate in Hz
k = 1;                      % dipole attenuation parameter
R0 = Rotate3D(0,0,0);       % dipole rotation matrices (tetax,tetay,tetaz)
Lambda = eye(3);
teta0 = -pi/3;              % initial phase of the ECG
n_bins = 250;
amp_deviations = 0.1; % percentage of beat-wise amplitude deviations
center_deviations = 0.1; % percentage of beat-wise wave center deviations
width_deviations = 0.1; % percentage of beat-wise wave width deviations
F_deviations = 0.1;  % percentage of beat-wise HR deviations

%//////////////////////////////////////////////////////////////////////////
% Sample beat types
tetai  = [-1.09  -0.83   -0.19     -.07  0 .06        0.22    1.2 1.42 1.68 2.9];
alphai = [0.03   .08    -0.13    .85 1.11 .75     0.06   0.1  0.17 0.39 .03];
bi     = [0.0906    0.1057    0.0453    0.0378    0.0332    0.0302    0.0378    0.6040 0.3020  0.1812 .5];

% tetai  = [-1.1  -0.9 -0.76       -0.11   -.01       0.065  0.8 1.58 2.9];
% alphai = [0.035 0.015 -0.019     0.32    .51     -0.32    0.04   0.08 .014];
% bi     = [0.07  .07  0.04        0.055    0.037    0.0604  0.450  0.3 .5];
%
% tetai  = [-1.1  -0.93 -0.7      -.4     -0.15    .095    1.05 1.25 1.55 2.8];
% alphai = [-0.03 -0.14 -0.035    .045     -0.4    .46    -.12 -.2 -.35 -.035];
% bi     = [.03  .12  .04         .4    .045       .05    .8 .4 .2 .4];
%
% alphai = [.1 -.2 1 -.3 .15];
% bi = [.3 .2 .2 .2 .5];
% tetai = pi*[-90 -10 0 15 100]/180;

%//////////////////////////////////////////////////////////////////////////
% Noise generation
noise =  cumsum(randn(1, N),2);
% % %     noise = NoiseGenerator(1,1,snr,N,fs,beta);

%//////////////////////////////////////////////////////////////////////////
% ECG calculation
[ECG, teta] = SingleChannelECGGeneratorStochastic(N, fs, F, F_deviations, alphai, amp_deviations, bi, width_deviations, tetai,center_deviations, teta0);

s = ECG + (sqrt(sum(ECG.^2,2))./sqrt(sum(noise.^2,2))/sqrt(10^(snr/10))*ones(1,size(ECG,2))).*noise;

[peaks, peak_indexes] = PeakDetection(s, F/fs);
[phase, ~] = PhaseCalculation(peaks); % phase calculation

%//////////////////////////////////////////////////////////////////////////
% find the beat covariances
phaseshift = pi/n_bins;
pphase = PhaseShifting(phase, phaseshift); % phase shifting to compensate half a bin shift

[M, end_of_beat_indexes] = ECGPhaseToMatrix(pphase, n_bins); % Convert ECG phase into a (phase x time) matrix

x_stacked = zeros(length(end_of_beat_indexes), n_bins);
for ll = 1 : length(end_of_beat_indexes)
    if ll > 1
        beat_start = end_of_beat_indexes(ll-1) + 1;
    else
        beat_start = 1;
    end
    beat_end = end_of_beat_indexes(ll);
    x_stacked(ll, :) = ECG(beat_start : beat_end) * (diag(max(1, sum(M(:, beat_start : beat_end), 2))) \ M(:, beat_start : beat_end))';
end
ECG_mean = mean(x_stacked, 1);
ECG_median = median(x_stacked, 1);

[~, ~, meanphase, ~, ECGSamplesPerBin] = MeanECGExtraction(ECG, phase, n_bins, 1); % mean ECG extraction

%     x_stacked_zero_mean = x_stacked - ECG_avg(ones(1, size(x_stacked, 1)), :);
%     K_x_ph = x_stacked_zero_mean' * x_stacked_zero_mean / size(x_stacked_zero_mean, 1);
K_x_ph = cov(x_stacked); % identical

%//////////////////////////////////////////////////////////////////////////
% data plotting
t = (0 : N-1)/fs;
figure;
hold on
plot(t,1000*s');
plot(t(peak_indexes), 1000*s(peak_indexes), 'bo', 'markersize', 14);
grid
xlabel('time(s)');
ylabel('Amplitude(mV)');
title('Synthetic multi-channel ECG with additive colored noise');

rr_interbvals = diff(peak_indexes);
disp(['HR range = ', num2str(mean(rr_interbvals)), '+-' num2str(std(diff(peak_indexes)))]);

figure
mesh(K_x_ph)