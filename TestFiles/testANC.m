% Sample code for active noise cancellation from audio signals
%
% Reza Sameni (C)
% Email: rsameni@shirazu.ac.ir
% Web: www.sameni.info
%
% Created October 2004
% Modified June 2018

clear;
close all;
clc

load SampleSpeech10kHz
fs = 10000; % Sampling rate of speeches (s1 and s2) is 10kHz.

s1 = s1(:).';   % row vector
s2 = s2(:).';   % row vector

s1 = s1 / max(abs(s1)); % normalization
s2 = s2 / max(abs(s2)); % normalization

A1 = 0.9; % Noise level
A2 = 0.8; % Noise level
f1 = 50; % Noise frequency 1 (Hz)
f2 = 95; % Noise frequency 2 (Hz)
mu1 = 0.001; % Adaptation factor
mu2 = 0.001; % Adaptation factor

n = 0:(length(s1)-1);

rand_phase1 = rand*2*pi;
noise1 = A1*cos(2*pi*f1*n/fs + rand_phase1);

rand_phase2 = rand*2*pi;
noise2 = A2*cos(2*pi*f2*n/fs + rand_phase2);

x1 = s1 + noise1; % Noisy speech signal f1 hum
x2 = s2 + noise2; % Noisy speech signal f2 hum

out1 = zeros(size(s1)); % Output vector 1
out2 = zeros(size(s2)); % Output vector 2

%   Noise cancelling channel 1
r1 = cos(2*pi*f1*n/fs);
r2 = sin(2*pi*f1*n/fs);
w = [1; 1];
for k = 1:length(s1)
    noise_estimated = w(1)*r1(k) + w(2)*r2(k);
    out1(k) = x1(k) - noise_estimated;
    w = w + mu1*out1(k)*[r1(k); r2(k)];
end

%   Noise cancelling channel 2
r1 = cos(2*pi*f2*n/fs);
r2 = sin(2*pi*f2*n/fs);
w = [1; 0];
for k = 1:length(s2)
    noise_estimated = w(1)*r1(k) + w(2)*r2(k);
    out2(k) = x2(k) - noise_estimated;
    w = w + mu2*out2(k)*[r1(k); r2(k)];
end

% Display results
t = n/fs;
figure
hold on;
plot(t, x1);
plot(t, out1, 'r');
plot(t, s1, 'g');
grid
legend('Noisy', 'ANC Output', 'Original');
xlabel('time (s)');
ylabel('Amplitude');
title('Channel 1');

figure
hold on;
plot(t, x2);
plot(t, out2, 'r');
plot(t, s2, 'g');
grid
legend('Noisy', 'ANC Output', 'Original');
xlabel('time (s)');
ylabel('Amplitude');
title('Channel 2');
