% Sample code for active noise cancellation from ECG signals
%
% Reza Sameni (C)
% Email: rsameni@shirazu.ac.ir
% Web: www.sameni.info
%
% Created May 2007
% Modified June 2018

clc;
clear;
close all;
load('SampleECG128Hz');
fs = 128;
ECG = y(:,2)';
N = length(ECG);
n = (0:N-1);

% Adding noise
rr = .1;
x = ECG + 1.21*sin(2*pi*50*n/fs) + (rr)*randn(1, N);

% Apply the adaptive line enhancer
delay = 3;
taps = 5;
mu = 0.005;
[ECG_estimate, Noise_estimate] = AdaptiveFilter(x, delay, taps, mu);

% Select indexes after transient response for SNR calculation
I = floor(2*length(ECG)/4):ceil(3*length(ECG)/4);

SignalPower = mean(ECG(I).^2);
PreNoisePower = mean((x(I) - ECG(I)).^2);
PostNoisePower = mean((ECG(I) - ECG_estimate(I)).^2);

Initial_SNR = 10*log10(SignalPower / PreNoisePower);
AF_SNR = 10*log10(SignalPower / PostNoisePower);

% Display results
t = n/fs;
figure;
subplot(211);
plot(t, x,'b');
hold on;
plot(t, ECG_estimate,'r');
plot(t, ECG, 'g');
grid;
legend('Noisy', 'ANC Output', 'Original');
xlabel('time (s)');
ylabel('Amplitude');

subplot(212);
plot(t, Noise_estimate);
grid;
legend('Noise Estimate');
xlabel('time (s)');
ylabel('Amplitude');

disp(['Pre Filtering SNR = ' num2str(Initial_SNR)]);
disp(['Post Filtering SNR = ' num2str(AF_SNR)]);
