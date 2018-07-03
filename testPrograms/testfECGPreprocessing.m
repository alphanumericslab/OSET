% test fetal ECG preprocessing
% Reza Sameni
% July 2018

clc
clear;
close all;

load('FOETAL_ECG.dat'); data = FOETAL_ECG(:,2:end)'; clear FOETAL_ECG; fs = 250;
% load patient165_s0323lre; data = data(:, 2:6)'; fs = 1000;

N = size(data, 2);

% make the data noisy
wgn = randn(size(data));
colored = cumsum( randn(size(data,1), size(data, 2)), 2);
x = data + 5*wgn + 1*colored;

ylp = LPFilter(x, 150/fs);
ybpf = ylp - LPFilter(ylp, 0.01/fs);

ywdn = zeros(size(x));
for k = 1:size(x,1)
%     ywdn(k, :) = wden(x(k, :), 'rigrsure', 's', 'sln', 8,'db4');
    ywdn(k, :) = wden(ybpf(k, :), 'rigrsure', 's', 'mln', 8,'db4');
end

t = (0 : N-1)/fs;

ch = 1;
figure
hold on;
plot(t, x(ch,:));
plot(t, data(ch, :), 'r');
plot(t, ylp(ch, :), 'g');
plot(t, ybpf(ch, :), 'c');
plot(t, ywdn(ch, :), 'm');
grid
legend('noisy', 'original', 'LPF', 'BPF', 'wdn');


figure
psd(data(ch, :), 1024, fs);
title('The original signal spectrum');
% figure
hold on
psd(x(ch, :), 1024, fs);
title('The noisy signal spectrum');

% PlotECG(data, 4, 'b', fs)