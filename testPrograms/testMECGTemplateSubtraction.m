% Sample code for testing naive maternal template subtraction technique
%
% Reza Sameni (C)
% Email: rsameni@shirazu.ac.ir
% Web: www.sameni.info
%
% Crated June 2018

clc;
clear;
close all;
load('SampleFECGDaISy');
data = cell2mat(data);
fs = cell2mat(fs);
L = size(data, 1);
N = size(data, 2);

mref = data(8, :);
fm = 1.0; % approximate maternal heart rate in Hz (BPM/60)
mpeaks = PeakDetection(mref, fm/fs);
mpeaks_logical = logical(mpeaks);

mref_beat_index = 2; % the maternal beat used as reference
before_peak = 31;
after_peak = 31;
I = find(mpeaks == 1);
mref_beat_peak_position = I(mref_beat_index);
mtemplates = data(:, (mref_beat_peak_position - before_peak) : (mref_beat_peak_position + after_peak)); % segment used as template

% replicate the maternal template at the maternal R-peaks
mECG = zeros(size(data));
for k = 1 : length(I)
    ind_start = I(k) - before_peak;
    ins_stop = I(k) + after_peak;
    mECG(:, ind_start : ins_stop) = mtemplates;
end

fECG = data - mECG;

% Diplay Results
t = (0:N-1)/fs;
PlotECG(data, L, 'b', fs);
PlotECG(mECG, L, 'r', fs);
PlotECG(fECG, L, 'g', fs);

figure
hold on;
plot(t, mref);
plot(t(mpeaks_logical), mref(mpeaks_logical), 'ro');
grid
xlabel('time(s)');
ylabel('Amplitude');

ch = 2;
figure;
hold on;
plot(t, data(ch, :), 'b');
plot(t, mECG(ch, :), 'r');
plot(t, fECG(ch, :), 'g');
grid
legend('Noisy', 'mECG', 'Residue');
xlabel('time (s)');
ylabel('Amplitude');
title('Naive Template Extraction');
