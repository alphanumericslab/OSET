% A test script for Laplacian Eigen-maps for ECG signals
% Reza Sameni, 2020
% The The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET.git

clc
clear
close all;

load('SampleECG2.mat'); data = data(:,2:end)';

%//////////////////////////////////////////////////////////////////////////
% algorithm parameters
fs = 1000.0;
f = 1.2;
f_cut = 3.0;
w_before_R = 0.3;
w_after_R = 0.4;
refch = 1;
detection_method = 'localpeaks';
correct_peaks = 0;
N = size(data, 1); % the number of channels
T = size(data, 2); % the number of samples
wlenL = round(w_before_R*fs); % before R-peak
wlenR = round(w_after_R*fs); % after R-peak %round(min(D)/2);
fmax = 2*f; % maximum expected frequency of the R-peaks
QRS_detector_template_width = 120e-3; % the matched filter window width used for R-peak detection
QRS_detector_template_std = 10e-3; % the matched filter standard deviation used for R-peak detection
beat_rejection_th = 0.5; % beat similarity rejection threshold
R_peak_rejection_th = 0.1; % R-peak matched filter rejection threshold
%//////////////////////////////////////////////////////////////////////////
% the denoising algorithm

% baseline wander removal
bsline = LPFilter(data, f_cut/fs);                  % baseline wander removal (may be replaced by other approaches)
x = data - bsline;

% R-peak detection
if(strcmp(detection_method,'matched'))
    tt = -QRS_detector_template_width/2:1/fs:QRS_detector_template_width/2;
    QRS_template = exp(-tt.^2/(2*QRS_detector_template_std^2));
    peaks = PeakDetection3(x(refch, :), fs, QRS_template, R_peak_rejection_th, fmax);
elseif(strcmp(detection_method, 'localpeaks'))
    peaks = PeakDetection(x(refch, :),f/fs);                  % peak detection
else
    error('unknown R-peak detection method')
end
peak_indexes = find(peaks);

segment_width = wlenL + wlenR + 1;
% refine the initial detected beats
if(correct_peaks)
    Y = zeros(length(peak_indexes), segment_width);
    for k = 1 : length(peak_indexes)
        start = max(peak_indexes(k)-wlenL, 1);
        stop = min(peak_indexes(k)+wlenR, T);
        xx = [zeros(1, start - (peak_indexes(k)-wlenL)), x(refch, start : stop), zeros(1, peak_indexes(k)+wlenR - T)];
        Y(k, :) = xx - x(refch, peak_indexes(k));
    end
    MN = RWAverage(Y);
    err = Y - MN(ones(length(peak_indexes),1),:);
    beat_error_vs_avg = std(err, [], 2)./std(Y, [], 2);
    Omit = beat_error_vs_avg > beat_rejection_th;
    peak_indexes(Omit) = [];
    %     peaks = zeros(1, T);
    %     peaks(peak_indexes) = 1; % update the peaks with the refined beats
end

if(isempty(peak_indexes))
    error('No R-peaks available. Try calling the function with correct_peaks = 0');
end

% extract the ECG windows arround each peak
MN_all = zeros(N, segment_width);
X = zeros(N, segment_width, length(peak_indexes));
C = zeros(N, N, length(peak_indexes));
t_segment = (0 : segment_width-1)/fs;
for k = 1 : length(peak_indexes)
    for ch = 1 : N
        start = max(peak_indexes(k)-wlenL, 1);
        stop = min(peak_indexes(k)+wlenR, T);
        xx = [zeros(1, start - (peak_indexes(k)-wlenL)), x(ch, start : stop), zeros(1, peak_indexes(k)+wlenR - T)];
        X(ch, :, k) = xx - x(ch, peak_indexes(k));
    end
    C(:, :, k) = X(:, :, k)*X(:, :, k)';
end

K = 7;
L = size(C, 3);
d = nan(K, L);
V = nan(K, K, L);
kappa = 1;
for n = K : L
    [V(:, :, n), d(:, n)] = LaplacianEigenmap(C(:, :, n-K+1:n), kappa);
end

wlen = 1000;
nrm = filtfilt(ones(1, wlen), wlen, sqrt(mean(x.^2, 1)));

PlotECG(x, ceil(N/3), 'b', fs, 'Original signals');

t = (0:T-1)/fs;
figure;
plot(t, x(refch, :));
hold on;
plot(t(peak_indexes), x(refch, peak_indexes), 'ro');
grid
xlabel('time(s)');
ylabel('Amplitude');
title('ECG R-peaks');

figure
plot(t(peak_indexes), log(d(2:5,: )'));
hold on
plot(t, log(nrm/max(nrm)), 'k');
grid
xlabel('time(s)');
ylabel('Laplacian eigenmaps');
title('ECG Spectral Embedding using Laplacian Eigenmaps');

% V1 = squeeze(V(2, :, :));
% V2 = squeeze(V(3, :, :));
% V3 = squeeze(V(4, :, :));
V1 = squeeze(V(:, 2, :));
V2 = squeeze(V(:, 3, :));
V3 = squeeze(V(:, 4, :));

figure
plot3(V1(1, :), V1(2, :), V1(3, :), 'b.');
hold on
plot3(V2(1, :), V2(2, :), V2(3, :), 'r.');
plot3(V3(1, :), V3(2, :), V3(3, :), 'g.');
grid
