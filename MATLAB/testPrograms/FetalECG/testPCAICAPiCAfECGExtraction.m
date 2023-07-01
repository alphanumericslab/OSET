% Comparing different multichannel fetal ECG extraction methods for test
% Reza Sameni, copyright 2018
%

clear;
close all;

load('FOETAL_ECG.dat'); data = FOETAL_ECG(:,2:end)'; time = FOETAL_ECG(:,1)'; clear FOETAL_ECG; fs = 250;

L1 = size(data,1);
L2 = size(data,2);

% BW removal
bl = LPFilter(data,1.5/fs);
x = data - bl;

% PCA decorrelation (in descending eigenvalue order)
Cx = cov(x');
Cx = (Cx + Cx')/2;
[V, D] = eig(Cx);
D = diag(D);
[B, I_pca] = sort(D, 'descend');
lambdas = D(I_pca);
s_pca = V(:, I_pca)'*x;

% JADE
[A_jade, s_jade] = jade(x);

% fastICA
s_fastica = fastica(x, 'verbose', 'off', 'displayMode', 'off');

% SOBI
[W_sobi, s_sobi] = sobi(x, [], 1000);

% PiCA based on maternal peaks
mref = x(6,:);
f_mat = 1.4; % the expected maternal heart rate in Hz
mpeaks = PeakDetection(mref, f_mat/fs); % maternal R-peak detection
Jm = find(mpeaks);
s_pica_m = PiCA(x, mpeaks);

% PiCA based on fetal peaks
fref = s_jade(5, :); % CAUTION: CHANNEL SELECTED BY VISUAL INSPECTION
f_fet = 2.2; % the expected fetal heart rate in Hz
fpeaks = PeakDetection(fref, f_fet/fs); % fetal R-peak detection
Jf = find(fpeaks);
s_pica_f = PiCA(x, fpeaks);

% PiCA based on maternal-fetal peaks
s_pica_mf = PiCA(x, mpeaks, fpeaks);

% ploting the results
figure
plot(time, mref);
hold on
plot(time(Jm), mref(Jm), 'ro');
grid
title('maternal R-peak detection');

figure
plot(time, fref);
hold on
plot(time(Jf), fref(Jf), 'ro');
grid
title('fetal R-peak detection');

PlotECG(x, 8, 'b', fs, 'Original data');
PlotECG(s_pca, 8, 'g', fs, 'PCA Results');
figure; semilogy(lambdas, 'bo'); grid; title('eigenvalue distribution plot');
PlotECG(s_jade, 8, 'r', fs, 'JADE Results');
PlotECG(s_sobi, 8, 'm', fs, 'SOBI Results');
PlotECG(s_fastica, 8, 'c', fs, 'FastICA Results');
PlotECG(s_pica_m, 8, 'm', fs, 'Maternal-based PiCA Results');
PlotECG(s_pica_f, 8, 'k', fs, 'Fetal-based PiCA Results');
PlotECG(s_pica_mf, 8, 'b', fs, 'Maternal-Fetal based PiCA Results');


