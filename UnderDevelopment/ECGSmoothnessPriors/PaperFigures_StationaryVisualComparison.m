clear
close all
clc

% load ECG1
% s = data;

% load ECG2
% s = data(14, 1:end);

% s = randn(1, 10000); fs = 1000;

% load ECG3
% s = data(2, 1:end);

addpath 'G:\Sameni\Documents\Papers\Journal 40 (Spline ECG Filtering)\CODES'
pth = 'J:\ECGData\Arrhythmia';
fs = 360;

wd = cd;
cd(pth);
fname = '233';
segstart = 11.2 + .6;
segstop = 17.2 + .375;
system(['rdsamp -r ' fname ' -f ' num2str(segstart) ' -t ' num2str(segstop) ' > sampledata.txt']);
data_all_channels = load([pth '\sampledata.txt'])';
s = data_all_channels(2,:);
s = s/norm(s);

% s = resample(s, 0.3*fs, fs);
% fs = 0.3*fs;
s = LPFilter(s - LPFilter(s, 2.0/fs), 100.0/fs);
s = s(:)';

figure
plot(s')
grid 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% warning('off', 'MATLAB:nearlySingularMatrix'); % THIS IS JUST FOR AVOIDING THE MATRIX SINGULARITY WARNINGS. A MORE NUMERICALLY STABLE VERSION OF THE CODE MIGHT BE NEEDED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
randn('seed', 0);
rand('seed', 0);

wlen = 100e-3; % window length (s)
f0 = 1.0; % Heart Beart (Hz)
DiffOrder = 2; % smoothness constraint order greater than 1
% lambda2 = 500*.5e-3*sqrt(fs^DiffOrder); % smoothness factor
% guardlen_l = 2*DiffOrder; % lower guard window lengh, >= DiffOrder
% guardlen_u = 2*DiffOrder; % upper guard window lengh, >= DiffOrder
snr = 3; % dB
% % % % adapt = 0; % adaptive(1), naive adaptation(2), or fixed(0) lambda
withwindow = 0; % with(1) or without(0) windowing
mode = 2; % 0..5

nvar = var(s)/10^(snr/10);
% noise =  NoiseGeneratorSimple('WHITE', nvar, fs, length(s));
% noise = noise.*cos(2*pi*(1:length(noise))*0.3/fs).*cos(2*pi*(1:length(noise))*0.21/fs);
% x = s + noise;
x = s + sqrt(nvar)*randn(size(s));


% lambda2 = 1e6*nvar;
% lambda2 = sqrt(nvar);
% [x_filtered1, x_filtered2] = ECGSmoothingPriors(x, DiffOrder, wlen*fs, lambda, guardlen_l, guardlen_u, adapt, withwindow);
% % % [x_filtered1, x_filtered2] = ECGSmoothnessPriorsDenoiser(x, DiffOrder, round(wlen*fs), 0, lambda2, adapt);
% [x_filtered1, x_filtered2] = ECGSmoothnessPriorsDenoiserBW(x, lambda2, mode, DiffOrder, round(wlen*fs));
[x_filtered1, x_filtered2, optim_gammas1, optim_gammas2, c1, e1, c2, e2, L_curveC1, L_curveE1, L_curveC2, L_curveE2] = ECGSmoothnessPriorsDenoiserBW(x, 1*(nvar), mode, DiffOrder, round(wlen*fs), 1e-8, 250, 1, 500);
[x_filtered11, x_filtered22, optim_gammas11, optim_gammas22, c11, e11, c22, e22, L_curveC11, L_curveE11, L_curveC22, L_curveE22] = ECGSmoothnessPriorsDenoiserBW(x, 2.0*(nvar), mode, DiffOrder, round(wlen*fs), 1e-8, 250, 1, 250);
% [x_KF,x_KS,Pbar,Phat,PSmoothed,Kgain,a, x_HI, PPhat, KKgain, margin] = ECGSmoothnessPriorsDenoiserKF(x, 2, 1*(nvar), 1);%, b, a, beta, var(x - BPFilter(x,4/fs,25/fs)), gamma, wlen, 'obsrv');
[x_KF,x_KS,Pbar,Phat,PSmoothed,Kgain,a] = ECGSmoothnessPriorsDenoiserKF(x, DiffOrder, 1*(nvar), .01*var(s), 1, 1, round(wlen*fs));%, b, a, beta, var(x - BPFilter(x,4/fs,25/fs)), gamma, wlen, 'obsrv');
x_LTI = ECGSmoothnessPriorsDenoiserLTI(x, .05*(nvar), 1, DiffOrder);
% x_LTI = ECGSmoothnessPriorsDenoiserLTI(x, max(optim_gammas2), 0, DiffOrder);

tt = (0:length(x)-1)/fs;

snr0 = 10*log10(sum(s.^2)/sum((s - x).^2));
snr1 = 10*log10(sum(s.^2)/sum((s - x_filtered1).^2));
snr2 = 10*log10(sum(s.^2)/sum((s - x_filtered2).^2));
snr3 = 10*log10(sum(s.^2)/sum((s - x_KF).^2));
snr4 = 10*log10(sum(s.^2)/sum((s - x_KS).^2));
snr5 = 10*log10(sum(s.^2)/sum((s - x_LTI).^2));

disp(['snr_input = ' num2str(snr0) '(dB)']);
disp(['snr_filtered1 = ' num2str(snr1) '(dB)']);
disp(['snr_filtered2 = ' num2str(snr2) '(dB)']);
disp(['snr_KF = ' num2str(snr3) '(dB)']);
disp(['snr_KS = ' num2str(snr4) '(dB)']);
disp(['snr_LTI = ' num2str(snr5) '(dB)']);

% peaks_0 = PeakDetection(s, f0/fs);
% I_0 = find(peaks_0);
% peaks_1 = PeakDetection(x_filtered1, f0/fs);
% I_1 = find(peaks_1);
% peaks_2 = PeakDetection(x_filtered2, f0/fs);
% I_2 = find(peaks_2);

% figure
% mesh(inv(eye(size(D,2)) + (lambda)^2*(D'*D)));
% mesh(inv(eye(size(Dd,2)) + (lambda)^2*(Dd'*Dd)));

% figure
% subplot(211);
% plot(tt, s, 'k');
% hold on
% plot(tt(I_1), s(I_1), 'bo');
% plot(tt(I_2), s(I_2), 'ro');
% grid
% subplot(212);
% HR_0 = 60*fs./diff(I_0);
% HR_1 = 60*fs./diff(I_1);
% HR_2 = 60*fs./diff(I_2);
% plot(tt(I_0(2:end)), HR_0, 'ko');
% hold on
% plot(tt(I_1(2:end)), HR_1, 'b');
% plot(tt(I_2(2:end)), HR_2, 'r');
% grid

% figure
% hold on
% plot(tt, x, 'c');
% plot(tt, s, 'k', 'linewidth', 2);
% plot(tt, x - x_filtered1, 'b');
% plot(tt, x - x_filtered2, 'r');
% grid

figure
hold on
plot(optim_gammas1, 'b');
plot(optim_gammas2, 'r');
grid

bias = .015;
% figure('position', [73 274 768 425]);  % create new figure with specified size  [500, 200, 500, 500] [95 283 561 425]
figure('Position', [304 259 617 425]);
hold on;
plot(tt, x, 'color', .5*ones(1,3), 'linewidth', 2);
plot(tt, s, 'k', 'linewidth', 2)
% plot(tt, x_filtered1 - bias, 'b', 'linewidth', 2);
plot(tt, x_filtered2 - 1*bias, 'b', 'linewidth', 2);
% plot(tt, x_filtered22 - 2*bias, 'r', 'linewidth', 2);
plot(tt, x_LTI - 2*bias, 'r', 'linewidth', 2);
% plot(tt, x_KF - 3*bias, 'm', 'linewidth', 2);
% plot(tt, x_KS - 4*bias, 'g', 'linewidth', 2);
% plot(tt, 350*optim_gammas2/max(optim_gammas2) - 3*bias, 'm', 'linewidth', 2);
% 
plot(tt, .003*optim_gammas2 - 2.8*bias, 'm', 'linewidth', 2);
% plot(tt, .00001*optim_gammas2 - 2.8*bias, 'm', 'linewidth', 2);

% bb = -bias*(0:5)'*ones(1, length(x));
% bb = -bias*(0:2)'*ones(1, length(x));
% plot(tt, bb, '--', 'color', .5*ones(1,3))
set(gca, 'YTickLabel', [])
set(gca, 'YTick', []);
set(gca, 'fontsize', 16);
xlabel('time(s)', 'fontsize', 16);
axis tight

figure;
subplot(121);
plot((c1), (e1), 'bo');
% hold on;
% plot((c2), (e2), 'ro');
grid
subplot(122);
plot((c11), (e11), 'bo');
% hold on;
% plot((c22), (e22), 'ro');
grid

   
figure
hold on
for i = 1 : size(L_curveC2, 2)
    plot(sqrt(squeeze(L_curveE2(1, i, :))), sqrt(squeeze(L_curveC2(1, i, :))), 'b', 'linewidth', 1);
end
plot(sqrt(e2), sqrt(c2), 'ro', 'linewidth', 1);
set(gca, 'fontsize', 16);
xlabel('$\|{D}_\theta {\theta}_k + {b}_k\|$', 'fontsize', 16, 'interpreter', 'latex');
ylabel('$\|{x}_k - {\theta}_k\|$', 'fontsize', 16, 'interpreter', 'latex');
% axis([0 round(4*median(sqrt(e2))/100)*100 0 round(3*median(sqrt(c2))/100)*100]);
% axis([0 4*median(sqrt(e2)) 0 3*median(sqrt(c2))]);
set(gca, 'box', 'on');
% set(gca, 'Position', [0.160714 0.142857 0.782143 0.782143]);
% set(gca, 'XScale', 'log')
% set(gca, 'YScale', 'log')


% set(gcf,'PaperPositionMode','auto');
% set(gcf, [426 201 550 581]);
% set(gca, 'Position', [0.13 0.011 0.775 0.815]);
% set(gca, 'OuterPosition', [0 1.38778e-017 1 1]);
% aa = axis;
% aa(1) = 0;
% aa(2) = tt(end);
% axis(aa);

% set(gca, 'YGrid', 'off')
% figure
% plot(a');
% grid
% 
% figure
% plot(Kgain');
% grid
% 
% figure
% hold on
% plot(squeeze(Pbar(1,1,:)), 'b');
% plot(squeeze(Phat(2,2,:)), 'r');
% plot(squeeze(PSmoothed(2,2,:)),'c');
% grid
% 
% figure
% hold on
% plot(squeeze(Pbar(1,1,:)), squeeze(Pbar(2,2,:)), 'b.');
% plot(squeeze(Phat(1,1,:)), squeeze(Phat(2,2,:)), 'r.');
% plot(squeeze(PSmoothed(1,1,:)), squeeze(PSmoothed(2,2,:)),'c.');
% grid
% 
% figure
% plot(margin');
% grid
