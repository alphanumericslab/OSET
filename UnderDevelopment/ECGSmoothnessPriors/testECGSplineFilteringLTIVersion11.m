clear
close all
clc

% load ECG1
% s = data;

load ECG2
s = data(5, 1:end);

% s = randn(1, 10000); fs = 1000;

% load ECG3
% s = data(2, 1:end);

% s = resample(s, 0.3*fs, fs);
% fs = 0.3*fs;
s = LPFilter(s - LPFilter(s, 1.0/fs), 100.0/fs);
s = s(:)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% warning('off', 'MATLAB:nearlySingularMatrix'); % THIS IS JUST FOR AVOIDING THE MATRIX SINGULARITY WARNINGS. A MORE NUMERICALLY STABLE VERSION OF THE CODE MIGHT BE NEEDED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
randn('seed', 0);
rand('seed', 0);

wlen = 150e-3; % window length (s)
f0 = 1.0; % Heart Beart (Hz)
DiffOrder = 2; % smoothness constraint order greater than 1
lambda = 5;%*.5e-3*sqrt(fs^DiffOrder); % smoothness factor
guardlen_l = 2*DiffOrder; % lower guard window lengh, >= DiffOrder
guardlen_u = 2*DiffOrder; % upper guard window lengh, >= DiffOrder
snr = 15; % dB
adapt = 1; % adapt(1) or not(0) lambda
withwindow = 0; % with(1) or without(0) windowing

nvar = var(s)/10^(snr/10);
x = s + sqrt(nvar)*randn(size(s));

[x_filtered1, x_filtered2] = ECGSmoothingPriors(x, DiffOrder, wlen*fs, lambda, guardlen_l, guardlen_u, adapt, withwindow);
% [x_KF,x_KS,Pbar,Phat,PSmoothed,Kgain,a, x_HI, PPhat, KKgain, margin] = KalmanSmoothingPriorsFilter(x, .1*nvar);%, b, a, beta, var(x - BPFilter(x,4/fs,25/fs)), gamma, wlen, 'obsrv');

%% Here's the LTI version
d = diff([zeros(1, DiffOrder) 1 zeros(1, DiffOrder)], DiffOrder)
% load DiffCoefs h
% d = conv(h, [1 -2 1]);

% phi = conv(d,d)
phi = conv(d,d(end:-1:1))
gg = phi*lambda^2;
midindex = (length(phi)+1)/2;
gg(midindex) = gg(midindex) + 1;
r = roots(gg)
r_abs = abs(r)
I_causal = find(r_abs < 0.9999);
r_causal = r(I_causal)
p_causal = poly(r_causal);
% gain = mean(poly(roots(g))./g);
% I_anticausal = find(r_abs > 1.001);
% r_anticausal = 1./r(I_anticausal)
% p_anticausal = poly(r_anticausal);

%% Here's the multiple difference LTI version
% d1 = diff([zeros(1, DiffOrder) 1 zeros(1, DiffOrder)], DiffOrder);
% d2 = diff([zeros(1, DiffOrder+1) 1 zeros(1, DiffOrder+1)], DiffOrder+1);
% d3 = diff([zeros(1, DiffOrder+2) 1 zeros(1, DiffOrder+2)], DiffOrder+2);
% lambda1 = lambda;%*sqrt(fs^DiffOrder);
% lambda2 = lambda*2;%^*sqrt(fs^(DiffOrder-1));
% lambda3 = lambda*4;%^*sqrt(fs^(DiffOrder-1));
% phi1 = conv(d1,d1);
% g1 = phi1*lambda1^2;
% phi2 = conv(d2,d2);
% g2 = phi2*lambda2^2;
% phi3 = conv(d3,d3);
% g3 = phi3*lambda3^2;
% gg = g3 + [0 g2 0] + [0 0 g1 0 0];
% midindexx = (length(gg)+1)/2;
% gg(midindexx) = gg(midindexx) + 1;
% rr = roots(gg)
% rr_abs = abs(rr)
% I_causal = find(rr_abs < 0.9999);
% r_causal = rr(I_causal)
% p_causal = poly(r_causal);

% figure
% stem(g);
% grid

% [Q,R] = deconv(g, 1)
% y = filter(1, g, x);

% z = filter(1, p_causal, x);
% y = filter(1, p_anticausal, z(end:-1:1));y = y(end:-1:1);

y = filtfilt(sum(p_causal), p_causal, x);
% y = filter([zeros(1, midindexx) 1], gg, x);

% y = y/lambda^2;
% y = y/std(y)*std(s);

% y = y*(1 + lambda^2*sum(phi)^2);

% y = .85*y;
% y = y/sqrt(mean(p_causal.^2));

tt = (0:length(x)-1)/fs;

snr0 = 10*log10(sum(s.^2)/sum((s - x).^2));
snr1 = 10*log10(sum(s.^2)/sum((s - x_filtered1).^2));
snr2 = 10*log10(sum(s.^2)/sum((s - x_filtered2).^2));
% snr3 = 10*log10(sum(s.^2)/sum((s - x_KF).^2));
% snr4 = 10*log10(sum(s.^2)/sum((s - x_KS).^2));
% snr5 = 10*log10(sum(s.^2)/sum((s - x_HI).^2));
snr6 = 10*log10(sum(s.^2)/sum((s - y).^2));

disp(['snr_input = ' num2str(snr0) '(dB)']);
disp(['snr_filtered1 = ' num2str(snr1) '(dB)']);
disp(['snr_filtered2 = ' num2str(snr2) '(dB)']);
% disp(['snr_KF = ' num2str(snr3) '(dB)']);
% disp(['snr_KS = ' num2str(snr4) '(dB)']);
% disp(['snr_HI = ' num2str(snr5) '(dB)']);
disp(['snr_LTI = ' num2str(snr6) '(dB)']);

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
% plot(tt, x - y, 'y');
% grid

sys = tf(gg, 1, 1/fs);
figure
pzplot(sys);

figure
freqz(gg, 1, [], fs);

figure
freqz(1, gg, [], fs);

figure
freqz(sum(p_causal), p_causal, [], fs);

figure
hold on;
plot(tt, x, 'c');
plot(tt, s, 'k', 'linewidth', 2)
plot(tt, x_filtered1, 'b', 'linewidth', 2);
plot(tt, x_filtered2, 'r', 'linewidth', 2);
plot(tt, y, 'y', 'linewidth', 2);
grid

% figure
% hold on
% plot(tt, x, 'c');
% plot(tt, s, 'k', 'linewidth', 2)
% % plot(tt, x_KF, 'b', 'linewidth', 2);
% % plot(tt, x_KS, 'r', 'linewidth', 2);
% % plot(tt, x_HI, 'm', 'linewidth', 2);
% plot(tt, y, 'y', 'linewidth', 2);
% grid

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