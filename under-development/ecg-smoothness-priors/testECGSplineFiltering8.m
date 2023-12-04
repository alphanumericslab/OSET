clear
close all
clc

load ECG2
s = data(14, 1:end);
% s = resample(s, 0.3*fs, fs);
% fs = 0.3*fs;
s = LPFilter(s - LPFilter(s, 0.3/fs), 40.0/fs);
s = s(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% warning('off', 'MATLAB:nearlySingularMatrix'); % THIS IS JUST FOR AVOIDING THE MATRIX SINGULARITY WARNINGS. A MORE NUMERICALLY STABLE VERSION OF THE CODE MIGHT BE NEEDED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
randn('seed', 0);
rand('seed', 0);

wlen = 150e-3; % window length (s)
f0 = 1.0; % Heart Beart (Hz)
DiffOrder = 3; % smoothness constraint order greater than 1
lambda = 10e-3*sqrt(fs^DiffOrder); % smoothness factor
guardlen_l = DiffOrder; % lower guard window lengh, >= DiffOrder
guardlen_u = DiffOrder; % upper guard window lengh, >= DiffOrder
snr = -10; % dB
adapt = 1; % adapt(1) or not(0) lambda

nvar = var(s)/10^(snr/10);
x = s + sqrt(nvar)*randn(size(s));

x = x(:);
knots1 = 1: round(wlen*fs) : length(x);
% add the first and last samples as knots to avoid end-point discontinuities
if(knots1(1) > 1)
    knots1 = [1 knots1];
end
if(knots1(end) < length(x))
    knots1 = [knots1 length(x)];
end

% first round
x_filtered1 = zeros(size(x));
for n = 1 : length(knots1) - 1,
    i = knots1(n) : knots1(n+1);
    ll = length(i);
    I = speye(ll);
    Dd = spdiags(ones(ll-2, 1)*diff([zeros(1, DiffOrder) 1 zeros(1, DiffOrder)], DiffOrder), 0:DiffOrder, ll-DiffOrder, ll);
    x_filtered1(i) = (I + lambda^2*(Dd'*Dd))\x(i);
end

% second round
knots2 = round(wlen*fs/2): round(wlen*fs) : length(x);
% II = II - round(wlen*fs/2);II(II < 1) = [];
if(knots2(1) > 1)
    knots2 = [1 knots2];
end
if(knots2(end) < length(x))
    knots2 = [knots2 length(x)];
end

x_filtered2 = zeros(size(x));
for n = 1 : length(knots1) - 1,
    i = knots2(n) : knots2(n+1); % I used this line to preserve the continuity of the curves
    
    ll = length(i);
    I = speye(ll);
    Dd = spdiags(ones(ll -2 + guardlen_l + guardlen_u, 1)*diff([zeros(1, DiffOrder) 1 zeros(1, DiffOrder)], DiffOrder), 0:DiffOrder, ll + guardlen_l + guardlen_u - DiffOrder, ll + guardlen_l + guardlen_u);
    D_l = Dd(:, 1:guardlen_l);
    D = Dd(:, guardlen_l+1:end-guardlen_u);
    D_u = Dd(:, end-guardlen_u+1:end);
    indxl = i(1)-guardlen_l : i(1)-1;
    indxl(indxl < 1) = 1;
    
    indxu = i(end) + 1 : i(end) + guardlen_u;
    indxu(indxu > length(x_filtered1)) = length(x_filtered1);
    
    l = x_filtered1(indxl);
    u = x_filtered1(indxu);
    if(adapt)
        lambda_adaptive = (lambda/(sqrt(length(x)/length(x(i)))*norm(x(i))/norm(x))); % adaptive version
    else
        lambda_adaptive = lambda; % fixed version
    end
    x_filtered2(i) = (I + lambda_adaptive^2*(D'*D))\(x(i) - lambda_adaptive^2*(D'*D_l(:, 1:length(l)))*l - lambda_adaptive^2*(D'*D_u(:, 1:length(u)))*u);
end

tt = (0:length(x)-1)/fs;

snr0 = 10*log10(sum(s.^2)/sum((s - x).^2));
snr1 = 10*log10(sum(s.^2)/sum((s - x_filtered1).^2));
snr2 = 10*log10(sum(s.^2)/sum((s - x_filtered2).^2));

disp(['snr0 = ' num2str(snr0) '(dB)']);
disp(['snr1 = ' num2str(snr1) '(dB)']);
disp(['snr2 = ' num2str(snr2) '(dB)']);

% peaks_0 = PeakDetection(s, f0/fs);
% I_0 = find(peaks_0);
% peaks_1 = PeakDetection(x_filtered1, f0/fs);
% I_1 = find(peaks_1);
% peaks_2 = PeakDetection(x_filtered2, f0/fs);
% I_2 = find(peaks_2);

figure
mesh(inv(eye(size(D,2)) + (lambda)^2*(D'*D)));
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

figure
hold on
plot(tt, x, 'c');
plot(tt, s, 'k', 'linewidth', 2);
plot(tt, x - x_filtered1, 'b');
plot(tt, x - x_filtered2, 'r');
grid

figure
hold on;
plot(tt, x, 'c');
plot(tt, s, 'k', 'linewidth', 2)
plot(tt, x_filtered1, 'b', 'linewidth', 2);
plot(tt, x_filtered2, 'r', 'linewidth', 2);
grid