% *************************************************************************
% * The main code representing the "Robust Statistical Framework for      *
% * Instantaneous EEG Phase and Frequency Analysis" including phase       *
% * calculation, SNR verification and post processing (Kalman Smoothing)  *
% * >>> Refer to the user manual and reference [2] for more detaiils.     *
% *************************************************************************
% 
% Dependencies: -The Cerebral Signal Phase Analysis Toolbox of Open Source 
%                Electrophysiological Toolbox
%               -functions 'BPFilter5.m', 'KFNotch.m', 'LPFilter.m' and 
%                'KalmanSmoother.m' from "General filtering and
%                processing tools" of OSET
% 
% Please make sure to reference BOTH the original studies [1-2] and the 
% OSET [3] to help others find these items.
% 
%     [1] Esmaeil Seraj, Reza Sameni. ”Robust Electroencephalogram Phase 
%         Estimation with Applications in Brain-computer Interface Systems” 
%         Physiological Measurements (2017)
%     [2] Reza Sameni and Esmaeil Seraj, “A Robust Statistical Framework 
%         for Instantaneous Electroencephalogram Phase and Frequency 
%         Analysis” Physiological Measurements (2017)                     
%     [3] R. Sameni, The Open-Source Electrophysiological Toolbox (OSET), 
%         version 3.1 (2014). URL http://www.oset.ir
%         Released under the GNU General Public License
%         Copyright (C) 2012  Reza Sameni
%         Shiraz University, Shiraz, Iran
%         reza.sameni@gmail.com 
% 
% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details. You should have received a copy of the
% GNU General Public License along with this program; if not, write to the
% Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
% MA  02110-1301, USA.
% 

close all
clear
clc

%%--------------initialization and loading real EEG signal---------------%%
dat = load('EEG.mat');x = dat.record1; fs = 160;
f0 = 10;
f0_dev = 1e-6;
f0_up_down_dev = 1e-2;
bw_base = 2.0;
bw_base_dev = 0.1;
dither_std = 1.0e-4;
order = 6;
Itr = 10;

transient = round(3.5*fs);       % defining the transient range

s = x(3, round(10*fs):round(47*fs));
N = length(s);

%%--------synthetic signal (sinusoidal plus noise) with known SNR--------%% 
%%-----(uncomment the following section to use the synthetic signal)-----%%
% N = round(37*fs);
% snrin = 10;                     % Define the in-band SNR here
% s0 = sin(2*pi*f0/fs*(1:N)) + 0.1*sin(2*pi*(f0*0.99)/fs*(1:N) - pi/7) + 0.3*sin(2*pi*(f0*1.01)/fs*(1:N) - pi/3);
% nvr = var(s0)/10^(snrin/10);
% nstd = sqrt(nvr);
% n0_temp = BPFilter5(randn(1, N), f0/fs, bw_base/fs, order);
% s = s0 + nstd*n0_temp/std(n0_temp);

%%-------------Uncomment if Needed (Not required within code)------------%%
% % test noise variance:
% x = randn(1, N);
% xx = BPFilter5(s, f0/fs, bw_base/fs, order);
% xx = xx/std(xx);
% xa = hilbert(xx);
% std(xa)

% s = 10*sin(2*pi*(f0 + cumsum(1*ones(1,N)/N))/fs.*(1:N));

% s_BP1 = BPFilter5(s, f0/fs, bw_base/fs, order);
% s_BP2 = BPFilter5(s, f0/fs, (bw_base + 0.1)/fs, order);
% s_BP3 = BPFilter5(s, f0/fs, (bw_base - 0.1)/fs, order);

% s = BPFilter5(s, f0/fs, bw_base/fs, order); % make the signal narrowband

%*************************************************************************%
%*************************************************************************%
%*************************************************************************%

%%---------------------------Phase Estimation----------------------------%%

s_BP = BPFilter5(s, f0/fs, bw_base/fs, order);
s_BP_up = BPFilter5(s, (f0 + f0_up_down_dev)/fs, bw_base/fs, order);
s_BP_down = BPFilter5(s, (f0 - f0_up_down_dev)/fs, bw_base/fs, order);

s_a = hilbert(s_BP);
s_a_up = hilbert(s_BP_up);
s_a_down = hilbert(s_BP_down);

% xx = zeros(Itr, length(s), 4);
xx_dith = zeros(Itr, length(s));
xx_f = zeros(Itr, length(s));
xx_bw = zeros(Itr, length(s));
xx_all = zeros(Itr, length(s));
ph = zeros(Itr, length(s));
ph2 = zeros(Itr, length(s));
freq = zeros(Itr, length(s));
bw = zeros(Itr, 1);
f = zeros(Itr, 1);
effective_power_bw = zeros(Itr, 1);
effective_power_all = zeros(Itr, 1);
n0 = randn(1, N);
n0_band_limited = hilbert(BPFilter5(n0, f0/fs, bw_base/fs, order));
n0_var = var(n0_band_limited);
ph_linear = 2*pi*f0/fs*(0:N-1);
for k = 1 : Itr,
    dither_narrow_band = BPFilter5(randn(1, N), f0/fs, bw_base/fs, order);
    xx_dith(k, :) = hilbert( BPFilter5(s, f0/fs, bw_base/fs, order) + dither_std*dither_narrow_band/std(dither_narrow_band)); % dithered signal
    
%     bw(k) = bw_base + bw_base_dev*(2*rand - 1);                            % deviations in both sides
    bw(k) = bw_base + bw_base_dev*rand;                                      % only wider deviations
    f(k) = f0 + f0_dev*(2*rand - 1);
    
    xx_f(k, :) = hilbert( BPFilter5(s , f(k)/fs, bw_base/fs, order) );       % f0 randomization
    
    noise_narrow_band_bw = hilbert(BPFilter5(n0, f0/fs, bw(k)/fs, order)) - n0_band_limited;
    noise_narrow_band_all = hilbert(BPFilter5(n0, f(k)/fs, bw(k)/fs, order)) - n0_band_limited;
    
    xx_bw(k, :) = hilbert( BPFilter5(s , f0/fs, bw(k)/fs, order));           % BW randomization % 
    xx_all(k, :) = hilbert( BPFilter5(s , f(k)/fs, bw(k)/fs, order) + dither_std*dither_narrow_band/std(dither_narrow_band)); % All randomization methods
    
    effective_power_bw(k) = var(noise_narrow_band_bw)/n0_var;
    effective_power_all(k) = var(noise_narrow_band_all)/n0_var;

    %     xx(k, :, 4) = hilbert( BPFilter5(s + dither_std*randn(1, N), f0/fs, bw_base/fs, order) ); % dithered signal
        ph(k, :) = unwrap(atan2(imag(xx_all(k, :)), real(xx_all(k, :))));
        freq(k, :) = fs*diff([ph(k, 1) ph(k, :)])/(2*pi);
        ph2(k, :) = ph(k, :) - ph_linear;
end

%%---------------------removing the transient effect---------------------%%
s = s(transient + 1 : end-transient);
s_BP = s_BP(transient + 1 : end-transient);
s_BP_up = s_BP_up(transient + 1 : end-transient);
s_BP_down = s_BP_down(transient + 1 : end-transient);
s_a = s_a(transient + 1 : end-transient);
s_a_up = s_a_up(transient + 1 : end-transient);
s_a_down = s_a_down(transient + 1 : end-transient);

ph = ph(:, transient + 1 : end-transient);
ph2 = ph2(:, transient + 1 : end-transient);
ph2 = ph2 - ph2(1);
freq = freq(:, transient + 1 : end-transient);

%*************************************************************************%
%*************************************************************************%
%*************************************************************************%

%%--------------SNR Estimation for All Randomization Methods-------------%%

ph_mean = mean(ph, 1);
ph_std = std(ph, [], 1);

ph2_mean = mean(ph2, 1);
ph2_std = std(ph2, [], 1);

freq_mean = mean(freq, 1);
freq_std = std(freq, [], 1);

xx_dith = xx_dith(:, transient + 1 : end-transient);
xx_dith_mean = mean(xx_dith, 1);
noise_dith = xx_dith - ones(Itr, 1)*xx_dith_mean;
noise_dith_normalized = noise_dith;
noise_dith_normalized_var = var(noise_dith_normalized, [], 1);
snr_dith = 10*log10((abs(xx_dith_mean).^2/2)./(noise_dith_normalized_var - dither_std^2/2) );

xx_f = xx_f(:, transient + 1 : end-transient);
xx_f_mean = mean(xx_f, 1);
noise_f = xx_f - ones(Itr, 1)*xx_f_mean;
noise_f_normalized = noise_f/n0_var;
noise_f_normalized_var = nanvar(noise_f_normalized, [], 1);
snr_f = 10*log10((abs(xx_f_mean).^2/2)./noise_f_normalized_var);

xx_bw = xx_bw(:, transient + 1 : end-transient);
xx_bw_mean = mean(xx_bw, 1);
noise_bw = xx_bw - ones(Itr, 1)*s_a;
noise_bw_normalized = noise_bw./sqrt(effective_power_bw(:, ones(1, length(noise_bw))));
noise_bw_normalized_var = nanvar(noise_bw_normalized, [], 1);

snr_bw = 10*log10((abs(s_a).^2/2)./nanmean(abs(noise_bw_normalized).^2));

xx_all = xx_all(:, transient + 1 : end-transient);
xx_all_mean = mean(xx_all, 1);
noise_all = xx_all - ones(Itr, 1)*s_a;
noise_all_normalized = noise_all./sqrt(effective_power_all(:, ones(1, length(noise_all))));
noise_all_normalized_var = nanvar(noise_all_normalized, [], 1);
snr_all = 10*log10((abs(xx_all_mean).^2/2)./(nanmean(abs(noise_all_normalized).^2) - dither_std^2/2));% - 10*log10(bw_base/bw_base_dev); % - 10*log10(fs/bw_base_dev)

%*************************************************************************%
%*************************************************************************%
%*************************************************************************%

%%-----------------Post Processing and Kalman Smoothing------------------%%

wlen = 50;
R = filtfilt(ones(1, wlen), wlen, freq_std).^2;
[yf,ys,Pbar,Phat,PSmoothed,Kgain, innovations] = KalmanSmoother(freq_mean, 1, [1 -1], 1e-5, R);
[yf2,ys2,Pbar2,Phat2,PSmoothed2,Kgain2, innovations2] = KalmanSmoother(freq_mean, 1, [1 -1], 1e-3, R);
% [yf3,ys3,Pbar3,Phat3,PSmoothed3,Kgain3, innovations3] = KalmanSmoother(freq_mean, 1, [1 -1], 1e-1, R);

[y1, y2,Pbar,Phat,PSmoothed,Kgain] = KFNotch(s_BP, f0, fs, 1, 1000*var(s_BP),1);

%*************************************************************************%
%*************************************************************************%
%*************************************************************************%

%%----------------------visualizing the results--------------------------%%

% xx_real_mean = mean(real(xx), 1);
% xx_imag_mean = mean(imag(xx), 1);
% % % xx_real_std = std(real(xx)./sqrt(bw(:, ones(1, N))/bw_base), [], 1);
% % % xx_imag_std = std(imag(xx)./sqrt(bw(:, ones(1, N))/bw_base), [], 1);
% xx_real_std = std(real(xx)./sqrt(effective_power_bw(:, ones(1, N))/n0_var), [], 1);
% xx_imag_std = std(imag(xx)./sqrt(effective_power_bw(:, ones(1, N))/n0_var), [], 1);

% % xx_mean = xx_real_mean + 1j*xx_imag_mean;
% xx_std = (xx_real_std + 1j*xx_imag_std)/sqrt(2);

% ph_mean = mean(ph, 1);
% ph_std = std(ph, [], 1);

% t = (0 : length(xx_f)-1)/fs;
t = (0 : length(xx_bw)-1)/fs;

figure
hold on
plot(s_BP);
plot(s_BP - y1, 'r');
plot(s_BP - y2, 'k');
grid
title('Smoothed s_{BP}');

dfft = 500;
figure
hold on
plot(10*log10(abs(fft(s_BP, dfft))));
plot(10*log10(abs(fft(s_BP - y1, dfft))), 'r');
plot(10*log10(abs(fft(s_BP - y2, dfft))), 'k');
grid
title('Smoothed s_{BP}');

%--Full view (remove the 'position' from figures if figures appear out/partially-out of your monitor screen)--%
position = [201 433 863 328];
inds = 1 : length(t);

% arround 7.5s
% position = [201 433 341 328];
% inds = round(7.1*fs) : round(7.8*fs);

% arround 11.5s
% position = [201 433 341 328];
% inds = round(11.4*fs) : round(12.0*fs);

% % arround 22.5s
% position = [201 433 341 328];
% inds = round(22.5*fs) : round(23.2*fs);

figure
plot(t(inds), 0.01*s(inds), 'b');
grid
xlabel('time(s)', 'fontsize', 16);
ylabel('Amplitude(mV)', 'fontsize', 16);
set(gca, 'fontsize', 16);
set(gcf, 'Position', position);
axis tight

figure
hold on;
plot(t(inds), s_BP_up(inds), 'g', 'linewidth', 2);
plot(t(inds), s_BP_down(inds), 'r', 'linewidth', 2);
plot(t(inds), s_BP(inds), 'b');
plot(t(inds), abs(s_a_up(inds)), 'g--', 'linewidth', 2);
plot(t(inds), abs(s_a_down(inds)), 'r--', 'linewidth', 2);
plot(t(inds), abs(s_a(inds)), 'k', 'linewidth', 2);
% % % plot(t, abs(s_a_up) - abs(s_a_down), 'm', 'linewidth', 3);
% % % plot(t, 10*(abs(s_a) - (abs(s_a_up) + abs(s_a_down))/2), 'c', 'linewidth', 3);
% % % plot(t(2:end), 100*diff(abs(s_a_up) - abs(s_a_down)), 'c', 'linewidth', 2);
% % % plot(t, (s_a - s_a_up)./s_a, 'm', 'linewidth', 2);
% % % plot(t, real(s_a_up - s_a_down), 'm', 'linewidth', 2);
% % % plot(t, imag(s_a_up - s_a_down), 'c', 'linewidth', 2);
% % % plot(t, angle(exp(1j*(s_a - s_a_up))), 'm', 'linewidth', 2);
% % % plot(t, angle(exp(1j*(s_a - s_a_down))), 'c', 'linewidth', 2);
grid;
xlabel('time(s)', 'fontsize', 16);
ylabel('Amplitude(mV)', 'fontsize', 16);
set(gca, 'fontsize', 16);
set(gcf, 'Position', position);
axis tight

% figure
% subplot(211)
% plot(t, real(noise_bw'))
% grid
% subplot(212)
% plot(t, imag(noise_bw'))
% grid

% % % figure
% % % subplot(211)
% % % plot(t, real(s_a))
% % % hold on
% % % plot(t, real(s(transient + 1 : end-transient)), 'r')
% % % grid
% % % subplot(212)
% % % plot(t, imag(s_a))
% % % hold on
% % % plot(t, imag(s(transient + 1 : end-transient)), 'r')
% % % grid

% figure
% hold on
% % plot(t, abs(noise_dith), 'color', 0.7*ones(1, 3))
% 
% % plot(t, abs(noise_f), 'color', 0.7*ones(1, 3));
% % plot(t, sqrt(noise_f_normalized_var), 'color', 0.2*ones(1, 3), 'linewidth', 2);
% 
% plot(t, abs(noise_bw_normalized), 'color', 0.7*ones(1,3));
% % plot(t, real(noise_bw_normalized), 'color', 0.7*[0 1 0]);
% % plot(t, imag(noise_bw_normalized), 'color', 0.7*[1 0 0]);
% plot(t, sqrt(noise_bw_normalized_var), 'color', 0.2*ones(1, 3), 'linewidth', 2);
% % plot(t, nstd, 'r', 'linewidth', 2);
% grid
% 
% figure
% hold on
% % plot(t, abs(xx_dith), 'color', 0.7*ones(1, 3))
% % plot(t, abs(xx_dith_mean), 'color', 0.3*ones(1, 3), 'linewidth', 2)
% % plot(t, abs(xx_f), 'color', 0.7*ones(1, 3))
% % plot(t, abs(xx_f_mean), 'color', 0.3*ones(1, 3), 'linewidth', 2)
% plot(t, abs(xx_bw), 'color', 0.7*ones(1, 3))
% plot(t, abs(s_a), 'color', 0.3*ones(1, 3), 'linewidth', 2)
% grid

% figure
% hold on
% plot(t, snr_dith - 0*max(snr_dith), 'b', 'linewidth', 2);
% plot(t, snr_f - 0*max(snr_f), 'r', 'linewidth', 2);
% plot(t, snr_bw - 0*max(snr_bw), 'k', 'linewidth', 2);
% plot(t, snr_all - 0*max(snr_all), 'y', 'linewidth', 2);
% grid;
% xlabel('time(s)', 'fontsize', 16);
% ylabel('SNR(dB)', 'fontsize', 16);
% set(gca, 'fontsize', 16);
% legend('dither', 'f0', 'bw', 'All methods');
% set(gcf, 'Position', position);
% axis tight

% a = axis;
% a(3) = 0;
% a(4) = 35;
% axis(a);

% smp1 = 5000;
% smp2 = 6000;
% smp3 = 4000;
% figure
% polar(ph(:, smp1), abs(xx_f(:, smp1)), 'bo');
% hold on
% polar(ph(:, smp2), abs(xx_f(:, smp2)), 'ro');
% polar(ph(:, smp3), abs(xx_f(:, smp3)), 'go');

figure
hold on
plot(t(inds), ph2(:,inds), 'color', 0.7*ones(1, 3), 'linewidth', 2);
plot(t(inds), ph2_mean(inds), 'color', 'k', 'linewidth', 2);
% plot(t, ph2_mean + ph2_std, 'b', 'linewidth', 1);
% plot(t, ph2_mean - ph2_std, 'r', 'linewidth', 1);
% plot(t, ph2_std, 'b', 'linewidth', 1);
grid;
xlabel('time(s)', 'fontsize', 16);
ylabel('Unwrapped Phase (rad)', 'fontsize', 16);
set(gca, 'fontsize', 16);
% a = axis;
% a(3) = 9;
% a(4) = 11;
% axis(a);
set(gcf, 'Position', position);
axis tight

figure
hold on
plot(t(inds), freq(:,inds), 'color', 0.7*ones(1, 3), 'linewidth', 2);
plot(t(inds), freq_mean(inds), 'color', 'k', 'linewidth', 2);
plot(t(inds), freq_mean(inds) + freq_std(inds), 'b--', 'linewidth', 1);
plot(t(inds), freq_mean(inds) - freq_std(inds), 'r--', 'linewidth', 1);
% plot(t, yf, 'color', 'm', 'linewidth', 2);
% plot(t(inds), ys(inds), 'color', 'c', 'linewidth', 2);
% plot(t, f0 + freq_std, 'b', 'linewidth', 1);
grid;
xlabel('time(s)', 'fontsize', 16);
ylabel('IF (Hz)', 'fontsize', 16);
set(gca, 'fontsize', 16);
axis tight
% a = axis;
% a(3) = 8;
% a(4) = 12;
% axis(a);
set(gcf, 'Position', position);

figure
hold on
plot(t(inds), freq(:,inds), 'color', 0.8*ones(1, 3), 'linewidth', 2);
plot(t(inds), freq_mean(inds), 'color', 'k', 'linewidth', 2);
% plot(t(inds), freq_mean(inds) + freq_std(inds), 'b--', 'linewidth', 1);
% plot(t(inds), freq_mean(inds) - freq_std(inds), 'r--', 'linewidth', 1);
% plot(t(inds), yf(inds), 'color', 'm', 'linewidth', 2);
% plot(t(inds), ys3(inds), 'color', 'c', 'linewidth', 2);
plot(t(inds), ys2(inds), 'color', 'b', 'linewidth', 2);
plot(t(inds), ys(inds), 'color', 'r', 'linewidth', 2);
% plot(t, f0 + freq_std, 'b', 'linewidth', 1);
grid;
xlabel('time(s)', 'fontsize', 16);
ylabel('IF (Hz)', 'fontsize', 16);
set(gca, 'fontsize', 16);
axis tight
a = axis;
% a(3) = 9.4;
% a(4) = 10.7;
% axis(a);
set(gcf, 'Position', [201 433 863 328]);
set(gca, 'Position', [0.0938586 0.164634 0.885284 0.760366]);

bias = 15;%15;
ndft = 20000;%length(s_BP);
ff = fs*(0:ndft-1)/ndft;
f1 = f0 - bw_base;%9.0;%7.2;%9.0
f2 = f0 + bw_base;%11.0;%8.3;%11.0
% plot_range = round(f1/fs*ndft):round(f2/fs*ndft);
% plot_range = 1:ndft;
plot_range = round(f1/fs*ndft):round(f2/fs*ndft);
window = hamming(length(s))';
S = fft(s.*window, ndft);
S_BP = fft(s_BP.*window, ndft);
S_STOP = fft((s - s_BP).*window, ndft);
% S_STOP = fft((s - 0.55*s_BP).*window, ndft);
figure
hold on
line([9.65 9.65], [0 65], 'linewidth', 3, 'linestyle', '--', 'color', 'k');
line([10.34 10.34], [0 65], 'linewidth', 3, 'linestyle', '--', 'color', 'k');
% line([7.53 7.53], [0 55], 'linewidth', 2, 'linestyle', '--', 'color', 'c');
% line([7.95 7.95], [0 55], 'linewidth', 2, 'linestyle', '--', 'color', 'c');
plot(ff(plot_range), 20*log10(abs(S(plot_range))) - bias, 'r', 'linewidth', 4);
plot(ff(plot_range), 20*log10(abs(S_BP(plot_range))) - bias, 'b', 'linewidth', 3);
plot(ff(plot_range), 20*log10(abs(S_STOP(plot_range))) - bias, 'g--', 'linewidth', 3);
% legend('S_1(f)', 'S_2(f)', 'S_3(f)');
% grid
xlabel('Frequency (Hz)', 'fontsize', 16);
ylabel('Normalized PSD (dB)', 'fontsize', 16);
set(gca, 'fontsize', 16);
a = axis;
a(1) = f1;
a(2) = f2;
a(3) = 10;
a(4) = 65;%55;
axis(a);
text(9.1, 62.5, 'Stop-band', 'fontsize', 16);
text(9.8, 62.5, 'Pass-band', 'fontsize', 16);
text(10.4, 62.5, 'Stop-band', 'fontsize', 16);
% text(7.21, 52.5, 'Stop-band', 'fontsize', 14, 'color', 'm');
% text(7.6, 52.5, 'Pass-band', 'fontsize', 14, 'color', 'm');
% text(8.0, 52.5, 'Stop-band', 'fontsize', 14, 'color', 'm');

% text(7.2, 20, 'Foreground+Background EEG', 'fontsize', 14);
% text(7.5, 20, 'Background EEG', 'fontsize', 14);
% text(8.0, 20, 'Foreground EEG', 'fontsize', 14);
set(gca, 'box', 'on');

% % % figure
% % % subplot(311);
% % % plot(t, Kgain)
% % % grid
% % % subplot(312);
% % % plot(t, squeeze(Pbar));
% % % hold on
% % % plot(t, squeeze(Phat), 'r');
% % % plot(t, squeeze(PSmoothed), 'k');
% % % grid
% % % subplot(313);
% % % plot(t, innovations);
% % % grid

% % % figure
% % % hold on
% % % plot(t, freq_mean, 'color', 'b', 'linewidth', 2);

% figure
% hold on
% plot(t, snr1 - 0*max(snr1), 'k', 'linewidth', 2);
% plot(t, snr2 - 0*max(snr2), 'b', 'linewidth', 2);
% plot(t, snr3 - 0*max(snr3), 'r', 'linewidth', 2);
% plot(t, snr4 - 0*max(snr4), 'g', 'linewidth', 2);
% grid;
% xlabel('time(s)', 'fontsize', 16);
% ylabel('SNR(dB)', 'fontsize', 16);
% set(gca, 'fontsize', 16);
% legend('All', 'f_0', 'bw', 'dither');
%
% figure
% plot(t, squeeze(abs(xx(:, :, 1))), 'color', 0.65*ones(1,3));
% hold on
% plot(t, abs(xx_mean_1), 'k', 'linewidth', 2);
% % plot(t, noise_normalized_std_1, 'r');
% plot(t, noise_normalized_std_2, 'b');
% plot(t, noise_normalized_std_3, 'r');
% plot(t, noise_normalized_std_4, 'g');
% % plot(t, abs(xx_mean) + noise_normalized_std, 'r');
% % plot(t, abs(xx_mean) - noise_normalized_std, 'g');
% grid
%
%
% 
% % figure
% % subplot(411);
% % subplot(412);
% % plot(t, imag(xx_mean + 2*xx_std), 'r');
% % hold on
% % plot(t, imag(xx_mean - 2*xx_std), 'g');
% % plot(t, imag(xx_mean));
% % grid
% %
% % subplot(413);
% % plot(t, snr);
% % grid
% %
% % subplot(414);
% % plot(t, ph_mean + 2*ph_std, 'r');
% % hold on
% % plot(t, ph_mean - 2*ph_std, 'g');
% % plot(t, ph_mean);
% % grid
% %
% % figure
% % subplot(311);
% % plot(t, abs(xx_mean));
% % hold on
% % plot(t, abs(xx_std), 'r');
% % grid
% %
% % subplot(312);
% % plot(t, snr);
% % grid
% %
% % subplot(313);
% % plot(t, ph_mean + 2*ph_std, 'r');
% % hold on
% % plot(t, ph_mean - 2*ph_std, 'g');
% % plot(t, ph_mean);
% % grid
%
% 
figure
spectrogram(s - LPFilter(s, 0.5/fs), hamming(256), 250, 512, fs,'yaxis');
%
% ndft = N;
% ff = fs*(0:ndft-1)/ndft;
% figure
% hold on
% for k = 1 : Itr,
%     %     psd(xx(k, :, 1), length(xx), fs);
%     plot(ff, 20*log10(abs(fft(xx(k, :, 1), ndft))));
% end
% grid
% plot(ff, 20*log10(abs(fft(s_a, ndft))), 'r');
% % psd(s_BP, 1024, fs);
%
% snr0 = 10*log10(mean(abs(xx_mean_1).^2/2)/mean(noise_normalized_var_1))
% 
