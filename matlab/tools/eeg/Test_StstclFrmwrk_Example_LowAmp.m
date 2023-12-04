% *************************************************************************
% * This code is the implementation of an example showing the effects     *
% * of calculating phase and frequency in low-amplitude analytic signal,  *
% * as represented in [2].                                                *
% * >>> Refer to the user manual and reference [2] for more detaiils.     *
% *************************************************************************
% 
% Dependencies: -The Cerebral Signal Phase Analysis Toolbox of Open Source 
%                Electrophysiological Toolbox
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

%%-initialization and synthetic signal generation
randn('seed', 0);

%%-synthetic signal (1)
% fs=800;
% snr = 25;
% N = round(3*fs);
% f0 = 30;
% t = (0:N-1)/fs;
% f = 15*sin(2*pi*.8*t) + f0;
% A = mod(t, 0.7);
% s = A.*sin(2*pi*cumsum(f)/fs);
% n_var = var(s)/10^(snr/10);
% x = s + sqrt(n_var)*randn(1, N);
% I1 = round(1.41*fs):round(1.50*fs);
% I2 = round(1.6*fs):round(1.8*fs);
% ind1 = 28;
% ind2 = 130;

%%-synthetic signal (2)
fs = 400;
snr = 25;
N = round(2.0*fs);
f0 = 20;
t = (0:N-1)/fs;
f = 1*sin(2*pi*.8*t) + f0;
A = sin(2*pi*1.1*t) + 0.7*sin(2*pi*3.3*t + pi/3);
s = A.*sin(2*pi*cumsum(f)/fs);
n_var = var(s)/10^(snr/10);
x = s + sqrt(n_var)*randn(1, N);
I1 = round(0.41*fs):round(0.43*fs);
I2 = round(0.25*fs):round(0.52*fs);
ind1 = 5;
ind2 = 19;

%%-phase and frequency estimation
xa = hilbert(x);
ph = atan2(imag(xa), real(xa));
ph = unwrap(ph);
theta = ph - 2*pi*f0*t;
fhat = fs*diff(ph([1 1:N]))/2/pi;

%%-visualizing results
figure;
subplot(311);
plot(t, x, 'linewidth', 1);
hold on;
plot(t, abs(xa), 'r--', 'linewidth', 2);
grid
ylabel('Amplitude', 'fontsize', 16);
set(gca, 'fontsize', 16);

subplot(312);
plot(t, theta, 'linewidth', 2);
grid
set(gca, 'fontsize', 16);
ylabel('IP (rad)', 'fontsize', 16);

subplot(313);
plot(t, fhat, 'r');
hold on;
plot(t, f, 'b--', 'linewidth', 2);
grid
xlabel('time(s)', 'fontsize', 16);
ylabel('IF (Hz)', 'fontsize', 16);
set(gca, 'fontsize', 16);
a = axis;
a(3) = -20;
a(4) = 75;
axis(a);

ff = fs*(0:N-1)/N;
figure
hold on;
plot(ff, 20*log10(abs(fft(xa))));
plot(ff, 20*log10(abs(fft(cos(ph)))), 'r');
grid

figure
polar(theta(I1), abs(xa(I1)));
hold on
polar(theta(I1(1:1:end)), abs(xa(I1(1:1:end))), 'ro');
set(gca, 'fontsize', 16);
axis 'tight'
line([0 abs(xa(I1(ind1)))*cos(theta(I1(ind1)))], [0 abs(xa(I1(ind1)))*sin(theta(I1(ind1)))], 'linewidth', 3, 'color', 'k');
line([0 abs(xa(I1(ind1+1)))*cos(theta(I1(ind1+1)))], [0 abs(xa(I1(ind1+1)))*sin(theta(I1(ind1+1)))], 'linestyle', '--','linewidth', 3, 'color', .5*ones(1,3));
title('Phasor Representation in Low-amplitude', 'fontsize', 12)

figure
polar(theta(I2), abs(xa(I2)));
hold on
polar(theta(I2(1:1:end)), abs(xa(I2(1:1:end))), 'ro');
set(gca, 'fontsize', 16);
axis 'tight'
line([0 abs(xa(I2(ind2)))*cos(theta(I2(ind2)))], [0 abs(xa(I2(ind2)))*sin(theta(I2(ind2)))], 'linewidth', 3, 'color', 'k');
line([0 abs(xa(I2(ind2+1)))*cos(theta(I2(ind2+1)))], [0 abs(xa(I2(ind2+1)))*sin(theta(I2(ind2+1)))], 'linestyle', '--', 'linewidth', 3, 'color', .5*ones(1,3));
title('Phasor Representation in High-amplitude', 'fontsize', 12)
