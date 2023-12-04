% *************************************************************************
% * Test program for Instantaneous Phase estimation using the Transfer-   *
% * Function Perturbation Phase estimation method (TFP) [2] through       *
% * analytic representation, IIR filters and forward-backward filtering.  * 
% * Refer to the User Guide for further details.                          *
% *************************************************************************
% 
% Dependencies: The Cerebral Signal Phase Analysis Toolbox of Open Source 
%               Electrophysiological Toolbox
%
% This program is provided by ESMAEIL SERAJ (esmaeil.seraj09@gmail.com). 
% 
% Open Source Electrophysiological Toolbox, version 3.1, 2014
% Released under the GNU General Public License
% Copyright (C) 2012  Reza Sameni
% Shiraz University, Shiraz, Iran
% reza.sameni@gmail.com
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

%%-loading the test signals-%%
sig = load('Z001.txt')';          % Ongoing EEG: change *fs* to 173.61Hz
[m, n] = size(sig);

%%-phase calculation-%%
fs = 173.61;                      % sampling frequency
pertnum = 20;
tt = (0:n-1)/fs;
f0 = 10;                % EEG alpha rhythms
bw_base = 4;            % EEG alpha rhythms
bw_base_dev = 0.1;      % bandwidth deviation
f0_dev = 1e-6;          % center frequency deviation
dither_std = 1e-4;      % dither level

% [phase_avg, freq_avg, amp_avg] = Phase_Ex_TFP(sig, fs);
[phase_avg, freq_avg, amp_avg] = Phase_Ex_TFP(sig, fs, f0, bw_base, [], [], dither_std, []);

%%-visualization-%%
figure
subplot(311)
plot(tt, phase_avg, 'linewidth', 1.5)
xlabel('time (sec)'); ylabel('inst. phase'); axis tight; grid on
title('Instantaneous Parameter Estimation Using TFP Method')
subplot(312)
plot(tt, freq_avg, 'linewidth', 1.5)
xlabel('time (sec)'); ylabel('inst. frequency'); axis tight; grid on
subplot(313)
plot(tt, amp_avg, 'linewidth', 1.5)
xlabel('time (sec)'); ylabel('inst. amplitude'); axis tight; grid on
%}