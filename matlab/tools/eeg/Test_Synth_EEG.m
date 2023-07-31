% *************************************************************************
% * Test program for generating synthetic EEG signal using Autoregressive *
% * (AR) model and innovation filter. Could be used in cross-validations  *
% * with real EEG signal. Refer to the User Guide for further details.    *
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

%---loading signal and initialization---%
sig = load('EEG.mat');
sig = sig.sig8;
fs = 160;                    % sampling frequency (Hz)
duration = 32;               % total signal duration (sec)
win = 4;                     % temporal window length to have a stationary signal (Sec)
AR_ord = 20;                 % order of AR model

eeg_synth = Synth_EEG(sig, fs, [], [], AR_ord);

%---visualization---%
figure;
subplot(211)
plot(sig)
legend('Real EEG')
grid on
subplot(212)
plot(eeg_synth)
legend('Synthetic EEG')
grid on

figure;
subplot(211)
pwelch(sig)
legend('Real EEG')
grid on
subplot(212)
pwelch(eeg_synth)
legend('Synthetic EEG')
grid on
