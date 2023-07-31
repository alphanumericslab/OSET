% *************************************************************************
% * Test program for Instantaneous Phase (IP) estimation using the        *
% * traditional analytic representation approach through FIR filtering    * 
% * and Hilbert Transform. Refer to the User Guide for further details.   *
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
ndft = 100;                                            
WS = 1;     
fs = 173.61;
max_f = 30;
tt = (0:n-1)/fs;
ff = 1:max_f;
[phase, inst_freq, amp] = Phase_Ex_Trad(sig, fs, [], max_f, []);

%%-visualization-%%
figure
subplot(2,2, [1, 3])
mesh(tt, ff, phase)
title('Unwrapped Instantaneous Phase Sequences')
xlabel('Time (Sec)')
ylabel('frequency channels')
subplot(222)
plot(tt, phase(8, :))
grid on
title('Unwrapped Instantaneous Phase of 8Hz')
xlabel('Time (Sec)')
ylabel('Radians')
subplot(224)
plot(tt, fs*[diff(phase(8, :)), 0]/2/pi)
grid on
title('Instantaneous Frequency of 8Hz')
xlabel('Time (Sec)')
ylabel('Hz')

figure
subplot(121)
mesh(tt, ff, inst_freq)
title('Instantaneous Frequency')
xlabel('Time (Sec)')
ylabel('Frequency (Hz)')
subplot(122)
mesh(tt, ff, amp)
title('Instantaneous Amplitudes')
xlabel('Time (Sec)')
ylabel('Frequency (Hz)')
