% *************************************************************************
% * Test program for calculating popular phase related quantities, i.e.   *
% * Phase Shift (PS), Phase Lock (PL), Phase Reset (PR), Phase            *
% * Difference (PD) and Instantaneous Frequency (IF) in SINGLE or         *
% * MULTI-Channel modes using phase sequences. Refer to the User Guide for*
% * further details.                                                      *
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

%%-loading test phase signals-%%
load('Z001_phase')
load('Z002_phase')

%%-parameter specification and phase feature calculation-%%
Th_PS = 0.02;
[m, n] = size(phase_sig1);
fs = 173.61;
tt = (0:n-1)/fs;

[PR, PS, PL, PDV, PD] = Phase_Features(phase_sig1, fs, Th_PS, phase_sig2);  % Multi-signal 
% [IF, PR, PS, PL] = Phase_Features(phase_sig1, fs, Th_PS);                 % Multi or Single -channel signal

%%-visualization-%% 
% Note that one of the presented figures must be used due to using multi or
%%single channel phase features calculation

% multi-channel
figure
subplot(511)
plot(tt, PD)
title('Phase Difference')
subplot(512)
plot(tt, PDV)
title('Phase Difference Variations')
subplot(513)
plot(tt, PL)
title('Phase Lock Events')
subplot(514)
plot(tt, PS)
title('Phase Shift Events')
subplot(515)
plot(tt, PR)
title('Phase Resseting')
xlabel('Time (Seconds)')

% % Single-channel
% figure
% subplot(411)
% plot(tt, IF)
% title('Instantaneous Frequency')
% subplot(412)
% plot(tt, PL)
% title('Phase Lock Events')
% subplot(413)
% plot(tt, PS)
% title('Phase Shift Events')
% subplot(414)
% plot(tt, PR)
% title('Phase Resseting')
% xlabel('Time (Seconds)')
