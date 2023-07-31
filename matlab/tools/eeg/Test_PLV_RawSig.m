% *************************************************************************
% * Test program for Calculating Phase Locking Value (PLV) matrix         *
% * (Pairwise) using raw signals. Refer to the User Guide for further     *
% * details.                                                              *
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

%%-loading the signals-%%
sig1 = load('Z001.txt')';       % Ongoing EEG: change *fs* to 173.61Hz
sig2 = load('Z002.txt')';       % Ongoing EEG: change *fs* to 173.61Hz
sig3 = load('Z003.txt')';       % Ongoing EEG: change *fs* to 173.61Hz
sig1 = [sig1; sig2; sig3];
[m, n] = size(sig1);

%%-parameter specification and phase feature calculation-%%
Method = 'ZPPP';                % available methods: 'ZPPP' and 'Trad'
WS = 1;                         % stop-band frequency
fs = 173.61;                    % sampling frequency
f0 = 8;                         % frequency of interest
ndft = 100;                     % number of frequency bins
pertnum = 100;                  % attempts to perturb filter's response

PLV = PLV_RawSig(Method, sig1, fs, f0, WS, [], []);
% PLV = PLV_RawSig(Method, sig2, fs, f0, [], [], [], sig3);

%%-visualization-%%
figure
imagesc(PLV)
xlabel('Channels')
ylabel('Channels')
title('Claculated Pairwise PLV Between Signals')
colorbar
