% *************************************************************************
% * Test program for calculating Phase Locking Value (PLV) matrix         *
% * (Pairwise) using phase sequences. Refer to the User Guide for further *
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

%%-loading test phase signals and calculating PLV matrix-%%
load('Z001_phase')
load('Z002_phase')
load('Z003_phase')
load('Z_phase')

% PLV = PLV_PhaseSeq(phase_sig1, phase_sig2, phase_sig3); % phase vectors
PLV = PLV_PhaseSeq(phase_sig);                            % phase matrix

%%-visualization-%%
figure
imagesc(PLV)
xlabel('Channels')
ylabel('Channels')
title('Claculated Pairwise PLV Between Signals')
colorbar
