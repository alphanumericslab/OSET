% *************************************************************************
% * This code calculates and illustrates the probability of signal        *
% * detectability versus required SNR by using the probability of false   *
% * alarm (represented in [2]).                                           *
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

pfa = [1e-1 1e-2 1e-3 1e-4 1e-5 1e-6 1e-8];      % false-alarm probabilities

rocpfa(pfa,'SignalType', 'NonfluctuatingNoncoherent', 'MinSNR',-5, 'MaxSNR', 20);
ylabel('Probability of Detection (p_d)','fontsize', 16);
xlabel('Signal-to-Noise Ratio (dB)','fontsize', 16);
set(gca, 'fontsize', 16);

