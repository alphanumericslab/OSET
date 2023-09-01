%
% Test program for generating single-channel ECGs with additive colored
% noise
%
% Dependencies: The synthetic ECG generator and noise generator package of
%   the Open Source ECG Toolbox
%
% Open Source ECG Toolbox, version 1.0, November 2006
% Released under the GNU General Public License
% Copyright (C) 2006  Reza Sameni
% Sharif University of Technology, Tehran, Iran -- LIS-INPG, Grenoble, France
% reza.sameni@gmail.com

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

%//////////////////////////////////////////////////////////////////////////
clc
clear
close all;

N = 3000;
fs = 1000;
t = 0:N-1/fs;

f = 1;                                          % approximate R-peak frequency

peaks = zeros(1,N);
I = 10:round(fs/f):N;
peaks(I) = 1;

[phase, phasepos] = phase_calculator(peaks);     % phase calculation

teta = pi;                                       % desired phase shift
pphase = phase_shifter(phase,teta);             % phase shifting


% L = length(OptimumParams)/3;% number of Gaussian kernels
% alphai = [.2 -.2 1 -.23 .3];
% bi = [.2 .2 .2 .2 .3];
% tetai = pi*[-60 -15 0 15 90]/180;

alphai = [.1 -.2 1 -.3 .15];
bi = [.3 .2 .2 .2 .5];
tetai = pi*[-90 -10 0 15 100]/180;

% teta0 = pi/2;
teta0 = 0;
[ECG, teta]= ecg_gen_gmm(pphase,teta0,alphai,bi,tetai);

%//////////////////////////////////////////////////////////////////////////
% data plotting
figure;
plot(t,ECG,'r','linewidth',2);
set(gca,'box','on','fontsize',16);
% grid
xlabel('time(s)','fontsize',16);
ylabel('Amplitude(mV)','fontsize',16);
title('Synthetic ECG','fontsize',16);
