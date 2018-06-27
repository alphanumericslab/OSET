function [ECG teta]= SingleChannelECGGenerator(teta,teta0,alphai,bi,tetai);
%
% [ECG teta]= DipoleGenerator2(teta,teta0,alphai,bi,tetai)
% Synthetic cardiac dipole generator using the 'direct form' of the
% dipole equations. Identical with the 'DipoleGenerator' function.
% Refer to references of the toolbox for further details.
%
% inputs:
% teta: ECG phase calculated from real ECG
% teta0: desired phase shift of the synthetic ECG
% alphai: structure contaning the amplitudes of Gaussian functions used for
%       modeling the x, y, and z coordinates of the cardiac dipole
% bi: structure contaning the widths of Gaussian functions used for
%       modeling the x, y, and z coordinates of the cardiac dipole
% tetai: structure contaning the phase of Gaussian functions used for
%       modeling the x, y, and z coordinates of the cardiac dipole
%
% output:
% ECG: synthetic ECG signal
% teta: vector containing the shifted ECG phase
%
%
% Open Source ECG Toolbox, March 2007
% Released under the GNU General Public License
% Copyright (C) 2007  Reza Sameni
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

N = length(teta);

teta = mod(teta + teta0 + pi, 2*pi) - pi;

dtetai = mod(teta(ones(length(tetai),1),:)' - tetai(ones(1,N),:) + pi , 2*pi) - pi;

ECG = [sum(alphai(ones(1,N),:) .* exp(-dtetai .^2 ./ (2*bi(ones(1,N),:) .^ 2)),2)]';
