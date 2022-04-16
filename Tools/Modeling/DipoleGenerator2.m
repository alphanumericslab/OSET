function [DIP teta]= DipoleGenerator2(N,fs,f,alphai,bi,tetai,teta0);
%
% [DIP teta]= DipoleGenerator2(N,fs,f,alphai,bi,tetai,teta0)
% Synthetic cardiac dipole generator using the 'direct form' of the
% dipole equations. Identical with the 'DipoleGenerator' function.
% Refer to references of the toolbox for further details.
%
% inputs:
% N: signal length
% fs: sampling rate
% f: heart rate (Hz)
% alphai: structure contaning the amplitudes of Gaussian functions used for
%       modeling the x, y, and z coordinates of the cardiac dipole
% bi: structure contaning the widths of Gaussian functions used for
%       modeling the x, y, and z coordinates of the cardiac dipole
% tetai: structure contaning the phase of Gaussian functions used for
%       modeling the x, y, and z coordinates of the cardiac dipole
% teta0: initial phase of the synthetic dipole
%
% output:
% DIP: structure contaning the x, y, and z coordinates of the cardiac dipole
% teta: vector containing the dipole phase
%
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

w = 2*pi*f;
dt = 1/fs;

teta = cumsum([teta0+pi,w*dt*ones(1,N-1)]);
teta = mod(teta,2*pi)-pi;

dtetaix = mod(teta(ones(length(tetai.x),1),:)' - tetai.x(ones(1,N),:) + pi , 2*pi) - pi;
dtetaiy = mod(teta(ones(length(tetai.y),1),:)' - tetai.y(ones(1,N),:) + pi , 2*pi) - pi;
dtetaiz = mod(teta(ones(length(tetai.z),1),:)' - tetai.z(ones(1,N),:) + pi , 2*pi) - pi;

X = sum(alphai.x(ones(1,N),:) .* exp(-dtetaix .^2 ./ (2*bi.x(ones(1,N),:) .^ 2)),2);
Y = sum(alphai.y(ones(1,N),:) .* exp(-dtetaiy .^2 ./ (2*bi.y(ones(1,N),:) .^ 2)),2);
Z = sum(alphai.z(ones(1,N),:) .* exp(-dtetaiz .^2 ./ (2*bi.z(ones(1,N),:) .^ 2)),2);

DIP.x = X';
DIP.y = Y';
DIP.z = Z';