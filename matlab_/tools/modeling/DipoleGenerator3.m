function [DIP teta]= DipoleGenerator3(HR,fs,alphai,bi,tetai,teta0)
%
% [DIP teta]= DipoleGenerator3(HR,fs,alphai,bi,tetai,teta0)
% Synthetic cardiac dipole generator using the 'differential form' of the
% dipole equations. Refer to references of the toolbox for further details.
%
% inputs:
% HR: heart rate vector in beats per minute (BPM)
% fs: sampling rate
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
% Open Source Electrophysiological Toolbox, version 2.1, May 2012
% Released under the GNU General Public License
% Copyright (C) 2012  Reza Sameni
% Shiraz University, Shiraz, Iran
% reza.sameni@gmail.com

dt = 1/fs;

X = [];
Y = [];
Z = [];
teta = teta0;
i = 1;
while(1)
    w = 2*pi*HR(i)/60;
    
    dtetaix = mod(teta(end) - tetai.x + pi , 2*pi) - pi;
    dtetaiy = mod(teta(end) - tetai.y + pi , 2*pi) - pi;
    dtetaiz = mod(teta(end) - tetai.z + pi , 2*pi) - pi;
    
    if(i==1),
        X = sum(alphai.x .* exp(-dtetaix .^2 ./ (2*bi.x .^ 2)));
        Y = sum(alphai.y .* exp(-dtetaiy .^2 ./ (2*bi.y .^ 2)));
        Z = sum(alphai.z .* exp(-dtetaiz .^2 ./ (2*bi.z .^ 2)));
    end
    
    X = [X , X(end) - dt*sum(w*alphai.x ./ (bi.x .^ 2) .* dtetaix .* exp(-dtetaix .^2 ./ (2* bi.x .^ 2)))];   % x state variable
    Y = [Y , Y(end) - dt*sum(w*alphai.y ./ (bi.y .^ 2) .* dtetaiy .* exp(-dtetaiy .^2 ./ (2* bi.y .^ 2)))];   % y state variable
    Z = [Z , Z(end) - dt*sum(w*alphai.z ./ (bi.z .^ 2) .* dtetaiz .* exp(-dtetaiz .^2 ./ (2* bi.z .^ 2)))];   % z state variable
    teta = [teta , mod(teta(end) + w*dt + pi , 2*pi) - pi];
    if(abs(teta(end)-teta(end-1))>pi)
        i = i + 1;
    end
    if(i > length(HR))
        break;
    end
    
end

DIP.x = X;
DIP.y = Y;
DIP.z = Z;