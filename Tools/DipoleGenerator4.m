function [DIP teta]= DipoleGenerator4(HR,fs,alphai,bi,tetai,teta0)
%
% [DIP teta]= DipoleGenerator4(HR,fs,alphai,bi,tetai,teta0)
% Synthetic cardiac dipole generator. This is a faster version of This is a
% faster version of 'DipoleGenerator4.m'
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
w = 2*pi*HR/60;
samples_per_beat = round(fs./(HR/60));
samples_per_beat(1) =  round((1 - (teta0+pi)/(2*pi))*samples_per_beat(1));
cumlen = cumsum(samples_per_beat);

N = sum(samples_per_beat);
ww = zeros(1,N);
ww(1:cumlen(1)) = repmat(w(1),1,samples_per_beat(1));
for i = 2:length(w)
%     ww(cumlen(i-1)+1:cumlen(i)) = repmat(w(i),1,samples_per_beat(i));
    ww(cumlen(i-1)+1:cumlen(i)) = w(i);
end
ww(cumlen(length(w)):end) = w(end);
    
teta = cumsum([teta0+pi,ww*dt]);
teta = mod(teta,2*pi)-pi;
teta = teta(2:end);

% figure
% plot(ww)

dtetaix = mod(teta(ones(length(tetai.x),1),:)' - tetai.x(ones(1,N),:) + pi , 2*pi) - pi;
dtetaiy = mod(teta(ones(length(tetai.y),1),:)' - tetai.y(ones(1,N),:) + pi , 2*pi) - pi;
dtetaiz = mod(teta(ones(length(tetai.z),1),:)' - tetai.z(ones(1,N),:) + pi , 2*pi) - pi;

X = sum(alphai.x(ones(1,N),:) .* exp(-dtetaix .^2 ./ (2*bi.x(ones(1,N),:) .^ 2)),2);
Y = sum(alphai.y(ones(1,N),:) .* exp(-dtetaiy .^2 ./ (2*bi.y(ones(1,N),:) .^ 2)),2);
Z = sum(alphai.z(ones(1,N),:) .* exp(-dtetaiz .^2 ./ (2*bi.z(ones(1,N),:) .^ 2)),2);

DIP.x = X';
DIP.y = Y';
DIP.z = Z';