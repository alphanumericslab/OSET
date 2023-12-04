function y = BPFilter5(x,fc,bw,order)
%
% y = BPFilter5(x,fc,bw,order),
% Zero-phase bandpass filter using CIC filter
%
% inputs:
% x: vector or matrix of input data (channels x samples)
% fc: normalized center frequency
% bw: normalized bandwidth
% order: CIC filter order
%
% output:
% y: vector or matrix of filtered data (channels x samples)
%
% Note:
% - fc and bw are the center frequency and bandwidth of the bandpass filter
% normalized by the sampling frequency
% - the filter performs forward-reverse filtering successively. It has 
% zero-phase for even CIC filter orders and a phase-lag equal to a single
% stage CIC filter for odd CIC filter orders
%
% Open Source ECG Toolbox, version 2.0, July 2010
% Released under the GNU General Public License
% Copyright (C) 2010  Reza Sameni
% Shiraz University, Shiraz, Iran
% reza.sameni@gmail.com

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.

% general parameters
M = size(x,1);
N = size(x,2);
n = 0:N-1;

% simple lowpass moving average filter used as template
L = round(2/bw);
h = ones(1,L)/L;

% the exponentials used for frequency shifting
w = exp(-1j*2*pi*fc*n);
v = exp(1j*2*pi*fc*n);

% the real part of the filter
x1 = x.*w(ones(M,1),:);
y1 = x1;
for k = 1:order,
    for i = 1:M,
        y1(i,:) = filter(h,1,y1(i,:));
    end
    y1 = y1(:,end:-1:1);
end
if(mod(order,2)==1)
    y1 = y1(:,end:-1:1);
end    
z1 = y1.*v(ones(M,1),:);

% the imaginary part of the filter
x2 = x.*v(ones(M,1),:);
y2 = x2;
for k = 1:order,
    for i = 1:M,
        y2(i,:) = filter(h,1,y2(i,:));
    end
    y2 = y2(:,end:-1:1);
end
if(mod(order,2)==1)
    y2 = y2(:,end:-1:1);
end    
z2 = y2.*w(ones(M,1),:);

y = real(z1 + z2);