function y = BPFilter6(x,fc,bw, factor)
%
% y = BPFilter6(x,fc,bw,order),
% Bandpass filter using an equiripple FIR filter
%
% inputs:
% x: vector or matrix of input data (channels x samples)
% fc: normalized center frequency
% bw: normalized bandwidth
% factor: controls the sharpness of the filter and its convergence (0 < factor < 1)
%
% output:
% y: vector or matrix of filtered data (channels x samples)
%
% Note:
% - fc and bw are the center frequency and bandwidth of the bandpass filter
% normalized by the sampling frequency
%
% Open Source ECG Toolbox, version 2.2, July 2014
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

% the equiripple FIR filter
[nn,fo,mo,w] = firpmord( [factor*bw/2 bw/2/factor], [1 0], [0.001 0.1], 1 );
h = firpm(nn,fo,mo,w);

% the exponentials used for frequency shifting
w = exp(-1j*2*pi*fc*n);
v = exp(1j*2*pi*fc*n);
x1 = x.*w(ones(M,1),:);
y1 = zeros(M,N);
for i = 1:M
    y1(i,:) = filter(h,1,x1(i,:));
end
z1 = y1.*v(ones(M,1),:);

x2 = x.*v(ones(M,1),:);
y2 = zeros(M,N);
for i = 1:M
    y2(i,:) = filter(h,1,x2(i,:));
end
z2 = y2.*w(ones(M,1),:);

y = real(z1 + z2);