function y = BPFilter3(x,fc,bw,wname)
%
% y = BPFilter3(x,fc,bw,@wname),
% Bandpass filter using windowed FFT filtering
%
% inputs:
% x: vector or matrix of input data (channels x samples)
% fc: normalized center frequency
% bw: normalized bandwidth
% @wname: any valid window function name, see help for Matlab's WINDOW
% function
%
% output:
% y: vector or matrix of filtered data (channels x samples)
%
% Note:
% - fc and bw are the center frequency and bandwidth of the bandpass filter
% normalized by the sampling frequency
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

N = size(x,2);

S = fft(x,N,2);

k = ceil((fc-bw/2)*N):floor((fc+bw/2)*N);
k(k<1) = [];
k(k>N) = [];

w = window(wname,length(k));
w = w(:)';

S(:,k) = (ones(size(S,1),1)*w).*S(:,k);
S(:,N-k+2) = (ones(size(S,1),1)*w).*S(:,N-k+2);

S1 = zeros(size(S));
S1(:,k) = S(:,k);
S1(:,N-k+2) = S(:,N-k+2);

y = real(ifft(S1,N,2));