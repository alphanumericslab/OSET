function y = BPFilter2(x,fl,fu,N)
%
% y = BPFilter2(x,fl,fu,N),
% Bandpass filter using FFT filtering with variable FFT points
%
% inputs:
% x: vector or matrix of input data (channels x samples)
% fl: normalized lower frequency
% fu: normalized upper frequency
% N: number of FFT points
%
% output:
% y: vector or matrix of filtered data (channels x samples)
%
% Note:
% - fl and fu are the lower and upper frequency ranges of the bandpass filter
% normalized by the sampling frequency
% - The filter does not perform any windowing on the data
%
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


S = fft(x,N,2);

k = 1:ceil(fl*N);
if(~isempty(k)),
    S(:,[k N-k(2:end)+2]) = 0;
end

k = floor(fu*N):ceil(N/2)+1;
if(~isempty(k)),
    S(:,[k N-k(1:end-1)+2]) = 0;
end

y = real(ifft(S,size(x,2),2));