function [peaks,mn,r] = PeakDetection3(x,fs,h,th,fmax)
%
% peaks = PeakDetection3(x,fs,h,th,fmax),
% R-peak detector based on a matched filter
%
% inputs:
% x: vector of input data
% fs: sampling rate
% h: template waveform
% th: detection threshold
% fmax: maximum expected frequency of the R-peaks
%
% output:
% peaks: vector of R-peak impulse train
%
%
% Open Source ECG Toolbox, version 2.0, March 2008
% Released under the GNU General Public License
% Copyright (C) 2008  Reza Sameni
% Sharif University of Technology, Tehran, Iran -- GIPSA-Lab, INPG, Grenoble, France
% reza.sameni@gmail.com

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.

N = length(x);
L = length(h);

h = h(end:-1:1);

w = floor(L/2);

r = filter(h,1,[x zeros(1,w-1)]);

r = r(w:N+w-1);

r(r<th*max(r)) = 0;

peaks = zeros(1,length(x));

if(nargout==1)
    wlen2 = round(fs/fmax);
    I = find(r>0);

    for i = 1:length(I)
        ind = max(1,I(i)-wlen2):min(N,I(i)+wlen2);
        if (max(r(ind))==r(I(i)))
            peaks(I(i)) = 1;
        end
    end
else
    wlen2 = round(fs/fmax);
    I = find(r>0);

    seg = zeros(length(I),L);
    for i = 1:length(I)
        ind = max(1,I(i)-wlen2):min(N,I(i)+wlen2);
        if (max(r(ind))==r(I(i)))
            peaks(I(i)) = 1;
            sg = x(max(1,I(i)-w+1):min(N,I(i)+w));
            seg(i,:) = [sg zeros(1,L-length(sg))];
        end
    end

    mn0 = mean(seg,1);

    % Robust weighted averaging
    noise = seg - mn0(ones(size(seg,1),1),:);
    vr = var(noise,[],2);
    sm = sum(1./vr);
    weight = 1./(vr*sm);
    mn = weight'*seg;

    mn = sqrt(sum(h.^2))*mn /sqrt(sum(mn.^2));
end