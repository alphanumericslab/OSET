function [T0,T1] = SynchPhaseTimes2(peaks)
%
% [T0,T1] = SynchPhaseTimes2(peaks)
% Calculation of synchronous times instants from beat to beat given a set of
% R-peaks, required for Pi-CA algorithm
%
% input:
% peaks: vector of R-peak pulse train
%
% outputs:
% T0: first (or reference) time instant vector
% T1: second time vector having synchronous phases with T0
%
% Open Source ECG Toolbox, version 2.0, March 2009
% Released under the GNU General Public License
% Copyright (C) 2009  Reza Sameni
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

I = find(peaks);
D = I(2:end) - I(1:end-1);

if (length(I)<3)
    T0 = [];
    T1 = [];
else
    start = I(1);
    stop = I(end-1);

    T1 = zeros(1,stop-start+1);
    k = 1;
    for t = start:stop
        T1(t-start+1) = I(k+1) + round( (t-I(k))*D(k+1)/D(k) );
        if(t>=I(k+1))
            k = k + 1;
        end
    end
    T0 = start:stop;
end