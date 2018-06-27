function b = BaseLine1(x,L,approach)
%
% b = BaseLine1(x,L,approach),
% Baseline wander extraction from biomedical recordings, using a single 
% stage of median or moving average filtering.
%
% inputs:
% x: vector or matrix of noisy data (channels x samples)
% L: averaging window length (in samples)
% approach:
%   'md': median filtering
%   'mn': moving average
%
% output:
% b: vector or matrix of baseline wanders (channels x samples)
%
%
% Open Source ECG Toolbox, version 1.0, November 2006
% Released under the GNU General Public License
% Copyright (C) 2006  Reza Sameni
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

N = size(x,2);
b = zeros(size(x));
flen = floor(L/2);

if (strcmp(approach,'mn'))      % moving average filter
    for j = 1:N,
        index = max(j-flen,1):min(j+flen,N);
        b(:,j) = mean(x(:,index),2);
    end
elseif (strcmp(approach,'md'))  % median filter
    for j = 1:N,
        index = max(j-flen,1):min(j+flen,N);
        b(:,j) = median(x(:,index),2);
    end
end