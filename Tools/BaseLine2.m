function [x2] = BaseLine2(x,L1,L2,approach)
%
% b = BaseLine2(x,L1,L2,approach),
% Baseline wander extraction from biomedical recordings, using two stages 
% of median or moving average filtering.
%
% inputs:
% x: vector or matrix of noisy data (channels x samples)
% L1: first stage averaging window length (in samples)
% L2: second stage averaging window length (in samples)
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
x1 = zeros(size(x));
x2 = zeros(size(x));

flen1 = floor(L1/2);
flen2 = floor(L2/2);

if (strcmp(approach,'mn'))
    for j = 1:N,
        index = max(j-flen1,1):min(j+flen1,N);
        x1(:,j) = mean(x(:,index),2);
    end

    for j = 1:N,
        index = max(j-flen2,1):min(j+flen2,N);
        x2(:,j) = mean(x1(:,index),2);
    end
elseif (strcmp(approach,'md'))
    for j = 1:N,
        index = max(j-flen1,1):min(j+flen1,N);
        x1(:,j) = median(x(:,index),2);
    end

    for j = 1:N,
        index = max(j-flen2,1):min(j+flen2,N);
        x2(:,j) = median(x1(:,index),2);
    end
end