function [skw, m1, s] = skew(x)
% [skw, mn, sd] = skew(x)
% Skewness of a matrix x over its second dimension
% input:
%   x: N channels times T samples data matrix
% outputs:
%   skw: skewness vector (N times 1)
%   mn: mean vector (N times 1)
%   sd: standard deviation vector (N times 1)
%
% Open Source ECG Toolbox, version 3.14, September 2019
% Released under the GNU General Public License
% Copyright (C) 2019  Reza Sameni
% reza.sameni@gmail.com

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.

m1 = mean(x, 2);
s = std(x, [], 2);
m3 = mean(x.^3, 2);

skw = (m3 - 3*m1.*s.^2 - m1.^3)./s.^3;
