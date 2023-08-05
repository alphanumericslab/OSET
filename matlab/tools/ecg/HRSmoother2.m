function y = HRSmoother2(x, FilterHR, wlen, varargin)
%
% HR = HRSmoother2(x, FilterHR, wlen, mode),
% Conditional heart rate smoothing
% replaces outliers defined by threshold FilterHR with the median value
% over window length wlen, or with NaN, depending on the selected mode.
%
% Open Source ECG Toolbox, version 2.1, September 2010
% Released under the GNU General Public License
% Copyright (C) 2010  Reza Sameni
% Shiraz University, Shiraz, Iran
% reza.sameni@gmail.com
%
% Revisions:
% February 2019: code simplification

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.

if(nargin>3 && ~isempty(varargin{1}))
    mode = varargin{1};
else
    mode = 'normal';
end

T = length(x);
x_med = zeros(1,T);
for j = 1 : T
    index = max(j-wlen, 1) : min(j+wlen, T);
    x_med(j) = median(x(index));
end

y = x;
J = find(abs(x - x_med) >= FilterHR);

if(strcmp(mode,'normal'))
    y(J) = x_med(J);
elseif(strcmp(mode,'nan'))
    y(J) = nan;
end

