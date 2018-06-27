function y = HRSmoother2(x, FilterHR, wlen, varargin)
%
% HR = HRSmoother2(x, FilterHR, wlen, mode),
% Conditional heart rate smoothing
%
% Open Source ECG Toolbox, version 2.1, September 2010
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

if(nargin>3 && ~isempty(varargin{1}))
    mode = varargin{1};
else
    mode = 'normal';
end

N = length(x);
y = zeros(1,N);
for j = 1:N,
    index = max(j-wlen,1):min(j+wlen,N);
    y(j) = median(x(index));
end
I = find(abs(y-x) < FilterHR);
J = find(abs(y-x) >= FilterHR);

if(strcmp(mode,'normal'))
    y(I) = x(I);
elseif(strcmp(mode,'nan'))
    y = x;
    y(J) = nan;
end

