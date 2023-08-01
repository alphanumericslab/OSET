function y = HRSmoother(x, ExcludedHR, wlen, varargin)
%
% HR = HRSmoother(x, ExcludedHR, wlen, RefLevel),
% Heart rate smoothing
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
    med = varargin{1};
else
    med = median(x);
end

N = length(x);
I = find(abs(x-med)>=ExcludedHR);
y = x;
for j = 1:length(I)
    index = max(I(j)-wlen,1):min(I(j)+wlen,N);
    y(I(j)) = median(x(index));
end
