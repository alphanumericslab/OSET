function data = RemoveChannels(data,I);
%
% y = RemoveChannels(x,I);
% Eliminate unwanted channels
%
% inputs:
% x: matrix of input data. Smaller dimension is considered as channels
% I: vector of unwanted channels
%
% output:
% y: matrix of wanted channels, rearranged in the form of: channels x samples
%
% Open Source ECG Toolbox, version 2.0, May 2009
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

if(size(data,1)>size(data,2))
    data = data';
end

% this is definitely not the best way to solve this problem!
I_ = [];
for i = 1:size(data,1),
    if(isempty(find(I==i)))
        I_ = [I_ i];
    end
end

data = data(I_,:);