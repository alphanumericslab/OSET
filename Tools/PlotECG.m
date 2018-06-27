function PlotECG(data,L,varargin)
%
% PlotECG(data,L,color,fs)
% Plot multichannel signals in several panels and figures
%
% inputs:
% data: the input data matrix
% L: number of panels per figure
% color: the color of the plots: 'b','r','g',etc. (blue by default)
% fs: sampling rate (1 by default)
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

%//////////////////////////////////////////////////////////////////////////
% input arguments
if(nargin>2  && ~isempty(varargin{1})),
    color = varargin{1};
else
    color = 'b';
end

if(nargin>3  && ~isempty(varargin{2})),
    fs = varargin{2};
    lbl = 'time(s)';
else
    fs = 1;
    lbl = 'index';
end

%//////////////////////////////////////////////////////////////////////////
if(size(data,1)>size(data,2))
    data = data';
end

L1 = size(data,1);
L2 = size(data,2);

t = (0:L2-1)/fs;

for i = 1:L1,
    if(mod(i,L)==1 || L==1)
        figure;
    end
    subplot(L,1,mod(i-1,L)+1);
    plot(t,data(i,:),color);
    ylabel(num2str(i));
    grid;
    if(mod(i,L)==0 || L==1)
        xlabel(lbl);
    end
end
