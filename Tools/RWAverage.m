function mn = RWAverage(x)
%
% mn = RWAverage(x)
% Robust weighted averaging
%
% inputs:
% x: an N x T matrix containing N ensembles of a noisy event-related signal of length T
%
% output:
% mn: the robust weighted average over the N rows of x
%
% Reference:
% J.M. Leski. Robust weighted averaging of biomedical signals. IEEE Trans. Biomed. Eng., 49(8):796-804, 2002.
%
% Open Source ECG Toolbox, version 2.0, June 2008
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

if(size(x,1)>1)
    mn0 = mean(x,1);
    noise = x - mn0(ones(size(x,1),1),:);
    vr = var(noise,[],2);
    sm = sum(1./vr);
    weight = 1./(vr*sm);
    mn = weight'*x;
else
    mn = x;
end
