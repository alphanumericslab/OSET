%  y = TrimmedFilter(x,type,L,k,h),
%  y = TrimmedFilter(x,type,L,k,l,h),
% 
%  A moving average, median, or trimmed mean filter.
% 
%  inputs:
% 	 x: input vector. x should be in 'double' format.
%    type:
%       'mean': a moving average
%       'median': a moving median
%       'trmean': trimmed-mean filter
%       'atrmean': asymmetric trimmed-mean filter
%       'wmedian': weighted median filter
% 	 L: moving window length (in samples)
% 	 k: for a trimmed filter, k is the number of samples excluded from
%       both sides of the probability density function in the window of
%       length L, before calculating the moving average. By default k = 0
%       which reduces to a moving average filter. For an asymmetric
%       trimmed-mean filter k is the number of samples excluded from the
%       beginning (lower side) of the probability density function.
%    l: for an asymmetric trimmed-mean filter, l is the number of samples
%       excluded from the end (upper side) of the probability density
%       function.
%    h: weighting window of length L, used for 'wmedian'. If k or l are 
%       non-zero, the probability density function of each window is
%       truncated by k and l samples from either side before calculating
%       the weighted median.
% 
%  output:
% 	 y: output column vector
% 
% 
% 	 Open Source ECG Toolbox, version 2.0, December 2008
% 	 Released under the GNU General Public License
% 	 Copyright (C) 2008  Reza Sameni
% 	 Sharif University of Technology, Tehran, Iran -- GIPSA-Lab, INPG, Grenoble, France
% 	 reza.sameni@gmail.com
% 	 
% 	 This program is free software; you can redistribute it and/or modify it
% 	 under the terms of the GNU General Public License as published by the
% 	 Free Software Foundation; either version 2 of the License, or (at your
% 	 option) any later version.
% 	 This program is distributed in the hope that it will be useful, but
% 	 WITHOUT ANY WARRANTY; without even the implied warranty of
% 	 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% 	 Public License for more details.