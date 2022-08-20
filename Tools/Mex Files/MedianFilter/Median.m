%  y = Median(x,N,L1,L2),
% 
%  Single or Double stage median filter
%  inputs:
% 	 x: input vector. x should be in 'double' format.
% 	 N: input signal length (length(x))
% 	 L1: first stage median window length (in samples)
% 	 L2: second stage median window length (in samples)
% 
%  output:
% 	 y: output column vector
% 
%  comment: consider using 'TrimmedFilter' instead, which has more options
%  and has been optimized for speed
% 
% 	 Open Source ECG Toolbox, version 2.0, March 2008
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