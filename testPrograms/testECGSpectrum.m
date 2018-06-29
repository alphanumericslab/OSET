%
% Test ECG Spectrum plot
%
% Dependencies: The baseline wander toolbox of the Open Source ECG Toolbox
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

clc
clear
close all;

% load('SampleECG2.mat'); data = data(:,2)';
load('SampleECG2.mat'); data = data(:,10)';

fs = 1000;
t = (0:length(data)-1)/fs;

% b = Median(data,length(data),round(fs*.3),round(fs*.6))';
b = LPFilter(data,5/fs);

x = data - b;

% % % [x1,x2] = KFNotch(x,50,fs,1e-6,var(x));

figure;
psd(x,500,fs);
% % % hold on
% % % psd(x2,500,fs);

figure;
plot(t,x,'b');
hold on
% % % plot(t,x2,'r');
grid;
xlabel('time(s)');

% fid = fopen('ECG2.txt','wt');
% fprintf(fid,'%f\n',x);
% fclose(fid);