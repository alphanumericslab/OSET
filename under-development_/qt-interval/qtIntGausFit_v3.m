function qtInt=qtIntGausFit_v3(ecg, fs, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% estimate the Q-T interval by fitting Gaussian functions.
% 
% Syntax:
% qtInt=qtParamsGausFit(ecg, fs)
% qtInt=qtParamsGausFit(ecg, fs, rpeaks)
% qtInt=qtParamsGausFit(ecg, fs, rpeaks, betaQ, betaT)
% 
% 
% Inputs:
% ecg: vector of single ECG data.
% fs: sampling frequency.
% rpeaks (or hr): If vector, it is index vector of R peak positions with the
%   same lenth of ecg columns. rpeaks has 1 or true values at R peaks and 0
%   or flase values elsewhere. If scalar, it is the average heart rate in Hz,
%   used in R peak detection.
% betaQ & betaT: (scalar) constant coefficient in (wave width)=2*beta*(Gaussian width), and 
%   wave offset/onset = Gaussian center) +/- beta*(Gaussian width).
% 
% Note: each of rpeaks, betaQ and betaT can be left empty ([]), so the
%   default values will be used.
% 
% Output:
% qtInt: Q-T intervals lendth.
% 
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2021  Davood Fattahi
% fattahi.d@gmail.com
% 
% 
% 
% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% initial values
options=struct; rPeaks=[]; beta.q=[]; beta.t=[];

%% factor of onset/offset detection
if nargin >=3
    rPeaks = varargin{1};
end
if nargin >=4
    beta.q = varargin{2};
end
if nargin ==5
    beta.t = varargin{3};
end
%% R peaks
if isempty(rPeaks)
    rPeaks = PeakDetection20(ecg,70/60/fs,.6);
end
if length(rPeaks)==length(ecg)
    rPeaks = find(rPeaks);
end


%% segments of interest
[~, Q,~, ~, T] = ecgWavesSoI_k(ecg, rPeaks, 1);
soi.q = Q(:,[1,3])/fs;
soi.t = T(:,[1,3])/fs;


%% initial points for gauss fitting optimization
p0=[];
% p0.q=[0; .02;  -.05];
% p0.t=[0; .05; .25];

%% lower and upper bounds for gauss fitting optimization
lb = [];
ub = [];


%% ooptimization options
options = struct('SpecifyObjectiveGradient',true);

[~, ~, ~, ~, qtInt] = qtParamsGausFit_cl(ecg, fs, rPeaks, p0, soi, lb, ub, beta, options);



