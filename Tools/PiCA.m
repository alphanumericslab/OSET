function [y,W,A] = PiCA(x,peaks,varargin)
%
% [y,W,A] = PiCA(x,peaks1,peaks2)
% Pseudo-Periodic Component Analysis
%
% input:
% x: input data array (channels x samples)
% peaks1: vector of first signal peaks impulse train (R-wave locations for ECG signals)
% peaks2 (optional): vector of second signal peaks impulse train (R-wave locations for ECG signals)
%
% outputs:
% y: pseudo-periodic components ranked in order of resemblance with the
% first to second (if available) desired signals.
% W: the decomposition matrix
% A: the mixing matrix
%
% comment: if the signal peaks2 is not available the components are ranked
% from the most similar to the least similar to the signal with
% corresponding to peaks1.
% 
% Reference:
% R. Sameni, C. Jutten, and M. B. Shamsollahi. Multichannel electrocardiogram
%  decomposition using periodic component analysis. IEEE Transactions on Biomedical Engineering,
%  55(8):1935-1940, Aug. 2008.
%
% Open Source ECG Toolbox, version 2.0, October 2008
% Released under the GNU General Public License
% Copyright (C) 2008  Reza Sameni
% Sharif University of Technology, Tehran, Iran -- LIS-INPG, Grenoble, France
% reza.sameni@gmail.com

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.

if(nargin>2),
    peaks2 = varargin{1};
    flag = 1;
else
    flag = 0;
end

% PM time calculation
[T0,T1] = SynchPhaseTimes2(peaks);

A = x(:,T0)*x(:,T1)';
B = x(:,T0)*x(:,T0)';
A = (A+A')/2;
B = (B+B')/2;

if flag == 0,
    [V,D] = eig(A,B,'chol');
elseif flag ==1,
    % PM time calculation for the second peaks sequence
    [T0,T1] = SynchPhaseTimes2(peaks2);
    AA = x(:,T0)*x(:,T1)';
    AA = (AA+AA')/2;
    [V,D] = eig(A-AA,B,'chol');
end

d = diag(D);
[YY,I] = sort(d);
I = I(end:-1:1);

W = V(:,I)';
A = pinv(W);

y = W*x;