function [y1,y2,Pbar,Phat,PSmoothed,Kgain] = KFNotch(x,f0,fs,varargin)
%
% [y1,y2,Pbar,Phat,PSmoothed,Kgain] = KFNotch(x,f0,fs,Q,R,gamma),
% Removing Power-Line noise using a linear Kalman filter and smoother
%
% inputs:
% x: vector of noisy signals contaminated with power-line
% f0: the notch frequency
% fs: sampling rate
% Q: covariance of the model errors (optional). By default: 1e-4*max(abs(x))
% R: covariance of the power-line noise (optional). By default: var(x)
%
% output:
% y1: vector of denoised signal using the Kalman filter
% y2: vector of denoised signal using the Kalman smoother
% Pbar: covariance matrix of the a priori error vector of the Kalman filter
% Phat: covariance matrix of the a posteriori error vector of the Kalman filter
% PSmoothed: covariance matrix of the a posteriori error vector of the Kalman smoother
% Kgain: vector of Kalman filter gain
% 
% practical notes:
% 1- For a normalized input signal, I usually take 1e-6<Q<1e-3
% 2- R is the variance of the signal of interest which is definitely less 
%   than var(x). Depending on the amplitude of the power-line noise you can
%   play with it in the range of, say, 1e-6*var(x)<R<0.1*var(x) and you can
%   find its exact value by changing its value and checking the spectrum on
%   notch frequency.
% 3- 0<gamma<1, when you set it to 1 the value of R is fixed to what you give
%   it; but by setting gamma<1 then R will be adapted in time. I usually get
%   reasonable adaptation with 0.7<gamma<.99
% 4- Note that, KFNotch only filters one single frequency and not its harmonics.
%   So if the harmonics are notable you need to apply KFNotch for each harmonic
%   separately.
%
% Reference:
% R. Sameni, "A Kalman Notch Filter for Removing Power-Line Noise from
% Biomedical Signals", Technical Report, GIPSA-LAB, October 15th 2007.
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

%//////////////////////////////////////////////////////////////////////////
% input arguments
if(nargin>3  && ~isempty(varargin{1})),
    Q = varargin{1};
else
    Q = 1e-4*max(abs(x));
end

if(nargin>4 && ~isempty(varargin{2})),
    R = varargin{2};
else
    R = var(x);
end

if(nargin>5  && ~isempty(varargin{3})),
    gamma = varargin{3};
else
    gamma = 1;
end

%//////////////////////////////////////////////////////////////////////////
% Kalman filter parameters
w = 2*pi*f0/fs;

Wmean = 0;
Vmean = 0;
X0 = [x(2) x(1)]';
P0 = 1000*Q;

L = 2;
A = [2*cos(w) -1 ; 1 0];
H = [1 0];
Xminus = X0;
Pminus = P0;
Samples = length(x);

VarWinlen = ceil(fs/10);
mem2 = zeros(VarWinlen,1) + R;
Xhat = zeros(2,Samples);
innovations = zeros(1,Samples);
Phat = zeros(2,2,Samples);
Xbar = zeros(2,Samples);
Pbar = zeros(2,2,Samples);
Kgain = zeros(2,Samples);
B = [1 0]';

%//////////////////////////////////////////////////////////////////////////
% Forward Filtering Stage
for k = 1 : Samples,
    % Store results
    Xbar(:,k) = Xminus';
    Pbar(:,:,k) = Pminus';

    % Measurement update (A posteriori updates)
    Yminus = H*Xminus + Vmean;

    K = Pminus*H'/(H*Pminus*H'+ R');                        % Kalman gain
    Pplus = (eye(L)-K*H)*Pminus*(eye(L)-K*H)'+K*R*K';       % Stabilized Kalman cov. matrix
    innovations(k) = x(k) - Yminus;
    Xplus = Xminus + K*(innovations(k));                    % A posteriori state estimate
    
    mem2 = [innovations(k).^2 ; mem2(1:end-1)];             % Observation covariance matrix update (for nonstationary signals)
    R = gamma*R + (1-gamma)*mean(mem2);

    Pplus = (Pplus + Pplus')/2;                             % Prevent nonpositive definiteness of the covariance matrix                           
    
    % Time update (A priori updates)
    Xminus = A*Xplus + Wmean;                               % State update
    Pminus = A*Pplus*A' + B*Q*B';                           % Cov. matrix update

    % Store results
    Xhat(:,k) = Xplus';
    Phat(:,:,k) = Pplus';
    Kgain(:,k) = K;
end

%//////////////////////////////////////////////////////////////////////////
% Backward Smoothing Stage
PSmoothed = zeros(size(Phat));
X = zeros(size(Xhat));
PSmoothed(:,:,Samples) = Phat(:,:,Samples);
X(:,Samples) = Xhat(:,Samples);
for k = Samples-1 : -1 : 1,
    S = Phat(:,:,k) * A' * inv(Pbar(:,:,k+1));
    X(:,k) = Xhat(:,k) + S * (X(:,k+1) - Xbar(:,k+1));
    PSmoothed(:,:,k) = Phat(:,:,k) - S * (Pbar(:,:,k+1) - PSmoothed(:,:,k+1)) * S';
    
    PSmoothed(:,:,k) = (PSmoothed(:,:,k) + PSmoothed(:,:,k)')/2;
end

y1 = x - Xhat(1,:);
y2 = x - X(1,:);

% Eigs = eig(A - A*Kgain(:,end)*H)
