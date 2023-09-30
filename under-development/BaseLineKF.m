function [baseline,baseline1,baseline2] = BaseLineKF(y,fc,winlen1,winlen2,approach)
%
% UNDER DEVELOPMENT
%
% b = BaseLineKF(x,fc),
% Baseline wander extraction from biomedical recordings, using a first
% order Kalman Smoother
%
% inputs:
% x: vector or matrix of input data (channels x samples)
% fc: -3dB cut-off frequency normalized by the sampling frequency
%
% output:
% b: vector or matrix of filtered data (channels x samples)
%
%
% Open Source ECG Toolbox, version 1.0, November 2006
% Released under the GNU General Public License
% Copyright (C) 2006  Reza Sameni
% Sharif University of Technology, Tehran, Iran -- LIS-INPG, Grenoble, France
% reza.sameni@gmail.com

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details. You should have received a copy of the
% GNU General Public License along with this program; if not, write to the
% Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
% MA  02110-1301, USA.

% if(nargin >= 3)
%     flag = varargin{1};
% else
%     flag = 1;
% end

% if(nargin > 2)
%     winlen1 = varargin{1};
%     winlen2 = varargin{2};
%     approach = varargin{3};
% end

baseline = zeros(size(y));
baseline1 = zeros(size(y));
baseline2 = zeros(size(y));
for i = 1:size(y,1)
    %     ymean = mean(y(i,:));
    %     y(i,:) = y(i,:) - ymean;

    % Parameters
    k = .7; % cut-off value
    alpha = (1-k*cos(2*pi*fc)-sqrt(2*k*(1-cos(2*pi*fc))-k^2*sin(2*pi*fc)^2))/(1-k);
    % bw = BaseLine2(y(i,:),winlen1,winlen2,approach);
    bw = lp_filter_zero_phase(y(i,:), fc);
    bw2 = filtfilt(1-alpha,[1 -alpha],bw);

    Q = var(bw);

    R = var(bw-bw2);

    Wmean = 0;
    %Vmean = mean(y(i,:)-bw);
    Vmean = 0;
    gamma = 1;
    % % %     X0 = y(i,1);
    X0 = bw(1);
    P0 = 1e6*Q;

    L = 1;
    A = alpha;
    H = 1;
    Xminus = X0;
    Pminus = P0;
    Samples = length(y(i,:));
    Xhat = zeros(1,Samples);
    Phat = zeros(1,Samples);
    Xbar = zeros(1,Samples);
    Pbar = zeros(1,Samples);
    B = (1-alpha);
    % % %     obs = y(i,:);
    obs = bw;
    %//////////////////////////////////////////////////////////////////////////
    % Filtering
    for k = 1 : Samples
        % Store results
        Xbar(k) = Xminus';
        Pbar(k) = Pminus';

        % Measurement update (A posteriori updates)
        Yminus = H*Xminus + Vmean;

        % Observation available at this sampling time
        K = Pminus*H'/(H*Pminus*H'+ R');                                    % Kalman gain
        Pplus   = (eye(L)-K*H)*Pminus*(eye(L)-K*H)'+K*R*K';                 % Stabilized Kalman cov. matrix
        innovations(k) = obs(k) - Yminus;
        Xplus = Xminus + K*(innovations(k));                                % A posteriori state estimate

        % Time update (A priori updates)
        Xminus = A*Xplus + Wmean;                                           % State update
        Pminus = A*Pplus*A' + B*Q*B';                                       % Cov. matrix update

        % Store results
        Xhat(k) = Xplus';
        Phat(k) = Pplus';
    end

    %//////////////////////////////////////////////////////////////////////////
    % Smoothing
    PSmoothed = zeros(size(Phat));
    X = zeros(size(Xhat));
    PSmoothed(Samples) = Phat(Samples);
    X(Samples) = Xhat(Samples);
    for k = Samples-1 : -1 : 1
        S = Phat(k) * A' * pinv(Pbar(k+1));
        X(k) = Xhat(k) + S * (X(k+1) - Xbar(k+1));
        PSmoothed(k) = Phat(k) - S * (Pbar(k+1) - PSmoothed(k+1)) * S';
    end

    %     b(i,:) = X + ymean;
    baseline(i,:) = X;
    baseline1(i,:) = bw;
    baseline2(i,:) = bw2;
end
