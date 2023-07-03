function [ARCoefs , ARCoefs2]= TimeVariantAR(data,order,x0,q,R,p0,alpha)
%
% [ARCoefsKF , ARCoefsKS]= TimeVariantAR(data,order,x0,q,R,p0,alpha);
% Time variant auto-regressive(AR) model estimated by Kalman Filter and Kalman Smoother
%
% inputs:
% data: template noise used for model training
% order: AR model order
% x0: a time-invariant set of AR-coefficients estimated by applying a
%     global AR-model estimation on the entire input signal
% q: AR coefficients covariance
% R: noise variance
% p0: covariance of the KF initial state
% alpha: KF forgetting factor (alpha=1 for standard KF)
%
% outputs:
% ARCoefsKF: AR coefficients estimated by a Kalman Filter
% ARCoefsKS: AR coefficients estimated by a Kalman Smoother
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

N = length(data);
A = 1;
Q = q*eye(order);
Wmean = zeros(order,1);
Vmean = mean(data);
ARCoefs = zeros(order,N);
Xminus = x0;
Pminus = p0*eye(order);
% P = zeros(N,1);

Pbar = zeros(order,order,N);
Phat = zeros(order,order,N);
Xhat = zeros(order,N);
Xbar = zeros(order,N);

% Filtering
for i = 1:N
    Pbar(:,:,i) = Pminus;
    Xbar(:,i) = Xminus;

    if(i<order+1)
        H = [-data(i-1:-1:1)' zeros(1,order-i+1)];
    else
        H = -data(i-1:-1:i-order)';
    end

    Yminus = H * Xminus + Vmean;
    inov = data(i)-Yminus;

    K = Pminus * H'/(H * Pminus * H' + alpha*R);
    Pplus = ( (eye(order) - K * H) * Pminus * (eye(order) - K * H)' + K * R * K' )/alpha;
    Xplus = Xminus + K*(inov);                                % A posteriori state estimate

    Xminus = A * Xplus + Wmean;                              % State update
    Pminus = A * Pplus * A' + Q;

    ARCoefs(:,i) = Xplus;
%     P(i) = max(diag(Pplus));
    Phat(:,:,i) = Pplus;
    Xhat(:,i) = Xplus;
end

% Smoothing
PSmoothed = zeros(size(Phat));
X = zeros(size(Xhat));
PSmoothed(:,:,N) = Phat(:,:,N);
X(:,N) = Xhat(:,N);
for k = N-1:-1:1
    S = Phat(:,:,k) * A' /Pbar(:,:,k+1);
    X(:,k) = Xhat(:,k) + S * (X(:,k+1) - Xbar(:,k+1));
    PSmoothed(:,:,k) = Phat(:,:,k) - S * (Pbar(:,:,k+1) - PSmoothed(:,:,k+1)) * S';
end

ARCoefs2 = X;