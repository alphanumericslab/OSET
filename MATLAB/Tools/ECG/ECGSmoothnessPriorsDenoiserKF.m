function [y1,y2,Pbar,Phat,PSmoothed,Kgain,alpha] = ECGSmoothnessPriorsDenoiserKF(x, order, R, Q, gamma, Qadapt, wlen)
% % % % % % %
% % % % % % % [y1,y2,Pbar,Phat,PSmoothed,Kgain] = ECGSmoothnessPriorsDenoiserKF(x,b,a,q,r,gamma,wlen,mode),
% % % % % % % ECG Kalman filter using smoothing priors
% % % % % % %
% % % % % % % inputs:
% % % % % % % x: vector of noisy signals contaminated with lowpass noise
% % % % % % % b: numerator of the model transfer function
% % % % % % % a: denominator of the model transfer function
% % % % % % % q: covariance of the model errors. By default: var(x)
% % % % % % % r: covariance of the power-line noise. By default: var(x)
% % % % % % % gamma: nonstationarity addaptation parameter. 0 <= gamma <= 1 and by default gamma = 1
% % % % % % % wlen: nonstationarity addaptation window length. by default wlen = length(x)/10
% % % % % % %
% % % % % % % output:
% % % % % % % y1: vector of denoised signal using the Kalman filter
% % % % % % % y2: vector of denoised signal using the Kalman smoother
% % % % % % % Pbar: covariance matrix of the a priori error vector of the Kalman filter
% % % % % % % Phat: covariance matrix of the a posteriori error vector of the Kalman filter
% % % % % % % PSmoothed: covariance matrix of the a posteriori error vector of the Kalman smoother
% % % % % % % Kgain: vector of Kalman filter gain
% % % % % % %
% % % % % % % Open Source ECG Toolbox, version 2.0, June 2015
% % % % % % % Released under the GNU General Public License
% % % % % % % Copyright (C) 2008  Reza Sameni
% % % % % % % Sharif University of Technology, Tehran, Iran -- GIPSA-LAB, INPG, Grenoble, France
% % % % % % % reza.sameni@gmail.com
% % % % % %
% % % % % % % This program is free software; you can redistribute it and/or modify it
% % % % % % % under the terms of the GNU General Public License as published by the
% % % % % % % Free Software Foundation; either version 2 of the License, or (at your
% % % % % % % option) any later version.

%//////////////////////////////////////////////////////////////////////////
% input arguments

% % % % % % % if(nargin>3  && ~isempty(varargin{1})),
% % % % % % %     Q = varargin{1};
% % % % % % % else
% % % % % % %     Q = var(x);
% % % % % % % end
% % % % % % %
% % % % % % % if(nargin>4 && ~isempty(varargin{2})),
% % % % % % %     R = varargin{2};
% % % % % % % else
% % % % % % %     R = var(x);
% % % % % % % end
% % % % % % %
% % % % % % % if(nargin>5  && ~isempty(varargin{3})),
% % % % % % %     gamma = varargin{3};
% % % % % % % else
% % % % % % %     gamma = 1;
% % % % % % % end
% % % % % % %
% % % % % % % if(nargin>6  && ~isempty(varargin{4})),
% % % % % % %     RVarWinlen = varargin{4};
% % % % % % % else
% % % % % % %     RVarWinlen = round(length(x)/10);
% % % % % % % end
% % % % % % %
% % % % % % % if(nargin>7  && ~isempty(varargin{5})),
% % % % % % %     mode = varargin{5};
% % % % % % % else
% % % % % % %     mode = 'obsvr';
% % % % % % % end
% % % % % % %
% % % % % % % if(length(b)>=length(a))
% % % % % % %     error('model is not strictly proper');
% % % % % % % end
% % % % % % % %//////////////////////////////////////////////////////////////////////////
% % % % % % % % Kalman filter parameters
% % % % % % %
% % % % % % % % normalize the coefficients
% % % % % % % b = b/a(1);
% % % % % % % a = a/a(1);
% % % % % % %
% % % % % % % b = b(:);
% % % % % % % b = [zeros(length(a)-2,1); b]; % make the model stricktly proper
% % % % % % % a = a(:);
% % % % % % %
% % % % % % % if(strcmp(mode,'contr'))   % controller model
% % % % % % %     order = length(a)-1;
% % % % % % %     A = [-a(2:end)'; eye(order-1) zeros(order-1,1)];
% % % % % % %     B = [1 ; zeros(order-1,1)];
% % % % % % %     H = b';
% % % % % % %
% % % % % % %     %     [bb aa] = ss2tf(A, B, H, 1)
% % % % % % % else                % observer model
% % % % % % %     order = length(a)-1;
% % % % % % %     A = [-a(2:end) [eye(order-1); zeros(1,order-1)]];
% % % % % % %     B = b;
% % % % % % %     H = [1 zeros(1,order-1)];
% % % % % % %
% % % % % % %     %     [bb aa] = ss2tf(A, B, H, length(b))
% % % % % % % end

% order = 2;
% fs = 1000;
% Q = 1;%2e-8*sqrt(fs^order);%1e-2;
% R = 1000;
RVarWinlen = wlen;%round(length(x)/10);
QVarWinlen = wlen;%round(length(x)/10);
% gamma = 0.9;
mode = 'obsvr';
% Qadapt = 1;

% Kalman filter parameters
if(strcmp(mode,'contr')) % controller model
    b = 1;
    a = diff([zeros(1, order) 1 zeros(1, order)], order);
    
    b = b(:);
    b = [zeros(length(a)-2,1); b]; % make the model stricktly proper
    a = a(:);
    order = length(a)-1;
    A = [-a(2:end)'; eye(order-1) zeros(order-1,1)];
    B = [1 ; zeros(order-1,1)];
    H = b';
    %     [bb aa] = ss2tf(A, B, H, 1)
elseif(strcmp(mode,'obsvr')) % observer model
    b = 1;
    a = diff([zeros(1, order) 1 zeros(1, order)], order);
    
    b = b(:);
    b = [zeros(length(a)-2,1); b]; % make the model stricktly proper
    a = a(:);
    order = length(a)-1;
    A = [-a(2:end) [eye(order-1); zeros(1,order-1)]];
    B = b;
    H = [1 zeros(1,order-1)];
    %     [bb aa] = ss2tf(A, B, H, length(b))
elseif(strcmp(mode,'adhoc')) % ad hoc model
    coefs = diff([zeros(1, order-1) 1 zeros(1, order)], order);
    coefs = -coefs(end:-1:1);
    A = [zeros(order - 1, 1)  eye(order-1) ; zeros(1, order - length(coefs)) coefs];
    B = [zeros(order - 1,1) ; 1];
    % H = [zeros(1, order - 1) 1];
    H = [1 zeros(1, order - 1)];
else
    error('Unknown mode parameter');
end

% A
% B
% H

X0 = [x(1); zeros(order-1,1)];
P0 = 100*B*Q*B';

Xminus = X0;
Pminus = P0;
N = length(x);

alpha = ones(1, N);
mem1 = ones(1, QVarWinlen);
mem2 = zeros(RVarWinlen,1) + R;
Xhat = zeros(order, N);
innovations = zeros(1,N);
Phat = zeros(order,order,N);
Xbar = zeros(order,N);
Pbar = zeros(order,order,N);
Kgain = zeros(order,N);
Q0 = Q;
%//////////////////////////////////////////////////////////////////////////
% Forward Filtering Stage
for k = 1 : N
    % Store results
    Xbar(:,k) = Xminus;
    Pbar(:,:,k) = Pminus;
    
    % Measurement update (A posteriori updates)
    Yminus = H*Xminus;
    M = H*Pminus*H'+R;
    
    K = Pminus*H'/M;                        % Kalman gain
    Pplus = (eye(order)-K*H)*Pminus*(eye(order)-K*H)'+K*R*K';       % Stabilized Kalman cov. matrix
    
    innovations(k) = x(k) - Yminus;
    Xplus = Xminus + K*innovations(k);                    % A posteriori state estimate
    
    if(gamma < 1)
        mem2 = [innovations(k).^2 ; mem2(1:end-1)];             % Observation covariance matrix update (for nonstationary signals)
        R = gamma*R + (1 - gamma)*mean(mem2);
    end
    if(Qadapt)
        indx = max(k - QVarWinlen, 1) : k;
        Q = Q0*(sqrt(length(x)/length(x(indx)))*norm(x(indx))/norm(x))^2;
    end
    
    % Monitoring the innovation variance
    mem1 = [innovations(k)^2/M , mem1(1:end-1)];
    alpha(k) = mean(mem1);
    
    % Time update (A priori updates)
    Xminus = A*Xplus;                                       % State update
    Pminus = A*Pplus*A' + B*Q*B';                           % Cov. matrix update
    
    % Store results
    Xhat(:,k) = Xplus;
    Phat(:,:,k) = Pplus;
    Kgain(:,k) = K;
end

%//////////////////////////////////////////////////////////////////////////
% Backward Smoothing Stage
PSmoothed = zeros(size(Phat));
X = zeros(size(Xhat));
PSmoothed(:,:,N) = Phat(:,:,N);
X(:,N) = Xhat(:,N);
for k = N-1 : -1 : 1
    S = Phat(:,:,k) * A' / Pbar(:,:,k+1);
    X(:,k) = Xhat(:,k) + S * (X(:,k+1) - Xbar(:,k+1));
    PSmoothed(:,:,k) = Phat(:,:,k) - S * (Pbar(:,:,k+1) - PSmoothed(:,:,k+1)) * S';
end

y1 = H * Xhat;      % Filtering results
y2 = H * X;         % Smoothing results

% %//////////////////////////////////////////////////////////////////////////
% % H-infinity filter
% % iinnovations = zeros(1,N);
% % XXbar = zeros(order,N);
% % PPbar = zeros(order,order,N);
% XXhat = zeros(order, N);
% PPhat = zeros(order,order,N);
% KKgain = zeros(order,N);
% S = 1;
% % % % H = [zeros(1, order - 1) 1]; % ????
% LL = H;
% theta = 1;
% XX = [x(1) ; zeros(order-1,1)];
% P = 100*B*Q*B';
% margin = zeros(order,N);
% for k = 1 : N,
%     if(gamma < 1)
%         mem2 = [innovations(k).^2 ; mem2(1:end-1)];             % Observation covariance matrix update (for nonstationary signals)
%         R = gamma*R + (1 - gamma)*mean(mem2);
%     end
%     if(Qadapt)
%         indx = max(k - QVarWinlen, 1) : k;
%         Q = Q0*(sqrt(length(x)/length(x(indx)))*norm(x(indx))/norm(x))^2;
%     end
%     
%     Sbar = LL'*S*LL;
%     M = eye(order) - theta*Sbar*P + H'/(R)*H*P;
%     K = (P /(M))*H'/(R);
%     XX = A*XX + A * K * (x(k) - H*XX);
%     P = A*(P/(M))*A' + B*Q*B';
%     
%     [~, DD] = eig(inv(P) - theta*Sbar + H'/(R)*H);
%     margin(:, k) = sort(diag(DD));
%     
%     % Store results
%     XXhat(:,k) = XX;
%     PPhat(:,:,k) = P;
%     KKgain(:, k) = K;
% end
% y3 = H * XXhat; % H-infinity filtering results