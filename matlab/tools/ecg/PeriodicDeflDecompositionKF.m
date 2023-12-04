function [output, criterion] = PeriodicDeflDecompositionKF(dat,Itr,peaks, bins, f0, fs)
%
% [output criterion] = PeriodicDeflDecompositionKF(dat,Itr,peaks, bins, f0, fs)
% (maternal) cardiac signal supression by deflation using periodic component
% analysis and extended Kalman filtering for nonlinear denoising.
%
% Ref: An implementation of the following paper:
% Sameni, Reza, Christian Jutten, and Mohammad B. Shamsollahi.
% "A deflation procedure for subspace decomposition."
% IEEE Transactions on Signal Processing 58, no. 4 (2010): 2363-2374.
%
% Open Source ECG Toolbox, version 3.14, November 2019
% Released under the GNU General Public License
% Copyright (C) 2019  Reza Sameni
% Shiraz University, Shiraz, Iran
% reza.sameni@gmail.com
% ECG phase calculation

phase = PhaseCalculation(peaks);

% PM time calculation
[T0,T1] = SynchPhaseTimes2(peaks);

% BW removal
bl = LPFilter(dat, f0/fs);
dat = dat - bl;

criterion = zeros(1, Itr);
for i = 1:Itr
    % periodic component analysis stage
    A = dat(:,T0)*dat(:,T1)';
    B = dat(:,T0)*dat(:,T0)';
    
    A = (A+A')/2;
    B = (B+B')/2;
    
    criterion(i) = trace(A)/trace(B);

    [V,D] = eig(A, B, 'chol');
    
    d = diag(D);
    [~, I] = sort(d);
    I = I(end:-1:1);
    
    W = V(:,I)';
    dat = W*dat;
    
    % remove the ECG from the most peridic component
    x = dat(1,:);
    [ECGmean,ECGsd,meanphase] = MeanECGExtraction(x, phase, bins,1); % mean ECG extraction
    
    OptimumParams = ECGBeatFitter(ECGmean,ECGsd,meanphase);           % ECG beat fitter GUI
    
    %//////////////////////////////////////////////////////////////////////////
    N = length(OptimumParams)/3;% number of Gaussian kernels
    JJ = find(peaks);
    fm = fs./diff(JJ);          % heart-rate
    w = mean(2*pi*fm);          % average heart-rate in rads.
    wsd = std(2*pi*fm,1);       % heart-rate standard deviation in rads.
    
    y = [phase ; x];
    
    X0 = [-pi 0]';
    P0 = [(2*pi)^2 0 ;0 (10*max(abs(x))).^2];
    Q = diag( [ (.1*OptimumParams(1:N)).^2 (.05*ones(1,N)).^2 (.05*ones(1,N)).^2 (wsd)^2 , (.03*mean(ECGsd(1:round(length(ECGsd)/10))))^2] );
    R = [(w/fs).^2/12 0 ;0 (mean(ECGsd(1:round(length(ECGsd)/10)))).^2];
    Wmean = [OptimumParams w 0]';
    Vmean = [0 0]';
    Inits = [OptimumParams w fs];
    
    InovWlen = ceil(.5*fs);     % innovations monitoring window length
    tau = [];                   % Kalman filter forgetting time. tau=[] for no forgetting factor
    gamma = 1;                  % observation covariance adaptation-rate. 0<gamma<1 and gamma=1 for no adaptation
    RadaptWlen = ceil(fs/2);    % window length for observation covariance adaptation
    
    %//////////////////////////////////////////////////////////////////////////
    [Xekf, Phat, Xeks, PSmoothed, ak] = EKSmoother(y,X0,P0,Q,R,Wmean,Vmean,Inits,InovWlen,tau,gamma,RadaptWlen,0);
    
    % returning back the residue
    dat(1,:) = x - Xeks(2,:);
    
    dat = pinv(W)*dat;
end

output = dat + bl;