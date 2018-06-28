% Test program for fetal ECG cancellation using EKF and EKS
%
% Dependencies: The baseline wander and ECG filtering toolboxes of the Open Source ECG Toolbox
%
% Sample code for testing naive maternal template subtraction technique
%
% Reza Sameni (C)
% Email: rsameni@shirazu.ac.ir
% Web: www.sameni.info
%
% Crated 2007
% Modified June 2018


%/////////////////////////////////////////////////////////////////////////
clc
clear;
close all;

%/////////////////////////////////////////////////////////////////////////
% initialization
fs = 1600;          % sampling rate
N = 16000;           % segment length
n0 = 1;             % start of segment
t = (0:N-1)/fs;     % time vector
L = 6;              % no. of channels

%/////////////////////////////////////////////////////////////////////////
% abdominal leads
load('CC20060830_ch_ab1');
data = ch_ab1(n0:n0+N-1);      clear ch_ab1;

% reference maternal ECG channel
load('CC20060830_ch_mat');
ref = ch_mat(n0:n0+N-1);      clear ch_mat;

%/////////////////////////////////////////////////////////////////////////
% baseline wander removal of the reference channel
b = LPFilter(ref,.5/fs);
ref = ref - b;

%/////////////////////////////////////////////////////////////////////////
% R-peak detection (may be replaced by better algorithms but requires the
% precise detection of the maternal R-peak)

f = 1.5;    % approximate maternal heart rate (the algorithm is not very sensitive to its precise value)
flag = 1;   % detect positive peaks
peaks = PeakDetection(ref,f/fs,flag);

%/////////////////////////////////////////////////////////////////////////
% maternal ECG phase calculation

[phase, phasepos] = PhaseCalculation(peaks);     % phase calculation
teta = 0;                                       % phase shift
pphase = PhaseShifting(phase,teta);             % phase shifting
bins = fs/4;                                    % number of phase bins

%/////////////////////////////////////////////////////////////////////////
 % baseline wander removal of the data channel (may be replaced by other approaches)

% bsline = LPFilter(data,.5/fs);
bsline = BaseLine2(data,fs*.2,fs*.6,'md');
x = data - bsline;

%/////////////////////////////////////////////////////////////////////////
% calculation of the ensemble average of the maternal ECG contaminated over
% the abdominal lead

[ECGmean,ECGsd,meanphase] = MeanECGExtraction(x,pphase,bins,1);     % mean ECG extraction

%/////////////////////////////////////////////////////////////////////////
% GUI interface for training the parameters of the maternal ECG ensemble average

ECGBeatFitter(ECGmean,ECGsd,meanphase,'OptimumParams');               % ECG beat fitter GUI
% [OptimumParams,mdl,error,approach] = ECGBeatFitterAuto(ECGmean,meanphase);    % Automatic ECG model training (not always reliable)

%//////////////////////////////////////////////////////////////////////////
% Kalman filter parameters

N = length(OptimumParams)/3;            % number of Gaussian kernels
JJ = find(peaks);
fm = fs./diff(JJ);                      % heart-rate
w = mean(2*pi*fm);                      % average heart-rate in rads.
wsd = std(2*pi*fm,1);                   % heart-rate standard deviation in rads.

y = [phase ; x];                         % observation vector

X0 = [-pi 0]';                          % initial state vector
P0 = [(2*pi)^2 0 ;0 (10*max(abs(x))).^2];   % initial state vector covariance matrix
Q = diag( [ (.1*OptimumParams(1:N)).^2 (.05*ones(1,N)).^2 (.05*ones(1,N)).^2 (wsd)^2 , (.01*mean(ECGsd(1:round(end/10))))^2] ); % process noise covariance matrix
R = [(w/fs).^2/12 0 ;0 (mean(ECGsd(1:round(end/10)))).^2];  % observation noise covariance matrix
Wmean = [OptimumParams w 0]';           % process noise vector mean
Vmean = [0 0]';                         % observation noise vector mean

Inits = [OptimumParams w fs];           % initial parameters for the function
InovWlen = ceil(.5*fs);                 % innovations monitoring window length
tau = [];                               % Kalman filter forgetting time. tau=[] for no forgetting factor
gamma = .7;                             % observation covariance adaptation-rate. 0<gamma<1 and gamma=1 for no adaptation
RadaptWlen = ceil(fs/2);                % window length for observation covariance adaptation

[Xekf,Phat,Xeks,PSmoothed,ak] = EKSmoother(y,X0,P0,Q,R,Wmean,Vmean,Inits,InovWlen,tau,gamma,RadaptWlen,1);  % The EKF and EKS

ECG_ekf = Xekf(2,:);                       % the ECG extracted by the EKF method
ECG_eks = Xeks(2,:);                       % the ECG extracted by the EKS method

Denoised_ekf = data - Xekf(2,:)-bsline;   % the denoised signal (Maternal ECG removed), by EKF
Denoised_eks = data - Xeks(2,:)-bsline;   % the denoised signal (Maternal ECG removed), by EKS

%//////////////////////////////////////////////////////////////////////////
% Display Results

figure;
plot(t,ref);
hold on;
plot(t,-peaks*100,'r');
grid;
legend('Reference Channel used for R-wave detection','detected R-waves');
xlabel('time(s)');
ylabel('Amplitude(mV)');

figure;
subplot(211);
plot(t, data(1,:),'b');
hold on;
plot(t, ECG_ekf,'r');
plot(t, ECG_eks,'m');
grid;
legend('Original Signal','Maternal ECG estimated by EKF','Maternal ECG estimated by EKS');
xlabel('time(s)');
ylabel('Amplitude(mV)');

subplot(212);
plot(t,data(1,:),'b');
hold on;
plot(t, Denoised_ekf,'r');
plot(t, Denoised_eks,'m');
plot(t, bsline,'y','linewidth',2);
grid;
legend('Original Signal','Fetal ECG by EKF','Fetal ECG by EKS','Baseline wander estimate');
xlabel('time(s)');
ylabel('Amplitude(mV)');
