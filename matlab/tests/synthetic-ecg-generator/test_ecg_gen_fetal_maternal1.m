%
% Test program for generating synthetic maternal and fetal ECG mixtures
% with realistic ECG noises and extracting independent components using the
% JADE ICA algorithm.
% The results may be compared in different SNRs for the maternal ECG and
% the fetus. When applying ICA in high input SNRs, ICA tends to extract
% three maternal and three fetal components, which correspond with the
% three degrees of freedom in the dipole model. However, as the input noise
% is increased, the weak fetal components vanish in noise and there are
% traces of the fetal components in the maternal ones. This implies that
% ICA fails to separate the maternal and fetal sub-spaces in low input
% SNRs.
%
% Dependencies: The synthetic ECG generator and noise generator package of
%   the Open Source ECG Toolbox. For ICA decomposition, the JADE ICA is
%   also required.
%
% Open Source ECG Toolbox, version 1.0, November 2006
% Released under the GNU General Public License
% Copyright (C) 2006  Reza Sameni
% Sharif University of Technology, Tehran, Iran -- LIS-INPG, Grenoble, France
% reza.sameni@gmail.com

%//////////////////////////////////////////////////////////////////////////
clc
close all;
clear;
% randn('state',2);

%//////////////////////////////////////////////////////////////////////////
decompose = 1; % 1 for decomposition using ICA, 0 for simple results

%//////////////////////////////////////////////////////////////////////////
% General parameters
N = 10000;       % # of signal samples
NN = 10000;      % # of signal samples to plot <= N
fs = 1000;       % desired sampling rate

% Noise parameters
snr = 20;
w_bw = 5;       % weight of baseline wander noise in the generated noise
w_em = 1;       % weight of electrode movement noise in the generated noise
w_ma = 1;       % weight of muscle artifact noise in the generated noise

% Heart locations
pos_m = [-25 7 20];     % Maternal heart location
pos_f = [-15 -4 2];     % Fetal heart location

% Electrode pair locations
ElecPos = [-5 -7 7 ; -5 -7 -7 ; -5 7 7 ; -5 7 -7 ; -5 -1 -5 ; -10 10 18; -10 0 15; -10 10 15];
ElecNeg = [0 0 0 ; 0 0 0 ; 0 0 0 ; 0 0 0 ; 0 0 0 ; -35 10 18; -10 10 15; -10 10 24];
NumCh = size(ElecPos,1);

% Maternal dipole parameters
F_m = .9;                     % maternal heart rate
k_m = 1;                      % maternal dipole attenuation parameter
R_m = rotation_matrix_3d(0,0,0);        % maternal dipole rotation matrices (tetax,tetay,tetaz)
Lambda_m = eye(3);

teta0_m = pi/3;

tetai_m.x  = [-1.09  -0.83   -0.19     -.07  0 .06        0.22    1.2 1.42 1.68];
alphai_m.x = [0.03   .08    -0.13    .85 1.11 .75     0.06   0.1  0.17 0.39];
bi_m.x     = [0.0906    0.1057    0.0453    0.0378    0.0332    0.0302    0.0378    0.6040 0.3020    0.1812];

tetai_m.y  = [-1.1  -0.9 -0.76       -0.11   -.01       0.065  0.8      1.58];
alphai_m.y = [0.035 0.015 -0.019     0.32    .51     -0.32    0.04   0.08];
bi_m.y     = [0.07  .07  0.04        0.055    0.037    0.0604  0.450  0.3];

tetai_m.z  = [-1.1  -0.93 -0.7      -.4     -0.15    .095    1.05 1.25 1.55];
alphai_m.z = [-0.03 -0.14 -0.035    .045     -0.4    .46    -.12 -.2 -.35];
bi_m.z     = [.03  .12  .04         .4    .045       .05    .8 .4 .2];

% Fetal dipole parameters
F_f = 2.2;                          % fetal heart rate
k_f = .1;                            % fetal dipole attenuation parameter
R_f = rotation_matrix_3d(-3*pi/4,0,-pi/2);     % fetal rotation matrices (tetax,tetay,tetaz)
Lambda_f = eye(3);

teta0_f = -pi/2;

tetai_f.x  = [-0.7    -0.17    0       0.18     1.4];
alphai_f.x = .1*[0.07     -0.11   1.3     0.07   0.275];
bi_f.x     = [.1       .03     .045     0.02    0.3];

tetai_f.y  = [-0.9     -0.08   0       0.05        1.3];
alphai_f.y = .1*[0.04     0.3     .45     -0.35       0.05];
bi_f.y     = [.1       .05      .03    .04         .3];

tetai_f.z  = [-0.8      -.3     -0.1        .06     1.35];
alphai_f.z = .1*[-0.14    .03     -0.4        .46     -0.1];
bi_f.z     = [.1       .4      .03         .03     .3];

%//////////////////////////////////////////////////////////////////////////
% Noise generation
M = 10800;
% original noise template
template =  biosignal_noise_gen(5,1,0,M,360,[w_bw,w_em,w_ma],1000);

% parameters required for estimating the AR coefficients using a Kalman Filter (KF)
order = 12;                         % AR model order for modeling the ECG noise
[a0,e] = aryule(template,order);    % a global AR model
q = (.05*max(abs(a0(2:end))))^2;    % AR coefficients covariance
R = 1;                              % unit variance noise
p0 = 1e6*q;                         % covariance of the KF initial state
alpha = 1;                          % KF forgetting factor

% time-variant AR parameter estimation using KF and Kalman Smoother (KS)
[Ahat,Asmoothed] = time_variant_ar_kalman(template,order,a0(2:end)',q,R,p0,alpha);

% generating different instances of ECG noise using the time-variant AR parameters
noise =  zeros(NumCh,N);
for j = 1:NumCh
    x = randn(1,M);
    y1 =  zeros(M,1);
    for i = order+1:M-1
        y1(i) = (sqrt(1)*x(i)-Ahat(:,i)'*y1(i-1:-1:i-order))/1;         % KF
    end
    % resampling the noise matrix to the desired sampling rate
    n1 = resample((y1-mean(y1))/std(y1),fs,360);
    noise(j,:) = n1(101:N+100);
end

%//////////////////////////////////////////////////////////////////////////
% ECG calculation
H_m = zeros(NumCh, 3);
H_f = zeros(NumCh, 3);
for i = 1:NumCh
    for j = 1:3
        H_m(i,j) = k_m* ((ElecPos(i,j)-pos_m(j))/sqrt(sum((ElecPos(i,:)-pos_m).^2))^3 - (ElecNeg(i,j)-pos_m(j))/sqrt(sum((ElecNeg(i,:)-pos_m).^2))^3);
        H_f(i,j) = k_f* ((ElecPos(i,j)-pos_f(j))/sqrt(sum((ElecPos(i,:)-pos_f).^2))^3 - (ElecNeg(i,j)-pos_f(j))/sqrt(sum((ElecNeg(i,:)-pos_f).^2))^3);
    end
end
[DIP_m, teta_m] = vcg_gen_direct_sum(N,fs,F_m,alphai_m,bi_m,tetai_m,teta0_m);
[DIP_f, teta_f] = vcg_gen_direct_sum(N,fs,F_f,alphai_f,bi_f,tetai_f,teta0_f);
VCG_m = R_m*Lambda_m*[DIP_m.x ; DIP_m.y ; DIP_m.z];
VCG_f = R_f*Lambda_f*[DIP_f.x ; DIP_f.y ; DIP_f.z];
s0 = H_m*VCG_m + H_f*VCG_f;
s = s0 + (sqrt(sum(s0.^2,2))./sqrt(sum(noise.^2,2))/sqrt(10^(snr/10))*ones(1,size(s0,2))).*noise;

%//////////////////////////////////////////////////////////////////////////
% Mixture decomposition and plotting
time = (0:N-1)/fs;
if(decompose == 0)
    % plotting results
    figure;
    plot(time,1000*s');
    grid
    xlabel('time(s)');
    ylabel('Amplitude(mV)');
    title('Synthetic maternal and fetal ECG mixtures with realistic ECG noises.');
    
else
    %//////////////////////////////////////////////////////////////////////////
    % Independent component analysis using the JADE algorithm
    
    % fastICA
    % [s2, A, W] = fastica(s, 'displayMode', 'off');
    
    % JADE ICA
    W =  jadeR(s);
    s2 = real(W*s);
    A = pinv(W);
    
    % plotting results
    chindex = 1:NumCh;
    
    %//////////////////////////////////////////////////////////////////////////
    for i = 1:size(s,1)
        figure;
        plot(time(1:NN),1000*s(i,1:NN),'k','Linewidth',.5);
        xlabel('time(s)','FontSize',10);
        ylabel(['Ch_',num2str(chindex(i)),'(mV)'],'FontSize',10);
        grid;
        set(gca,'FontSize',10);
    end
    
    %//////////////////////////////////////////////////////////////////////////
    for i = 1:size(s2,1)
        figure;
        plot(time(1:NN),s2(i,1:NN),'k','Linewidth',.5);
        xlabel('time(s)','FontSize',10);
        ylabel(['IC_',num2str(chindex(i))],'FontSize',10);
        grid;
        set(gca,'FontSize',10);
    end
    
    %//////////////////////////////////////////////////////////////////////////
    VV1 = s(3,:)-s(1,:);
    VV2 = s(4,:)-s(2,:);
    VV3 = s(6,:);
    h = figure;
    plot3(1000*VV1,1000*VV2,1000*VV3,'k','Linewidth',.5);
    grid;
    view([30,15])
    
    xlabel('Ch_3-Ch_1 (mV)','FontSize',10);
    ylabel('Ch_4-Ch_2 (mV)','FontSize',10);
    zlabel('Ch_6 (mV)','FontSize',10);
    set(gca,'FontSize',10);
end