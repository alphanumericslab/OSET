%
% Test program for generating synthetic maternal and twin fetal ECG
% mixtures with realistic ECG noises, and extracting independent components
% using the JADE ICA algorithm.
% The results may be compared in different SNRs for the maternal ECG and
% the fetuses. In this example, since there are eight recording channels
% and nine degrees of freedom in the cardiac dipoles (three for each
% cardiac dipole), ICA is not able to separate all the dipole components.
% This leads to components which are mixtures of the signals from both
% fetuses.
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

%//////////////////////////////////////////////////////////////////////////
clc
close all;
clear;
% randn('state',2);

%//////////////////////////////////////////////////////////////////////////
decompose = 1; % 1 for decomposition using ICA, 0 for simple results

%//////////////////////////////////////////////////////////////////////////
% General parameters
N = 5000;       % # of signal samples
NN = 5000;      % # of signal samples to plot <= N
fs = 500;       % desired sampling rate

% Noise parameters
snr = 0;
w_bw = 5;       % weight of baseline wander noise in the generated noise
w_em = 5;       % weight of electrode movement noise in the generated noise
w_ma = 5;       % weight of muscle artifact noise in the generated noise

% Heart locations
pos_m = [-25 7 20];     % Maternal heart location
pos_f1 = [-15 -4 2];     % Fetal heart location
pos_f2 = [-13 4 1];     % Fetal heart location

% Electrode pair locations
ElecPos = [-5 -7 7 ; -5 -7 -7 ; -5 7 7 ; -5 7 -7 ; -5 -1 -5 ; -10 10 18; -10 0 15; -10 10 15];
ElecNeg = [0 0 0 ; 0 0 0 ; 0 0 0 ; 0 0 0 ; 0 0 0 ; -35 10 18; -10 10 15; -10 10 24];
NumCh = size(ElecPos,1);

% Maternal dipole parameters
F_m = .9;                     % maternal heart rate
k_m = 1;                      % maternal dipole attenuation parameter
R_m = Rotate3D(0,0,0);        % maternal dipole rotation matrices (tetax,tetay,tetaz)
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

% First fetus dipole parameters
F_f1 = 2.1;                          % fetal heart rate
k_f1 = .1;                            % fetal dipole attenuation parameter
R_f1 = Rotate3D(-3*pi/4,0,-pi/2);     % fetal rotation matrices (tetax,tetay,tetaz)
Lambda_f1 = eye(3);

teta0_f1 = -pi/2;

tetai_f1.x  = [-0.7    -0.17    0       0.18     1.4];
alphai_f1.x = .1*[0.07     -0.11   1.3     0.07   0.275];
bi_f1.x     = [.1       .03     .045     0.02    0.3];

tetai_f1.y  = [-0.9     -0.08   0       0.05        1.3];
alphai_f1.y = .1*[0.04     0.3     .45     -0.35       0.05];
bi_f1.y     = [.1       .05      .03    .04         .3];

tetai_f1.z  = [-0.8      -.3     -0.1        .06     1.35];
alphai_f1.z = .1*[-0.14    .03     -0.4        .46     -0.1];
bi_f1.z     = [.1       .4      .03         .03     .3];

% Second fetus dipole parameters
F_f2 = 1.9;                          % fetal heart rate
k_f2 = .1;                            % fetal dipole attenuation parameter
R_f2 = Rotate3D(3*pi/4,pi/20,-pi/2.5);     % fetal rotation matrices (tetax,tetay,tetaz)
Lambda_f2 = eye(3);

teta0_f2 = 0;

tetai_f2.x  = [-0.71    -0.16    0       0.19     1.41];
alphai_f2.x = .1*[0.05     -0.13   1.32     0.06   0.27];
bi_f2.x     = [.12       .04     .04     0.01    0.2];

tetai_f2.y  = [-0.9     -0.08   0       0.05        1.31];
alphai_f2.y = .1*[0.03     0.32     .44     -0.34       0.02];
bi_f2.y     = [.12       .03      .02    .03         .2];

tetai_f2.z  = [-0.9      -.32     -0.11        .05     1.3];
alphai_f2.z = .1*[-0.13    .03     -0.41       .47     -0.09];
bi_f2.z     = [.1       .35      .025         .03     .29];

%//////////////////////////////////////////////////////////////////////////
% Noise generation
M = 10800;
% original noise template
template =  NoiseGenerator(5,1,0,M,360,[w_bw,w_em,w_ma],1000);

% parameters required for estimating the AR coefficients using a Kalman Filter (KF)
order = 12;                         % AR model order for modeling the ECG noise
[a0,e] = aryule(template,order);    % a global AR model
q = (.05*max(abs(a0(2:end))))^2;    % AR coefficients covariance
R = 1;                              % unit variance noise
p0 = 1e6*q;                         % covariance of the KF initial state
alpha = 1;                          % KF forgetting factor

% time-variant AR parameter estimation using KF and Kalman Smoother (KS)
[Ahat,Asmoothed] = TimeVariantAR(template,order,a0(2:end)',q,R,p0,alpha);

% generating different instances of ECG noise using the time-variant AR parameters
noise =  zeros(NumCh,N);
for j = 1:NumCh,
    x = randn(1,M);
    y1 =  zeros(M,1);
    for i = order+1:M-1,
        y1(i) = (sqrt(1)*x(i)-Ahat(:,i)'*y1(i-1:-1:i-order))/1;         % KF
    end
    % resampling the noise matrix to the desired sampling rate
    n1 = resample((y1-mean(y1))/std(y1),fs,360);
    noise(j,:) = n1(101:N+100);
end

%//////////////////////////////////////////////////////////////////////////
% ECG calculation
for i = 1:NumCh,
    for j = 1:3,
        H_m(i,j) = k_m* ((ElecPos(i,j)-pos_m(j))/sqrt(sum((ElecPos(i,:)-pos_m).^2))^3 - (ElecNeg(i,j)-pos_m(j))/sqrt(sum((ElecNeg(i,:)-pos_m).^2))^3);
        H_f1(i,j) = k_f1* ((ElecPos(i,j)-pos_f1(j))/sqrt(sum((ElecPos(i,:)-pos_f1).^2))^3 - (ElecNeg(i,j)-pos_f1(j))/sqrt(sum((ElecNeg(i,:)-pos_f1).^2))^3);
        H_f2(i,j) = k_f2* ((ElecPos(i,j)-pos_f2(j))/sqrt(sum((ElecPos(i,:)-pos_f2).^2))^3 - (ElecNeg(i,j)-pos_f2(j))/sqrt(sum((ElecNeg(i,:)-pos_f2).^2))^3);
    end
end
[DIP_m teta_m] = DipoleGenerator2(N,fs,F_m,alphai_m,bi_m,tetai_m,teta0_m);
[DIP_f1 teta_f1] = DipoleGenerator2(N,fs,F_f1,alphai_f1,bi_f1,tetai_f1,teta0_f1);
[DIP_f2 teta_f2] = DipoleGenerator2(N,fs,F_f2,alphai_f2,bi_f2,tetai_f2,teta0_f2);

VCG_m = R_m*Lambda_m*[DIP_m.x ; DIP_m.y ; DIP_m.z];
VCG_f1 = R_f1*Lambda_f1*[DIP_f1.x ; DIP_f1.y ; DIP_f1.z];
VCG_f2 = R_f2*Lambda_f2*[DIP_f2.x ; DIP_f2.y ; DIP_f2.z];

s0 = H_m*VCG_m + H_f1*VCG_f1 + H_f2*VCG_f2;
s = s0 + (sqrt(sum(s0.^2,2))./sqrt(sum(noise.^2,2))/sqrt(10^(snr/10))*ones(1,size(s0,2))).*noise;

%//////////////////////////////////////////////////////////////////////////
% Mixture decomposition and plotting
time = [0:N-1]/fs;
if(decompose == 0),
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
    for i = 1:size(s,1),
        h = figure;
        plot(time(1:NN),1000*s(i,1:NN),'k','Linewidth',.5);
        xlabel('time(s)','FontSize',10);
        ylabel(['Ch_',num2str(chindex(i)),'(mV)'],'FontSize',10);
        grid;
        set(gca,'FontSize',10);
    end

    %//////////////////////////////////////////////////////////////////////////
    for i = 1:size(s2,1),
        h = figure;
        plot(time(1:NN),s2(i,1:NN),'k','Linewidth',.5);
        xlabel('time(s)','FontSize',10);
        ylabel(['IC_',num2str(chindex(i))],'FontSize',10);
        grid;
        set(gca,'FontSize',10);
    end

    %//////////////////////////////////////////////////////////////////////////
%     VV1 = s(3,:)-s(1,:);
%     VV2 = s(4,:)-s(2,:);
%     VV3 = s(6,:);
    VV1 = s(6,:);
    VV2 = s(7,:);
    VV3 = s(8,:);

    h = figure;
    plot3(1000*VV1,1000*VV2,1000*VV3,'k','Linewidth',.5);
    grid;
    view([30,15])

%     xlabel('Ch_3-Ch_1 (mV)','FontSize',10);
%     ylabel('Ch_4-Ch_2 (mV)','FontSize',10);
%     zlabel('Ch_6 (mV)','FontSize',10);
    xlabel('Ch_6 (mV)','FontSize',10);
    ylabel('Ch_7 (mV)','FontSize',10);
    zlabel('Ch_8 (mV)','FontSize',10);
    set(gca,'FontSize',10);
end