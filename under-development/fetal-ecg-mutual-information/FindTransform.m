clear;
close all;

load('FOETAL_ECG.dat');
data = FOETAL_ECG(:,2:end)';
time = FOETAL_ECG(:,1)';
clear FOETAL_ECG;
I = 1251:1300;
%//////////////////////////////////////////////////////////////////////////
N = length(I);
fs = 250;
%time = [0:N-1]/fs;

% mother params
F_m = 1.359;
%teta0_m = -pi/3.05;
teta0_m = -pi/3.4;
teta_iso = -1.4;


%//////////////////////////////////////////////////////////////////////////
% tetai_m.x  = [-0.7    -0.17    0       0.18     1.4];     % [Parameters for SingleECG2.mat]
% alphai_m.x = 100*[0.08     -0.12   1.4     0.08   0.585];
% bi_m.x     = 1.51*[.1       .03     .045     0.02    0.3];
% 
% tetai_m.y  = [-0.9     -0.08   0       0.05        1.3];
% alphai_m.y = 100*[0.04     0.3     .45     -0.35       0.09];
% bi_m.y     = 1.51*[.1       .05      .03    .04         .3];
% 
% tetai_m.z  = [-0.8      -.3     -0.1        .06     1.35];
% alphai_m.z = 100*[-0.14    .03     -0.4        .46     -0.6];
% bi_m.z     = 1.51*[.1       .4      .03         .03     .3];

%//////////////////////////////////////////////////////////////////////////
tetai_m.x  = [-0.55 -0.14 0 0.1 1.1];     % [Parameters for SingleECG1.mat]
alphai_m.x = 100*[0.09 -0.2 1.32 -0.05 0.48];
bi_m.x     = 1.51*[.1 .05 .04 .08 .15];

tetai_m.y  = [-0.7 -0.08 0 0.05 1];
alphai_m.y = 100*[0.03 0.3 .45 -0.35 0.08];
bi_m.y     = 1.51*[.1 .05 .03 .04 .25];

tetai_m.z  = [-0.6 -.3 -0.1 .06 1];
alphai_m.z = 100*[-0.15 .03 -0.4 .46 -0.6];
bi_m.z     = 1.51*[.1 .4 .03 .03 .2];

% % % tetai_m.x  = [-0.9 -0.24 0 0.15 1.5];     % [Parameters for SingleECG1.mat]
% % % alphai_m.x = 100*[0.09 -0.3 1.32 -0.3 0.4];
% % % bi_m.x     = 1.51*[.07 .05 .09 .18 .15];
% % % 
% % % tetai_m.y  = [-0.7 -0.08 0 0.15 1];
% % % alphai_m.y = 100*[0.03 0.3 .6 -0.35 0.08];
% % % bi_m.y     = 1.51*[.1 .05 .03 .04 .25];
% % % 
% % % tetai_m.z  = [-0.9 -.3 -0.1 .06 1]+.1;
% % % alphai_m.z = 100*[-0.15 .06 -0.4 .46 -0.6];
% % % bi_m.z     = 1.51*[.1 .4 .08 .05 .2];
% % % 
%//////////////////////////////////////////////////////////////////////////
[DIPm tetam] = DipoleGenerator(N,fs,F_m,alphai_m,bi_m,tetai_m,teta0_m,teta_iso);
s = [DIPm.x;DIPm.y;DIPm.z];

%H = data(1:5,I)/s;
H = data(1:5,I)/data(6:8,I);
s2 = H*s;
%//////////////////////////////////////////////////////////////////////////
L = 1;
close all;
%for i = 1:size(data,1),
for i = [1:5],
    figure(floor((i-1)/L)+1);
    subplot(L,1,mod(i-1,L)+1);
    plot(s2(i,:),'b');
    hold on;
    plot(data(i,I),'r');
    %plot(s(:,:)','g');
    grid;
end

% % % figure;
% % % plot(time,data(1,:));
% % % grid

%//////////////////////////////////////////////////////////////////////////
% % % [DIPm tetam] = DipoleGenerator(10000,fs,F_m,alphai_m,bi_m,tetai_m,teta0_m,teta_iso);
% % % s = [DIPm.x;DIPm.y;DIPm.z];
% % % figure
% % % plot(s');
% % % grid;
