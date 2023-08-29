% clc
clear;
close all;

% % % load('FOETAL_ECG.dat');
% % % data = FOETAL_ECG(:,2:end)';
% % % H1 = data(1:5,:)/data(6:8,:)
% % % clear FOETAL_ECG;
% % % clear data;
% % %
% % % load patient165_s0323lre;
% % % data = data(:,2:16)';
% % % H1 = data(1:12,:)/data(13:15,:)

%//////////////////////////////////////////////////////////////////////////
km = 1;
kf = .1;

%elec = [-2 5 5; -2 -5 5; -2 -5 -5; -2 5 -5 ; -6 5 5; -6 -5 5; -6 -5 -5; -6 5 -5];% ; -25 7 20];
elec = [-6 -10 2;-6 -10 -2;-3 -5 2;-3 -5 -2;-3 5 2;-3 5 -2;-6 10 2;-6 10 -2];
chindex = [1 2 3 4 5 6 7 8];

pos_f = [-13 -4 1];
pos_m = [-15 7 20];

r_f = sqrt(sum(pos_f.^2));
r_m = sqrt(sum(pos_m.^2));
for i = 1:size(elec,1),
    for j = 1:3,
        H_m(i,j) = km* ((elec(i,j)-pos_m(j))/sqrt(sum((elec(i,:)-pos_m).^2))^3 - (0-pos_m(j))/r_m^3);
        H_f(i,j) = kf* ((elec(i,j)-pos_f(j))/sqrt(sum((elec(i,:)-pos_f).^2))^3 - (0-pos_f(j))/r_f^3);
    end
end

FrankPos = [-5 -8 -2; -3 0 -1; -3 -5 -5];
FrankNeg = [-10 -8 -2; -3 -15 -1; -3 -5 5];
for i = 1:3,
    for j = 1:3,
        Frank_m(i,j) = km* ((FrankPos(i,j)-pos_m(j))/sqrt(sum((FrankPos(i,:)-pos_m).^2))^3 - (FrankNeg(i,j)-pos_m(j))/sqrt(sum((FrankNeg(i,:)-pos_m).^2))^3);
        Frank_f(i,j) = kf* ((FrankPos(i,j)-pos_f(j))/sqrt(sum((FrankPos(i,:)-pos_f).^2))^3 - (FrankNeg(i,j)-pos_f(j))/sqrt(sum((FrankNeg(i,:)-pos_f).^2))^3);
    end
end

%//////////////////////////////////////////////////////////////////////////
N = 2500;
fs = 500;

% mother params
% F_m = 1.359;
teta0_m = -pi/3.05;
teta_iso = -1.4;

tetai_m.x  = [-1.09  -0.83   -0.19     -.07  0 .06        0.22    1.2 1.42 1.68];% 2.9];
alphai_m.x = [0.03   .08    -0.13    .85 1.11 .75     0.06   0.1  0.17 0.39];% .03];
bi_m.x     = [0.0906    0.1057    0.0453    0.0378    0.0332    0.0302    0.0378    0.6040 0.3020    0.1812];% .5];

tetai_m.y  = [-1.1  -0.9 -0.76       -0.11   -.01       0.065  0.8      1.58];% 2.9];
alphai_m.y = [0.035 0.015 -0.019     0.32    .51     -0.32    0.04   0.08];% .014];
bi_m.y     = [0.07  .07  0.04        0.055    0.037    0.0604  0.450  0.3];% .5];

tetai_m.z  = [-1.1  -0.93 -0.7      -.4     -0.15    .095    1.05 1.25 1.55];% 2.8];
alphai_m.z = [-0.03 -0.14 -0.035    .045     -0.4    .46    -.12 -.2 -.35];% -.035];
bi_m.z     = [.03  .12  .04         .4    .045       .05    .8 .4 .2];% .4];

% mother params
F_m = 1;
teta0_m = pi;

% % tetai_m.x  = [-0.7    -0.17    0       0.18     1.4];     % [Parameters for SingleECG2.mat]
% % alphai_m.x = [0.08     -0.12   1.4     0.08   0.585];
% % bi_m.x     = [.1       .03     .045     0.02    0.3];
% %
% % tetai_m.y  = [-0.9     -0.08   0       0.05        1.3];
% % alphai_m.y = [0.04     0.3     .45     -0.35       0.09];
% % bi_m.y     = [.1       .05      .03    .04         .3];
% %
% % tetai_m.z  = [-0.8      -.3     -0.1        .06     1.35];
% % alphai_m.z = [-0.14    .03     -0.4        .46     -0.6];
% % bi_m.z     = [.1       .4      .03         .03     .3];

% fetal params
F_f = 1.9;
teta0_f = pi/3;

tetai_f.x  = [-0.7    -0.17    0       0.18     1.4];
alphai_f.x = .05*[0.07     -0.11   1.3     0.07   0.275];
bi_f.x     = [.1       .03     .045     0.02    0.3];

tetai_f.y  = [-0.9     -0.08   0       0.05        1.3];
alphai_f.y = .03*[0.04     0.3     .45     -0.35       0.09];
bi_f.y     = [.1       .05      .03    .04         .3];

tetai_f.z  = [-0.8      -.3     -0.1        .06     1.35];
alphai_f.z = .04*[-0.14    .03     -0.4        .46     -0.2];
bi_f.z     = [.1       .4      .03         .03     .3];


%//////////////////////////////////////////////////////////////////////////
[DIPm tetam] = DipoleGenerator(N,fs,F_m,alphai_m,bi_m,tetai_m,teta0_m);%,teta_iso);
[DIPf tetaf] = DipoleGenerator(N,fs,F_f,alphai_f,bi_f,tetai_f,teta0_f);%,teta_iso);
% s = H_f*[DIPm.x;DIPm.y;DIPm.z];


%Rf = Rotate3D(pi,pi/7,pi/2);   % Fetal rotation matrices (tetax,tetay,tetaz)
%Rf = Rotate3D(pi/8,pi,pi/2);   % Fetal rotation matrices (tetax,tetay,tetaz)
Rf = Rotate3D(pi,0,-pi/2);   % Fetal rotation matrices (tetax,tetay,tetaz)
Rm = Rotate3D(0,0,0);           % Maternal rotation matrices (tetax,tetay,tetaz)
BWf = 0*[.001*sin(2*pi*[1:N]*.14/fs-pi/6) ; .001*sin(2*pi*[1:N]*.14/fs) ; .003*sin(2*pi*[1:N]*.14/fs)];
BWm = 0*[.05*sin(2*pi*[1:N]*.13/fs+pi/7) ; .01*sin(2*pi*[1:N]*.13/fs) ; .01*sin(2*pi*[1:N]*.13/fs)];

%BW = [.15 .35 .15 .2 .5 .11 .4 .3]'*.15*sin(2*pi*[1:N]*.14/fs-pi/6);
%BW = D(:,1)*.01*sin(2*pi*[1:N]*.14/fs - randn(1,N));
BW = H_m(:,1)*.1*sin(2*pi*[1:N]*.14/fs);

VCGm = Rm*diag([1,1,1])*[DIPm.x ; DIPm.y ; DIPm.z];
VCGf = Rf*diag([1,1,1])*[DIPf.x ; DIPf.y ; DIPf.z];

%s = D*[s.x ; s.y ; s.z] + .01*randn(size(D,1),N);
%s = D*[VCGm+BWm] + D0*[VCGf+BWf] + 1e-6*randn(size(D,1),N);

%s = H_m*[VCGm+BWm] + H_f*[VCGf+BWf] + BW + 20e-4*randn(size(H_m,1),N);
%s = H_m*[VCGm+BWm] + H_f*[VCGf+BWf] + BW + filter(1,[1 -1],3e-6*randn(size(H_m,1),N)')' + 1e-6*randn(size(H_m,1),N);
s = H_m*[VCGm+BWm] + H_f*[VCGf+BWf] + BW + filter(1,[1 -1],5e-7*randn(size(H_m,1),N)')' + 0e-7*randn(size(H_m,1),N);

BW2 = Frank_m(:,2)*.001*sin(2*pi*[1:N]*.14/fs);
Frank = Frank_m*[VCGm+BWm] + Frank_f*[VCGf+BWf] + BW2 + filter(1,[1 -1],5e-7*randn(size(Frank_m,1),N)')' + 10e-7*randn(size(Frank_m,1),N);
Frank = 1000*Frank;
%Frank = [VCGm+BWm]
plot3(Frank(1,:),Frank(2,:),Frank(3,:));
grid;

figure;
plot(Frank');
grid;


%//////////////////////////////////////////////////////////////////////////
W =  jadeR(s);
s2 = real(W*s);
A = pinv(W);

%//////////////////////////////////////////////////////////////////////////
L = 1;
close all;
for i = 1:size(s,1),
    %figure(floor((i-1)/L)+1);
    %subplot(L,1,mod(i-1,L)+1);
    h = figure;
    plot(s2(i,:),'b','Linewidth',1);
    xlabel('n');
    ylabel(['IC',num2str(chindex(i))]);
    %title(['Ch',num2str(chindex(i))]);
    grid;

    % % %     set(h,'PaperUnits','inches');
    % % %     %set(h,'Resize','off');
    % % %     %set(h,'PaperPositionMode','auto');
    % % %     set(h,'PaperPosition',[.01 .01 3.5 2])
    % % %     %print('-djpeg',['C:\testplot',num2str(chindex(i)),'.jpg']);
    % % %     print('-dpng',['C:\testplot',num2str(chindex(i)),'.png']);

    %     a = axis;
    %     a(3) = -1;
    %     a(4) = 1;
    %     axis(a);
    %     hold on;
    %    plot(data(i,:));
    %    plot(100*s(1,:)','r');
end

%//////////////////////////////////////////////////////////////////////////



