clc
clear;
close all;
% addpath('C:\Reza\work\ECGKalmanFilter\ToolBox');
% addpath('C:\Reza\work\MILCA');
%//////////////////////////////////////////////////////////////////////////
load LOUETTE_07-Oct-2005_30s.txt -mat; data = LOUETTE1(1:1200,:)';clear LOUETTE1
%//////////////////////////////////////////////////////////////////////////

% maternal abdomen as a spheroid

scale = 30;
a = 1.2*scale;
b = 1*scale;
c = 1*scale;

%phi = [0:.2:pi pi]';
phi = [pi/2.5:.2:pi-pi/2.5]';
step = (2*pi)/18;
teta = 0:step:2*pi-step;

x = a*sin(phi)*cos(teta);
y = b*sin(phi)*sin(teta);
z = c*cos(phi)*ones(size(teta));

x = -x;
y = -y;

S = [x(:) y(:) z(:)]';

navalindex = 39;

naval = S(:,navalindex);

% S = [S(:,1:navalindex-1) S(:,navalindex+1:end)];
% x = [x(:,1:navalindex-1) x(:,navalindex+1:end)];
% y = [y(:,1:navalindex-1) y(:,navalindex+1:end)];
% z = [z(:,1:navalindex-1) z(:,navalindex+1:end)];

S = S - naval*ones(1,size(S,2));
x = x - naval(1);
y = y - naval(2);
z = z - naval(3);
figure;
h = mesh(x,y,z);

% % figure;
% % hh = mesh(x,y,z,[MI [0 0 0]']);


set(h,'LineWidth',2);
colormap([0 0 1]);
axis equal
hold on;
plot3(S(1,:),S(2,:),S(3,:),'ks','LineWidth',2);
plot3(0,0,0,'ro','LineWidth',3);
xlabel('X');
ylabel('Y');
zlabel('Z');

%//////////////////////////////////////////////////////////////////////////
km = 1;
kf = .3;

%elec = [-2 5 5; -2 -5 5; -2 -5 -5; -2 5 -5 ; -6 5 5; -6 -5 5; -6 -5 -5; -6 5 -5];% ; -25 7 20];
%elec = [-6 -10 2;-6 -10 -2;-3 -5 2;-3 -5 -2;-3 5 2;-3 5 -2;-6 10 2;-6 10 -2];

elec = S';

% rf = find(S(1,:)==0 & S(2,:)==0 & S(3,:)==0);
% elec = [S(:,1:rf-1) S(:,rf+1:end)]';

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

% FrankPos = [-5 -8 -2; -3 0 -1; -3 -5 -5];
% FrankNeg = [-10 -8 -2; -3 -15 -1; -3 -5 5];
FrankPos = [-10 20 20; -20 20 20; -30 20 30];       % Maternal Frank leads
FrankNeg = [-60 20 20; -20 -20 20; -30 20 10];     % Maternal Frank leads
for i = 1:3,
    for j = 1:3,
        Frank_m(i,j) = km* ((FrankPos(i,j)-pos_m(j))/sqrt(sum((FrankPos(i,:)-pos_m).^2))^3 - (FrankNeg(i,j)-pos_m(j))/sqrt(sum((FrankNeg(i,:)-pos_m).^2))^3);
        Frank_f(i,j) = kf* ((FrankPos(i,j)-pos_f(j))/sqrt(sum((FrankPos(i,:)-pos_f).^2))^3 - (FrankNeg(i,j)-pos_f(j))/sqrt(sum((FrankNeg(i,:)-pos_f).^2))^3);
    end
end

%//////////////////////////////////////////////////////////////////////////
N = 500;
fs = 250;

% mother params
F_m = .9;
teta0_m = -pi/2;

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
% F_m = .9;
% teta0_m = pi;

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
F_f = 2.2;
teta0_f = -pi/2.5;

tetai_f.x  = [-0.7    -0.17    0       0.18     1.4];
alphai_f.x = .07*[0.07     -0.11   1.3     0.07   0.275];
%alphai_f.x = .2*[0.07     -0.11   1.3     0.07   0.275];
bi_f.x     = [.1       .03     .045     0.02    0.3];

tetai_f.y  = [-0.9     -0.08   0       0.05        1.3];
alphai_f.y = .08*[0.04     0.3     .45     -0.35       0.09];
%alphai_f.y = .3*[0.04     0.3     .45     -0.35       0.09];
bi_f.y     = [.1       .05      .03    .04         .3];

tetai_f.z  = [-0.8      -.3     -0.1        .06     1.35];
alphai_f.z = .05*[-0.14    .03     -0.4        .46     -0.2];
%alphai_f.z = .4*[-0.14    .03     -0.4        .46     -0.2];
bi_f.z     = [.1       .4      .03         .03     .3];


%//////////////////////////////////////////////////////////////////////////
[DIPm tetam] = DipoleGenerator(N,fs,F_m,alphai_m,bi_m,tetai_m,teta0_m);
[DIPf tetaf] = DipoleGenerator(N,fs,F_f,alphai_f,bi_f,tetai_f,teta0_f);
% s = H_f*[DIPm.x;DIPm.y;DIPm.z];


%Rf = Rotate3D(pi,pi/7,pi/2);   % Fetal rotation matrices (tetax,tetay,tetaz)
%Rf = Rotate3D(pi/8,pi,pi/2);   % Fetal rotation matrices (tetax,tetay,tetaz)
Rf = Rotate3D(pi,0,-pi/2);   % Fetal rotation matrices (tetax,tetay,tetaz)
Rm = Rotate3D(0,0,0);           % Maternal rotation matrices (tetax,tetay,tetaz)
BWf = 1*[.001*sin(2*pi*[1:N]*.14/fs-pi/6) ; .001*sin(2*pi*[1:N]*.14/fs) ; .003*sin(2*pi*[1:N]*.14/fs)];
BWm = 1*[.1*sin(2*pi*[1:N]*.13/fs+pi/7) ; .01*sin(2*pi*[1:N]*.13/fs) ; .01*sin(2*pi*[1:N]*.13/fs)];

%BW = [.15 .35 .15 .2 .5 .11 .4 .3]'*.15*sin(2*pi*[1:N]*.14/fs-pi/6);
%BW = D(:,1)*.01*sin(2*pi*[1:N]*.14/fs - randn(1,N));
BW = H_m(:,1)*.001*sin(2*pi*[1:N]*.14/fs);

VCGm = Rm*diag([1,1,1])*[DIPm.x ; DIPm.y ; DIPm.z];
VCGf = Rf*diag([1,1,1])*[DIPf.x ; DIPf.y ; DIPf.z];

%s = D*[s.x ; s.y ; s.z] + .01*randn(size(D,1),N);
%s = D*[VCGm+BWm] + D0*[VCGf+BWf] + 1e-6*randn(size(D,1),N);

%s = H_m*[VCGm+BWm] + H_f*[VCGf+BWf] + BW + 20e-4*randn(size(H_m,1),N);
%s = H_m*[VCGm+BWm] + H_f*[VCGf+BWf] + BW + filter(1,[1 -1],5e-7*randn(size(H_m,1),N)')' + 1e-7*randn(size(H_m,1),N);

s = H_m*[VCGm+BWm] + H_f*[VCGf+BWf] + BW + filter(1,[1 -1],1e-6*randn(size(H_m,1),N)')' + 1e-5*randn(size(H_m,1),N);

s_mref = H_m*VCGm;
s_fref = H_f*VCGf;


% s = 1e-7*rand(size(H_m,1),N);
% s = filter(1,[1 -1],5e-6*randn(size(H_m,1),N)')';

% Num = 0.0001*[1.0000    2.0000    1.0000];
% Den = [1.0000   -1.9684    0.9688];

% Num = 1;
% Den = [1 -1 1];
% s = filter(Num,Den,5e-6*randn(size(H_m,1),N)')';
% s = s - mean(s,2)*ones(1,size(s,2));
% s = data;


BW2 = Frank_m(:,2)*.001*sin(2*pi*[1:N]*.14/fs);
Frank = Frank_m*[VCGm+BWm] + Frank_f*[VCGf+BWf] + BW2 + filter(1,[1 -1],1e-6*randn(size(Frank_m,1),N)')' + 2e-5*randn(size(Frank_m,1),N);
% Frank = 1000*Frank;
%Frank = [VCGm+BWm]

% figure;
% plot3(Frank(1,:),Frank(2,:),Frank(3,:));
% grid;
% 
figure;
plot(Frank');
grid;


figure;
plot(s');
grid;

%//////////////////////////////////////////////////////////////////////////
p = randperm(size(s,2));
[yy,pp] = sort(p);
s = s(:,p);

W =  jadeR(s,10);
s2 = real(W*s);
A = pinv(W);

s = s(:,pp);
s2 = s2(:,pp);

% [s2, A, W] = fastica(s,'displayMode', 'off','lastEig', 7);%,'numOfIC',30);

% % % %//////////////////////////////////////////////////////////////////////////
plot_multichannel_data(s2,2,'b');
% % % L = 1;
% % % close all;
% % % for i = 1:size(s,1),
% % %     %figure(floor((i-1)/L)+1);
% % %     %subplot(L,1,mod(i-1,L)+1);
% % %     h = figure;
% % %     plot(s2(i,:),'b','Linewidth',1);
% % %     xlabel('n');
% % %     ylabel(['IC',num2str(chindex(i))]);
% % %     %title(['Ch',num2str(chindex(i))]);
% % %     grid;
% % % 
% % %     % % %     set(h,'PaperUnits','inches');
% % %     % % %     %set(h,'Resize','off');
% % %     % % %     %set(h,'PaperPositionMode','auto');
% % %     % % %     set(h,'PaperPosition',[.01 .01 3.5 2])
% % %     % % %     %print('-djpeg',['C:\testplot',num2str(chindex(i)),'.jpg']);
% % %     % % %     print('-dpng',['C:\testplot',num2str(chindex(i)),'.png']);
% % % 
% % %     %     a = axis;
% % %     %     a(3) = -1;
% % %     %     a(4) = 1;
% % %     %     axis(a);
% % %     %     hold on;
% % %     %    plot(data(i,:));
% % %     %    plot(100*s(1,:)','r');
% % % end
% % % 
% % % %//////////////////////////////////////////////////////////////////////////

figure
plot(s2','b');

% hold on;
% plot(1000*s','r');

% y = filter(1,[1 -1],randn(length(s),1));
% %y = filter([0 1],[-1 1],y);
% figure;
% plot(y');
% grid;
% 

% figure
% psd(s2(3,:),length(s));
% hold on;
% psd(s(4,:),length(s));
% psd(s2(5,:),length(s));
% psd(s(5,:),length(s));

MIPerRefM = zeros(1,size(s,1));
MIPerRefF = zeros(1,size(s,1));
MI3D = zeros(1,size(s,1));
MI = zeros(3,size(s,1));
KN = 10;
for i = 1:size(s,1);
    i
    
    normalized = (s(i,:)-mean(s(i,:)))/std(s(i,:));
    normalizedref = (Frank - mean(Frank,2)*ones(1,length(Frank)))./(std(Frank,1,2)*ones(1,length(Frank)));
    
    mnormalizedref = (s_mref(i,:)-mean(s_mref(i,:)))/std(s_mref(i,:));
    fnormalizedref = (s_fref(i,:)-mean(s_fref(i,:)))/std(s_fref(i,:));
    
    MI3D(i) = MIxnyn(normalizedref,normalized,KN);
    MIPerRefM(i) = MIxnyn(mnormalizedref,normalized,KN);
    MIPerRefF(i) = MIxnyn(fnormalizedref,normalized,KN);
    MI(1,i) = MIxnyn(normalizedref(1,:),normalized,KN);
    MI(2,i) = MIxnyn(normalizedref(2,:),normalized,KN);
    MI(3,i) = MIxnyn(normalizedref(3,:),normalized,KN);
end
MI3D(navalindex) = min(MI3D);
MIPerRefM(navalindex) = min(MIPerRefM);
MIPerRefF(navalindex) = min(MIPerRefF);
MI(1,navalindex) = min(MI(1,:));
MI(2,navalindex) = min(MI(2,:));
MI(3,navalindex) = min(MI(3,:));


figure;
plot(MI3D);
grid;

hold on
plot(MI','r');

figure;
plot(MI3D);
grid;

hold on
plot(MI','r');

figure;
plot(MIPerRefM,'b');
hold on
plot(MIPerRefF,'r');
grid;

figure;
h = surf(x,y,z,reshape(MIPerRefF,size(z)));
set(h,'LineWidth',2);
%colormap([0 0 1]);
shading interp
colormap jet
axis equal
hold on;
plot3(S(1,:),S(2,:),S(3,:),'ks','LineWidth',2);
plot3(0,0,0,'ro','LineWidth',3);
xlabel('X');
ylabel('Y');
zlabel('Z');






