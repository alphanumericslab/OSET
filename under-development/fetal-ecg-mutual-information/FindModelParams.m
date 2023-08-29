clear;
close all;
load patient165_s0323lre;
time = data(:,1)';
data = data(:,14:16)';

N = length(data);
fs = 1000;
Len = 1000;
%//////////////////////////////////////////////////////////////////////////
S = fft(data,N,2);
f = [0:N-1]*fs/N;
S(:,1) = 0;
k1 = 2:ceil(.5*N/fs);
S(:,k1) = 0; S(:,N-k1+2) = 0;
k2 = floor(150*N/fs):ceil(500*N/fs);
S(:,k2) = 0; S(:,N-k2+2) = 0;
data = real(ifft(S,N,2));
%//////////////////////////////////////////////////////////////////////////
F_m = fs*58 /(59248-839)
%//////////////////////////////////////////////////////////////////////////
% Peak search
phase = ones(N,1);
peaks = zeros(N,1);
cntr = 0;
for i = 1:N,
    if(data(1,i)>1 & max(data(1,max(i-5,1):min(i+5,N)))==data(1,i))
        peaks(i) = 1;
        cntr = cntr + 1;
    end
end
I = find(peaks==1);
AllECG = zeros(length(I)-1,size(data,1),Len);
for i = 1:length(I)-1;
    m = I(i+1) - I(i);
    phase(I(i)+1:I(i+1)-1) = 2*pi/m : 2*pi/m : 2*pi-2*pi/m;
    ind = max(I(i)-Len/2+1,1):min(I(i)+Len/2,N);
    AllECG(i,:,end/2-length(ind)/2+1:end/2+length(ind)/2) = data(:,ind);
end
phase(1:I(1)) = 2*pi/I(1):2*pi/I(1):2*pi;
m = length(phase(I(end)+1:end));
phase(I(end)+1:end) = 2*pi/m:2*pi/m:2*pi;
phase = mod(phase,2*pi);
I = find(phase>pi);
phase(I) = phase(I) - 2*pi;

%//////////////////////////////////////////////////////////////////////////
mnECG = squeeze(mean(AllECG,1));
bs = mean(mnECG(:,200:300),2)*ones(1,length(mnECG));
% bs = mean(mnECG(:,300:400),2)*ones(1,length(mnECG));
mnECG = mnECG-bs;

Nrr = fs/F_m;
ph = 2*pi/Nrr : 2*pi/Nrr : Len/Nrr*2*pi;
ph = ph - ph(end/2); % R-peak phase;
phase = mod(ph,2*pi);
I = find(ph>pi);
ph(I) = ph(I) - 2*pi;
I = find(ph<-pi);
ph(I) = ph(I) + 2*pi;

%J = 200:1100;
plot(ph,mnECG');
grid


% % % %time = [0:N-1]/fs;
% mother params
F_m = F_m; % already calculated
teta0_m = ph(1);

%//////////////////////////////////////////////////////////////////////////
tetai_m.x  = [-1.09  -0.83   -0.19     -.07  0 .06        0.22    1.2 1.42 1.68 2.9];
alphai_m.x = [0.03   .08    -0.13    .85 1.11 .75     0.06   0.1  0.17 0.39 .03];
bi_m.x     = [0.0906    0.1057    0.0453    0.0378    0.0332    0.0302    0.0378    0.6040 0.3020    0.1812 .5];
% OK

tetai_m.y  = [-1.1  -0.9 -0.76       -0.11   -.01       0.065  0.8      1.58 2.9];
alphai_m.y = [0.035 0.015 -0.019     0.32    .51     -0.32    0.04   0.08 .014];
bi_m.y     = [0.07  .07  0.04        0.055    0.037    0.0604  0.450  0.3 .5];
% OK

tetai_m.z  = [-1.1  -0.93 -0.7      -.4     -0.15    .095    1.05 1.25 1.55 2.8];
alphai_m.z = [-0.03 -0.14 -0.035    .045     -0.4    .46    -.12 -.2 -.35 -.035];
bi_m.z     = [.03  .12  .04         .4    .045       .05    .8 .4 .2 .4];

teta_iso = -pi;
[DIPm tetam] = DipoleGenerator(Len,fs,F_m,alphai_m,bi_m,tetai_m,teta0_m,teta_iso);
s = [DIPm.x;DIPm.y;DIPm.z];
% H = data(1:5,I)/s;
% s2 = H*s;
%//////////////////////////////////////////////////////////////////////////
L = 3;
close all;
%for i = 1:size(data,1),
for i = [1:3],
    figure(floor((i-1)/L)+1);
    subplot(1,L,mod(i-1,L)+1);
%     plot(s2(i,:),'b');
    plot(ph,mnECG(i,:),'b');
    hold on;
    plot(ph,s(i,:)','r');
    grid;
end
legend('Original ECG','Synthetic ECG')
L = 3;
close all;
for i = [1:3],
    figure(floor((i-1)/L)+1);
    subplot(1,L,mod(i-1,L)+1);
    plot(ph,100*abs(mnECG(i,:)-s(i,:)),'b');
    grid;
end
I = 20:980;
figure;
plot3(s(1,I),s(2,I),s(3,I));
xlabel('X');
ylabel('Y');
zlabel('Z');
grid;

%//////////////////////////////////////////////////////////////////////////
teta_iso = -1.4;
teta0_m = -pi;
[DIPm tetam] = DipoleGenerator(10000,fs,F_m,alphai_m,bi_m,tetai_m,teta0_m,teta_iso);
s = [DIPm.x;DIPm.y;DIPm.z];
figure
plot(s');
grid;
% hold on;
% plot(tetam,'m');
