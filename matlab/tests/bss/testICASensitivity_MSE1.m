% Test semsitibity of different ECG segments to ICA dim reduction
% Reza Sameni, Copyright 2008

clear;
close all;

fs = 1000;
data = load('s0010_re.txt');
% data = load('s0014lre.txt');


time = data(:,1)';
data = data(:,2:end)';

data = data - LPFilter(data,1.5/fs);
data = LPFilter(data,150/fs);

I = 1:10*fs;
data = data(:,I);

N = length(data);

% % % mn = mean(data,2)*ones(1,N);
% % % data = data - mn;
% % % sd = std(data,[],2)*ones(1,N);
% % % data = data./sd;
% % % clear mn sd;
%//////////////////////////////////////////////////////////////////////////
refchannel = data(15,:);
peaks = PeakDetection(refchannel,1/fs);


QRSstart = -60;
QRSstop = 65;

Tstart = 80;
Tstop = 345;

Pstart = -220;
Pstop = -90;

window = zeros(1,N);
window1 = zeros(1,N);
window2 = zeros(1,N);
window3 = zeros(1,N);

II = find(peaks);
for i = II
    window(max(i+Pstart,1):min(i+Tstop,N)) = 1;     % ECG waves
    window1(max(i+QRSstart,1):min(i+QRSstop,N)) = 1;     % QRS waves
    window2(max(i+Tstart,1):min(i+Tstop,N)) = 1;           % T waves
    window3(max(i+Pstart,1):max(i+Pstop,1)) = 1;          % P waves
end

I = find(window);    %   AllECG indexes
I1 = find(window1);  %   QRS indexes
I2 = find(window2);  %   T indexes
I3 = find(window3);  %   QRS indexes

indx = 1:15;
L = length(indx);
% % % SNR = zeros(size(data,1),L);
% % % SNR1 = zeros(size(data,1),L);
% % % SNR2 = zeros(size(data,1),L);
% % % SNR3 = zeros(size(data,1),L);
SNR = zeros(1,L);
SNR1 = zeros(1,L);
SNR2 = zeros(1,L);
SNR3 = zeros(1,L);
%//////////////////////////////////////////////////////////////////////////
for i = 1:L
    W = jadeR(data);
    s = W*data;
    s(end-i+1:end,:) = 0;
    data2 = pinv(W)*s;

    d = data(:,I); d = d(:);
    d2 = data2(:,I); d2 = d2(:);
    P = std(d);
    N = std(d-d2);
    SNR(i) = 20*log10(N./P);

    d = data(:,I1); d = d(:);
    d2 = data2(:,I1); d2 = d2(:);
    P = std(d);
    N = std(d-d2);
    SNR1(i) = 20*log10(N./P);

    d = data(:,I2); d = d(:);
    d2 = data2(:,I2); d2 = d2(:);
    P = std(d);
    N = std(d-d2);
    SNR2(i) = 20*log10(N./P);

    d = data(:,I3); d = d(:);
    d2 = data2(:,I3); d2 = d2(:);
    P = std(d);
    N = std(d-d2);
    SNR3(i) = 20*log10(N./P);

% % %     P = std(data,[],2);
% % %     N = std(data-data2,[],2);
% % %     SNR(:,i) = 20*log10(P./N);
% % % 
% % %     P1 = std(data(:,I1),[],2);
% % %     N1 = std(data(:,I1)-data2(:,I1),[],2);
% % %     SNR1(:,i) = 20*log10(P1./N1);
% % % 
% % %     P2 = std(data(:,I2),[],2);
% % %     N2 = std(data(:,I2)-data2(:,I2),[],2);
% % %     SNR2(:,i) = 20*log10(P2./N2);
% % % 
% % %     P3 = std(data(:,I3),[],2);
% % %     N3 = std(data(:,I3)-data2(:,I3),[],2);
% % %     SNR3(:,i) = 20*log10(P3./N3);
end
% % % SNR = mean(SNR,1);
% % % SNR1 = mean(SNR1,1);
% % % SNR2 = mean(SNR2,1);
% % % SNR3 = mean(SNR3,1);

shift = 60;
I = 1:L-1;
figure;
bar(indx(I),[SNR(I) ; SNR1(I) ; SNR2(I) ; SNR3(I)]' + shift,'BarWidth',1);
legend('Total ECG','QRS','T','P');
grid;
axis tight
xlabel('Number of Dimensions Removed','FontSize',16);
ylabel('Normalized MSE (dB)','FontSize',16);
set(gca,'Box','On','FontSize',16);
set(gcf,'Position',[438 350 890 420]);
tic = get(gca,'YTick');
set(gca,'YTickLabel',tic-shift);
colormap(summer(4));

% % % figure;
% % % plot(SNR,indx,'m');
% % % hold on;
% % % plot(SNR1,indx,'b');
% % % plot(SNR2,indx,'r');
% % % plot(SNR3,indx,'g');
% % % legend('Total','QRS','T','P');
% % % grid;


% % % figure;
% % % plot(data','b');
% % % hold on;
% % % plot(refchannel.*window1,'go');
% % % plot(refchannel.*window2,'mo');
% % % plot(refchannel.*window3,'ro');
% % % plot(refchannel.*peaks,'ko');
% % % grid

%//////////////////////////////////////////////////////////////////////////
% L = 1;
% close all;
% for i = 1:size(data,1),
%     figure(floor((i-1)/L)+1);
%     subplot(L,1,mod(i-1,L)+1);
%     plot(data(i,:));grid;
%     hold on;
%     plot(data2(i,:),'g');
%     plot(window1,'r');
%     plot(window2,'m');
%     plot(window3,'c');
% end

