% test baseline wander removal using median filter
% Reza Sameni, 2005

clear;
close all;
% load LOUETTE_07-Oct-2005_40s.txt -mat;data = LOUETTE_40s;clear LOUETTE_40s
% load fortict_20-Jan-2006.datfortict_20-Jan-2006.dat_acq.mat;
% load fortict4_20-Jan-2006.datfortict4_20-Jan-2006.dat_acq.mat;
%load Bonhomme3_20-Jan-2006.datBonhomme3_20-Jan-2006.dat_acq.mat;
%load patient165_s0323lre;
load FOETAL_ECG.dat; data = FOETAL_ECG(:,2:end)';clear FOETAL_ECG;


%//////////////////////////////////////////////////////////////////////////
%remove undesired channels
%I_ = [7,8,16,18,19,30,31,41,51,52,57,58,59]; % LOUETTE_07-Oct-2005_40s
%I_ = [7,8,15,16,19,31,41,52];
%I_ = [7,8,16,19,41,50];
% I_ = [];
% I = [];
% for i = 1:16,
%     if(isempty(find(I_==i)))
%         I = [I i];
%     end
% end
% J = 10001:15000;
% data = data(J,I);
% data = data';
%//////////////////////////////////////////////////////////////////////////
% denoising
L = 2;
close all;
datalen = length(data);
sdn = zeros(size(data));
sdn2 = zeros(size(data));
sdn4 = zeros(size(data));
bl = zeros(size(data));

xd = zeros(1,datalen);
xd2 = zeros(1,datalen);
xd3 = zeros(1,datalen);

flen0 = 100;
flen = 50;

for i = 1:size(data,1),
    %xd0  = wden(data(i,:),'heursure','s','mln',10,'sym5');
    for j = 1+flen0:datalen-flen0+1,
        %index = max(j-100,1):min(j+100,datalen);
        %xd(j) = median(data(i,index));
        xd(j) = median(data(i,j-flen0:j+flen0-1));
    end
    for j = 1+flen:datalen-flen+1,
        xd2(j) = median(data(i,j-flen:j+flen-1));
    end
    for j = 1:datalen,
        index = max(j-150,1):min(j+149,datalen);
        xd3(j) = median(xd2(index));
    end
%     for j = 1:datalen,
%         index = max(j-600,1):min(j,datalen);
%         xd4(j) = median(xd2(index));
%     end
    sdn(i,:) = data(i,:) - xd;
    sdn2(i,:) = data(i,:) - xd2;
    sdn3(i,:) = data(i,:) - xd3;
    bl(i,:) = xd3;
%     sdn4(i,:) = data(i,:) - xd4;
    % % %     figure(floor((i-1)/L)+1);
    % % %     subplot(L,1,mod(i-1,L)+1);
    % % %     plot(data(i,:));grid;
    % % %     hold on;
    % % %     plot(sdn(i,:),'r');
    % % %     plot(xd,'g');
    % % %     ylabel(num2str(i));
end
% % % % % %//////////////////////////////////////////////////////////////////////////
fs = 500;
L = length(data);

S = fft(data,L,2);

S(1,:) = 0;

k1 = 2:ceil(0.4*L/fs);
S(:,k1) = 0; S(:,L-k1+2) = 0;

% k = floor(488.4*L/fs):ceil(488.6*L/fs); % ADC noise
% S(:,k) = 0; S(:,L-k+2) = 0;
% 
% for i = 1:9,
%     k = floor((i*50-.1)*L/fs):ceil((i*50+.1)*L/fs);
%     S(:,k) = 0; S(:,L-k+2) = 0;
% end
sdn5 = real(ifft(S,L,2));% + 1e-3*randn(size(data));

clear S*;
clear k*;
% % %//////////////////////////////////////////////////////////////////////////
L = 1;
close all;
tag = ['A','B','C','D','E','F','G','H'];
for i = 1:size(data,1),
    figure(floor((i-1)/L)+1);
    subplot(L,1,mod(i-1,L)+1);
    plot(data(i,:));grid;
    hold on;
    plot(sdn(i,:),'r');
    plot(sdn2(i,:),'g');
    plot(sdn3(i,:),'m');
    plot(sdn5(i,:),'k');
    plot(bl(i,:),'c');
    title(['Channel: ',num2str(i),' (',tag(floor((i-1)/9)+1),num2str(mod(i-1,9)+1),')']);
    %legend('Original','Median','Mean','Lowpass','MMedian');
    legend('Original','Median 200','Median 100','Median 100-300','Lowpass','baseline');
end
% % % % %//////////////////////////////////////////////////////////////////////////
% % Second ICA
%
% %[s2, A2, W2] = fastica(xbar,'initGuess',A,'approach', 'symm','displayMode', 'off');
% [s2, A2, W2] = fastica(xbar,'approach', 'symm','displayMode', 'off');
% [s2, A2, W2] = fastica([sdn -sdn],'approach', 'symm','displayMode', 'off');%,'lastEig', 30);
% % % % % %
% % % % % % % [E, D] = pcamat(data');
% % % % % % % dat = data*E;
% % % L = 3;
% % % close all;
% % % %for i = 1:size(data,2),
% % % for i = 1:size(s2,1),
% % %     figure(floor((i-1)/L)+1);
% % %     subplot(L,1,mod(i-1,L)+1);
% % %     %plot(s2(i,:));grid;
% % %     plot(s2(i,1:end/2));grid;
% % %     title(['ICA Extracted Channel Number: ',num2str(i)]);
% % % end

% save ProcessData3_
