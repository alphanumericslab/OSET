% test plot adult ECG and VCG loops
% Reza Sameni, 2005

clear;
close all;

% % % fs = 1000;
% % % data = load('s0010_re.txt');
% % % % data = load('s0014lre.txt');
% % % time = data(:,1)';
% % % data = data(:,2:end)';
% % % I = 1:10*fs;
% % % data = data(:,I);

fs = 1000;
load patient165_s0323lre; data = data(101:50100,2:end)';

data = data(:,1:10*fs);

time = (0:length(data)-1)/fs;


% % % %//////////////////////////////////////////////////////////////////////////
% % % % baseline removal
% % % sdn = zeros(size(data));
% % % for i = 1:size(data,1),
% % %     sdn(i,:) = data(i,:) - LPFilter(Median(data(i,:),length(data),200,400)',50/fs);
% % % end
% % % %//////////////////////////////////////////////////////////////////////////
W =  jadeR(data,8);
A = pinv(W);
s = W*data;
%//////////////////////////////////////////////////////////////////////////

ss = s;
ss(1,:) = ss(1,:) - max(ss(1,:));
for i = 2:size(ss,1)
    ss(i,:) = ss(i,:) + min(ss(i-1,:)) - max(ss(i,:));
end


% I = 511:1110;
I = 1:length(time);
figure
plot(time(I)-time(I(1)),ss(:,I)','linewidth',1);
grid
set(gca,'Box','On','FontSize',16);
set(gca,'YTickLabel',[]);
xlabel('Time (s)','FontSize',16);
ylabel('Components');
zoom reset;
mn = mean(ss,2);
set(gca,'YTick',mn(end:-1:1));
set(gca,'YTickLabel',num2str((size(ss,1):-1:1)'));

% axis tight;
% % % L = size(s,1);
% % % for i = 1:L,
% % %     subplot(L,1,i);
% % %     plot(time(I)-time(I(1)),s(i,I)','linewidth',1);
% % %     set(gca,'Box','On','FontSize',16);
% % % %     axis tight;
% % %     grid
% % %
% % %     if (mod(i,L)~=0)
% % %         set(gca,'XTickLabel',[]);
% % %     else
% % %         xlabel('Time (s)','FontSize',16);
% % %     end
% % % end
% % % set(gcf,'Position',[438 350 860 420]);

mn = mean(data,2);

figure
plot3(data(13,:),data(14,:),data(15,:),'b.');
set(gca,'Box','On','FontSize',16);
grid;
hold on;
A = A*10;
for k = 1:8
    quiver3(mn(13),mn(14),mn(15),A(13,k),A(14,k),A(15,k),'g','LineWidth',2);
    quiver3(mn(13),mn(14),mn(15),-A(13,k),-A(14,k),-A(15,k),'g','LineWidth',2);
end
view([157 20]);
xlabel('X','FontSize',16);
ylabel('Y','FontSize',16);
zlabel('Z','FontSize',16);

% % % plot(window,'g','LineWidth',2);
% % % grid
% % % figure
% % % plot(s2');
% % % hold on;
% % % plot(window,'g','LineWidth',2);
% % % grid
% % % L = 1;
% % % close all;
% % % for i = 1:size(data,1),
% % %     figure(floor((i-1)/L)+1);
% % %     subplot(L,1,mod(i-1,L)+1);
% % %     %    if(i<7),plot(s(i,:));elseif(i==7),plot(s(8,:));elseif(i==8),plot(s(7,:));end
% % %     plot(sdn(i,:));
% % %     %     hold on;
% % %     %     plot(s2(i,:),'r');
% % %     %     plot(window,'g');
% % %     grid;
% % % end
%//////////////////////////////////////////////////////////////////////////
