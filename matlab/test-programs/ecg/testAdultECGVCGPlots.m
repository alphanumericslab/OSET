% Test Adults ECG/VCG scatter plots
%
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
time = (0:length(data)-1)/fs;





%//////////////////////////////////////////////////////////////////////////
% baseline removal
sdn = data - LPFilter(data, 50/fs);
% sdn = zeros(size(data));
% for i = 1:size(data,1),
%     sdn(i,:) = data(i,:) - LPFilter(Median(data(i,:),length(data),200,400)',50/fs);
% end

%//////////////////////////////////////////////////////////////////////////
W =  jadeR(sdn,6);
A = pinv(W);
s = W*sdn;
%//////////////////////////////////////////////////////////////////////////
I = 511:1110;
figure
axis tight;
plot(time(I)-time(I(1)),s(:,I)','linewidth',2);
set(gca,'Box','On','FontSize',16);
xlabel('time(s)','FontSize',16);
ylabel('Amplitude','FontSize',16);
grid

figure
plot3(sdn(13,:),sdn(14,:),sdn(15,:),'b.');
set(gca,'Box','On','FontSize',16);
grid;
hold on;
A = A*10;
for k = 1:5
    quiver3(0,0,0,A(13,k),A(14,k),A(15,k),'g','LineWidth',2);
    quiver3(0,0,0,-A(13,k),-A(14,k),-A(15,k),'g','LineWidth',2);
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
