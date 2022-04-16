% Simulated scatter plots and ICA
%
clear;
close all;

a = 2;
b = 2;
% t = asin(1/b) + [0:.01:18*pi 0];
% t = asin(1/b) + [0:.01:2*pi 0];
t = 0:.01:20*pi;

N = length(t);
r = a * (1-b*sin(t));% + .2*randn(1,length(t));
[x,y] = pol2cart(t,r);
z = 2*cos(t);
% z = t/10;

A0 = [-1 1 -2; 1 1 1; -.2 1 1];
xyz = A0*[x;y;z];

data = zeros(size(xyz));

data(1,:) = xyz(1,:) + .05*LPFilter(cumsum(randn(size(t))),.001);
data(2,:) = xyz(2,:) + .05*LPFilter(cumsum(randn(size(t))),.001);
data(3,:) = xyz(3,:) + .05*LPFilter(cumsum(randn(size(t))),.001);

% teta = [0:.01:16*pi];
% x = .5*cos(teta).^3+ .1*(rand(1,length(teta))-.5);
% y = .6*sin(teta).^3+ .2*(rand(1,length(teta))-.5);

% x = 15*cos(teta).^3+ .5*randn(1,length(teta));
% y = 16*sin(teta).^3+ .7*randn(1,length(teta));

% x = 15*(cos(teta).^2-.5)+ 2*randn(1,length(teta));
% y = 16*(sin(teta).^2-.5)+ 1*randn(1,length(teta));

% x = 15*cos(teta)+ 2*randn(1,length(teta));
% y = 16*sin(teta)+ 1*randn(1,length(teta));
% z = 0*teta;

% xyz = xyz(:,1:2:end);
% % % npts = length(x);
% % % figure;
% % % plot3(xyz(1,:),xyz(2,:),xyz(3,:),'ro','LineWidth',2);
% % % % text(xyz(1,:),xyz(2,:),xyz(3,:),[repmat('  ',npts,1), num2str([1:npts]')]);
% % % set(gca,'XTick',[],'YTick',[],'ZTick',[]);
% % % box on;
% % % hold on;
% % % fnplt(cscvn(xyz(:,[1:end 1])),'r',2);
% % % grid on;
% % % % plot3(x,1:npts,z*0-1,'g');
% % % % plot3(1:npts,y,z*0-1,'b');
% % % hold off;
% % % % figure
% % % % plot(t,r);
% % % % grid;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I = 20000:25000;
%dd = data(:,2:end)';

[s, A, W] = fastica(data,'approach', 'symm','displayMode', 'off','epsilon',.001);
%[s, A, W] = fastica(dd,'displayMode', 'off');


% % % I = 511:1110;
% % % figure
% % % axis tight;
% % % plot(time(I)-time(I(1)),s(:,I)','linewidth',2);
% % % set(gca,'Box','On','FontSize',16);
% % % xlabel('time(s)','FontSize',16);
% % % ylabel('Amplitude','FontSize',16);
% % % grid
% % % 
% % % figure
% % % plot3(sdn(13,:),sdn(14,:),sdn(15,:),'b.');
% % % set(gca,'Box','On','FontSize',16);
% % % grid;
% % % hold on;
% % % A = A*10;
% % % for k = 1:5,
% % %     quiver3(0,0,0,A(13,k),A(14,k),A(15,k),'g','LineWidth',2);
% % %     quiver3(0,0,0,-A(13,k),-A(14,k),-A(15,k),'g','LineWidth',2);
% % % end
% % % view([157 20]);
% % % xlabel('X','FontSize',16);
% % % ylabel('Y','FontSize',16);
% % % zlabel('Z','FontSize',16);


%plot3(data(I,14),data(I,15),data(I,16));

% % % figure
% % % plot3(s(1,:),s(2,:),s(3,:),'b.','LineWidth',2);
% % % set(gca,'Box','On','FontSize',16);

mn = mean(data,2);

figure
plot3(data(1,:),data(2,:),data(3,:),'b.','LineWidth',2);
set(gca,'Box','On','FontSize',16);



grid;
hold on;
A = A*2;
for k = 1:3
    quiver3(mn(1),mn(2),mn(3),A(1,k),A(2,k),A(3,k),2,'r','LineWidth',3);
    quiver3(mn(1),mn(2),mn(3),-A(1,k),-A(2,k),-A(3,k),2,'r','LineWidth',3);
end
view([-160 -20]);
xlabel('X','FontSize',16);
ylabel('Y','FontSize',16);
zlabel('Z','FontSize',16);

% % % for i = 0:10:90;
% % %     for j = -180:1:180;
% % %         view(j,i);
% % %         pause(.001);
% % %     end
% % % end

I = 1:5000;
figure;
plot(data(:,I)','linewidth',2);
grid;
set(gca,'Box','On','FontSize',16);
xlabel('Samples','FontSize',16);
ylabel('Amplitude','FontSize',16);

I = 1:4000;
figure;
plot(s(:,I)','linewidth',2);
grid;
set(gca,'Box','On','FontSize',16);
xlabel('Samples','FontSize',16);
ylabel('Amplitude','FontSize',16);
