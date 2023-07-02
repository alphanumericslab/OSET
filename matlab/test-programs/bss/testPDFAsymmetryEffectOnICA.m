% test scatter plot assymetry effect on ICA
% Reza Sameni, Copyright 2005
%
clear;
close all;

t = [0:.001:16*pi 0];

N = length(t);
%r = ceil(t/pi/2).*sin(2.001*t).^3;% + sin(5.00*t).^5;% + .01*rand(size(t)); % + sin(t)
r = sin(2.00*t).^4;% + sin(5.00*t).^5;% + .01*rand(size(t)); % + sin(t)

[x,y] = pol2cart(t,r);

d = [x ; y];
d = d + .04*randn(size(d));
d(2,:) = abs(d(2,:));
% I = find(d(2,:)>0);
% d = d(:,I);

mn = mean(d,2);

%[s, A, W] = fastica([d],'approach', 'symm','displayMode', 'off','epsilon',.0001);
%[s, A, W] = fastica([d],'displayMode', 'off');
W =  jadeR(d);
s = real(W*d);
A = pinv(W);


figure
plot(d(1,:),d(2,:),'b.','LineWidth',1);
axis([-1 1 -1 1]);
grid;
hold on;
for k = 1:2
    quiver(mn(1),mn(2),A(1,k),A(2,k),2,'k','LineWidth',3);
    quiver(mn(1),mn(2),-A(1,k),-A(2,k),2,'k','LineWidth',3);
end
plot(mn(1),mn(2),'ro','LineWidth',5);
set(gca,'Box','On','FontSize',16);
xlabel('x_1','FontSize',16);
ylabel('x_2','FontSize',16);
axis square

W =  jadeR([d -d]);
s = real(W*d);
A = pinv(W);


figure
plot(d(1,:),d(2,:),'b.','LineWidth',1);
hold on;
plot(-d(1,:),-d(2,:),'g.','LineWidth',1);
axis([-1 1 -1 1]);
grid;
for k = 1:2
    quiver(0,0,A(1,k),A(2,k),2,'k','LineWidth',3);
    quiver(0,0,-A(1,k),-A(2,k),2,'k','LineWidth',3);
end
plot(0,0,'ro','LineWidth',5);
set(gca,'Box','On','FontSize',16);
xlabel('x_1','FontSize',16);
ylabel('x_2','FontSize',16);
axis square

