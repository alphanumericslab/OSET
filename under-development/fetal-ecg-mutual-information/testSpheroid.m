clear;
close all;

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

S = [x(:) y(:) z(:)]';

naval = S(:,end/2+3);

S = S - naval*ones(1,size(S,2));
x = x - naval(1);
y = y - naval(2);
z = z - naval(3);

% I = find(S(3,:)<10);
% S = S(:,I);
% I = find(S(3,:)>-8);
% S = S(:,I);

% I = find(z<10);
% [k1,k2] = ind2sub(size(z),I);
% z = z(k1,k2);
% y = y(k1,k2);
% x = x(k1,k2);

% I = find(z>-8);
% [k1,k2] = ind2sub(size(z),I);
% z = z(k1,k2);
% y = y(k1,k2);
% x = x(k1,k2);


plot3(S(1,:),S(2,:),S(3,:));
grid
axis equal
%axis square

% C = hadamard(2^k);
figure;
h = mesh(x,y,z);
set(h,'LineWidth',2);
% alpha(.7);
%shading interp
%colormap pink
colormap([0 0 1]);
axis equal
hold on;
%plot3(S(1,:),S(2,:),S(3,:),'b');
plot3(S(1,:),S(2,:),S(3,:),'ks','LineWidth',2);
plot3(0,0,0,'ro','LineWidth',3);
xlabel('X');
ylabel('Y');
zlabel('Z');
%grid off;

%axis square

% % % k = 5;
% % % n = 2^k-1;
% % % theta = pi*(-n:2:n)/n;
% % % phi = (pi/2)*(-n:2:n)'/n;
% % %
% % % % % phi = [0:.3:pi]';
% % % % % step = (2*pi)/50;
% % % % % theta = 0:step:2*pi-step;
% % % % %
% % % X = a*cos(phi)*cos(theta);
% % % Y = b*cos(phi)*sin(theta);
% % % Z = c*sin(phi)*ones(size(theta));
% % % colormap([0 0 0;1 1 1])
% % % C = hadamard(2^k);
% % % surf(X,Y,Z,C)
% % % axis square