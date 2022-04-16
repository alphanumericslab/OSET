% fastICA test
% Reza Sameni, 2005

close all;


% teta = -4*pi:.01:4*pi;
% r = 1*sin(5*teta).^9 + 1*sin(teta);
% [s1,s2] = pol2cart(teta,r);
% N = length(teta);
Wp = .2;
Ws = .25;
Rp = .5;
Rs = 60;

[N, Wn] = ellipord(Wp, Ws, Rp, Rs);
[B,A] = ellip(N,Rp,Rs,Wn);

N = 1000;
n = 0:N-1;
s1 = sin(2*pi*.003*n).^2;
s2 = 1*ones(size(n));a = mod(floor(n/53),2); I = find(a==0); s2(I) = 0;
% s3 = 2*mod(n,N/7)/(N/7)-1;
% h = 2*mod(n,N/7)/(N/7)-1;
s3 = .3*(rand(size(n))-.5);

% s1 = filter(B,A,s1);
% s2 = filter(B,A,s2);
s3 = filter(B,A,s3);

% [TH,PHI,R] = CART2SPH(X,Y,Z)

% S = [s1 ; s2 ; s3 ; s4];
S = [s1 ; s2 ; s3];

% A = [1 2 3 1; 1 -1 4 .1; 3 2 5 -2 ; 4 -3 2 -1];
A = [1 2 3; 1 -1 4; 3 2 5];

X = A*S;

% X = [X -X(:,end:-1:1)];

% [Shat,Ahat,What] = fastica(X, 'displayMode', 'off','approach', 'symm','g','tanh','epsilon',1e-6);
[Shat,Ahat,What] = fastica(X, 'displayMode', 'off','approach', 'symm');%,'epsilon',1e-9,'maxFinetune',1000);
% Shat = (Shat(:,1:end/2) - Shat(:,end:-1:end/2+1))/2;

figure
plot(S');grid;
figure
plot(X');grid;
figure
plot(Shat');grid;
figure
plot3(X(1,:),X(2,:),X(3,:),'.');grid;
