clear;
close all;

N = 30000;
% fs = 393.4;
order = 12;
q = (.5e-3)^2;
R = (50)^2;%
p0 = (2e-0)^2;%
alpha = 1;

beta = .99;
gamma = 1;

load WholeBW;
d = data(3,1:N);
time = data(1,1:N)
% N = length(d);

d = 50*(d-mean(d))/std(d);
d = d - 1.1*min(d);
vmean = mean(d);

% d0 = 1*randn(N,1);
% AA = [1 -1 .01 .002];
% d = filter(1,AA,d0)';
%d = 5*(d-mean(d))/std(d);
% figure
% freqz(1,AA);

[a0,e] = aryule(d,order);
[Ahat,Asmoothed,P] = TimeVariantAR(d',order,a0(2:end)',q,R,p0,alpha,beta,gamma,vmean);

x = 1*randn(N,1)+mean(d);
% x = (x-mean(x))/std(x);
% x = d0;

y =  zeros(N,1);
y2 =  zeros(N,1);
sm =  zeros(N,1);

for i = order+1:N-1,
    y(i) = (sqrt(1)*x(i)-Ahat(:,i)'*y(i-1:-1:i-order))/1;
    y2(i) = (sqrt(1)*x(i)-Asmoothed(:,i)'*y2(i-1:-1:i-order))/1;
    sm(i) = 1/(1+sum(Ahat(:,i)));
end

% figure;
% plot(sqrt(P));
% grid;

figure;
plot(time,d);
%grid;
%figure;
hold on;
plot(time,y,'r');
plot(time,y2,'m--');
grid;

figure;
plot(time,d/100);
grid;
xlabel('time (sec.)');
ylabel('mV');

figure;
plot(time,y/100,'r');
grid;
xlabel('time (sec.)');
ylabel('mV');


% figure
% plot(sm,'m');
% grid;
