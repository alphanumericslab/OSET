% baseline wander removal based on sparsity assumption of the BW in the
% frequency domain;
% Reza Sameni, Copyright 2006
%

clear;
close all;
clc;
%load patient165_s0323lre; data = data(1:10000,6)';
load FOETAL_ECG.dat; data = FOETAL_ECG(1:end,5)';clear FOETAL_ECG;
%plot(data);grid;

fs = 500;
N = length(data);
%w = 2*pi*[.1 , .2 , .3 , fs-.1 , fs-.2 , fs-.3]/fs;
%w = 2*pi*[2 fs-2 3 fs-3 .7 fs-.7]/fs;

%K = [1:ceil(N*.5/fs) 50*N/fs];
K = [1:ceil(N*.5/fs)];
%w = 2*pi*[K N-K]'/N;
ff = [];
w = 2*pi*[ [K N-K]/N [ff fs-ff]/fs ]';

% f = .1:.1:.4;
% w = 2*pi*[f fs-f]/fs;

L = length(w);

noise = 0*(.7*sin(2*pi*[0:N-1]*.11/fs) + .3*cos(2*pi*[0:N-1]*.31/fs + pi/3) + .2*cos(2*pi*[0:N-1]*.21/fs));
data = data + noise;

T = zeros(L);
U = zeros(L,1);
j = sqrt(-1);

% % % for p = 1:L,
% % %     for k = 1:L,
% % %         coef = 0;
% % %         for m = 1:N,
% % %             coef = coef + sum((exp(j*w(p)*m)-exp(j*w(p)*[1:N])).*(exp(j*w(k)*m)-exp(j*w(k)*[1:N])));
% % %         end
% % %         T(p,k) = coef;
% % %     end
% % %
% % %     coef = 0;
% % %     for m = 1:N,
% % %         coef = coef + sum((exp(j*w(p)*m)-exp(j*w(p)*[1:N])).*(data(m)-data([1:N])));
% % %     end
% % %     C(p) = coef;
% % % end

% % % for p = 1:L,
% % %     numr = 0;
% % %     numi = 0;
% % %     den = 0;
% % %     for m = 1:N,
% % %         numr = numr + sum( (data(m)-data([1:N])).*(cos(w(p)*m)-cos(w(p)*[0:N-1])) );
% % %         numi = numi + sum( (data(m)-data([1:N])).*(sin(w(p)*m)-sin(w(p)*[0:N-1])) );
% % %         den = den + 2*sum( 1 - cos(w(p)*(m-[0:N-1])) );
% % %     end
% % %     A(p) = (numr + j*numi)/den;
% % % end
% A = inv(T)*C;

% % % cntr = 0;
% % % for m = 1:N,
% % %     %for n = m+1:N,
% % %     %for n = m+1:1170:N,
% % %     % % %     for n = m+1:183:N, %maternal
% % %     % % %         %for i = max(1,n-5):min(N,n+5),
% % %     % % %             C = cos(w*m)-cos(w*n);
% % %     % % %             S = sin(w*m)-sin(w*n);
% % %     % % %             %T = T + (C*C'+S*S');
% % %     % % %             U = U + (C-j*S)*(data(m)-data(n));
% % %     % % %             cntr = cntr + 1;
% % %     % % %         %end
% % %     % % %     end
% % %     %for n = m+1:115:N, %fetal
% % %     for n = max(1,m-50):min(N,m+50),
% % %         %for i = max(1,n-5):min(N,n+5),
% % %         C = cos(w*m)-cos(w*n);
% % %         S = sin(w*m)-sin(w*n);
% % %         T = T + (C*C'+S*S');
% % %         U = U + (C-j*S)*(data(m)-data(n));
% % %         cntr = cntr + 1;
% % %         %end
% % %     end
% % % 
% % % end

m = [0:N-1];
C = cos(w*m);
S = sin(w*m);

E = 2*N*C*C'-2*sum(C,2)*sum(C,2)'+2*N*S*S'-2*sum(S,2)*sum(S,2)';
F = 2*N*S*C'-2*sum(S,2)*sum(C,2)'-2*N*C*S'+2*sum(C,2)*sum(S,2)';
a = 2*N*C*data' - 2*sum(C,2)*sum(data);
b = 2*N*S*data' - 2*sum(S,2)*sum(data);

alpha = pinv(pinv(E)*F + pinv(F)*E) * (pinv(F)*a + pinv(E)*b);
beta = pinv(pinv(E)*F + pinv(F)*E) * (pinv(E)*a - pinv(F)*b);
A = alpha + j*beta;


%A = inv(T)*U;
%A = U/N^2;
%A = U/cntr/2;

baseline = 0;
for k = 1:L,
    baseline = baseline + A(k)*exp(j*w(k)*[0:N-1]);
end

DATA = fft(data,N);
B = zeros(N,1);
I = [K N-K]+1;
B(I) = DATA(I);
baseline2 = ifft(B,N);
% figure;
% plot(data);
% grid;
% hold on;
% plot(abs(baseline),'r');
% plot(real(baseline),'g');
% plot(imag(baseline),'m');

figure;
plot(data);
grid
hold on;
plot(data-real(baseline),'r');
plot(real(baseline),'m--');
plot(real(baseline2),'g--');



