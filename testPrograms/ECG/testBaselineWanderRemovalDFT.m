% test baseline wander removal by freq. domain filtering
% Reza Sameni, 2018
%

clear;
close all;
load patient165_s0322lre;
%load LOUETTE_07-Oct-2005_40s.txt -mat;data = LOUETTE_40s;clear LOUETTE_40s

fs = 1000;
%d = data(:,4:8);
% d = data(:,2:end);
d = val';
N = length(d);

S = fft(d,N,1);
SS = 10*log10(S.^2);
f = [0:N-1]*fs/N;
plot(f,SS);grid;

S(1,:) = 0;

k1 = 2:ceil(.4*N/fs);
S(k1,:) = 0; S(N-k1+2,:) = 0;

k2 = floor(49.5*N/fs):ceil(50.5*N/fs);
S(k2,:) = 0; S(N-k2+2,:) = 0;

k3 = floor(149.5*N/fs):ceil(150.5*N/fs);
S(k3,:) = 0; S(N-k3+2,:) = 0;

k4 = floor(249.5*N/fs):ceil(250.5*N/fs);
S(k4,:) = 0; S(N-k4+2,:) = 0;

k5 = floor(349.5*N/fs):ceil(350.5*N/fs);
S(k5,:) = 0; S(N-k5+2,:) = 0;

SS1 = 10*log10(S.^2);
dd = real(ifft(S,N,1));
rr = imag(ifft(S,N,1));
figure;
plot(f,SS1);grid;

figure;
plot(d);grid;
figure;
plot(dd);grid;

% L = 3;
% for i = 1:size(data,2),
%     figure(floor((i-1)/L)+1);
%     subplot(L,1,mod(i-1,L)+1);
%     plot(data(:,i));grid;
%     ylabel(num2str(i));
% end
% 
