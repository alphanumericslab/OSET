% test fECG extraction using an additional synthetic channel
% Reza Sameni, Copyright March 2006
% Modified June 2018
%

clear;
close all;
clc;
%load LOUETTE_07-Oct-2005_20s.txt -mat;data = LOUETTE_20s';clear LOUETTE_20s;
%load patient165_s0323lre;data = data(:,2:end)';
load FOETAL_ECG.dat; data = FOETAL_ECG(:,2:end)';clear FOETAL_ECG; fs = 250;

%//////////////////////////////////////////////////////////////////////////
% ICA 0
W0 =  jadeR(data);
s0 = W0*data;
A0 = pinv(W0);
%//////////////////////////////////////////////////////////////////////////
N = length(data);
% fetal peak detection
%th = [4.8,-3.8,5,3.2,0,0,3,2.683]; % from s0
th = 3;
ch = s0(7,:);
phase = zeros(1,length(data));
cntr = 0;
for i = 1:N
    %if(ch(i)>th & max(ch(max(i-10,1):min(i+10,N)))==ch(i)),
    if( (th>0 && ch(i)>th && max(ch(i-10:i+10))==ch(i)) || (th<0 && ch(i)<th && min(ch(i-10:i+10))==ch(i)) ),
        inds = i-0:i+0;
        %phase(inds) = 1000*(10-abs(inds-i));%100*exp(-(inds-i).^2/100);
        phase(inds) = -1;
        cntr = cntr + 1;
    end
end
plot(ch);
hold on;
plot(phase,'r');

%//////////////////////////////////////////////////////////////////////////
% ICA 1
data2 = [data; phase]
W1 =  jadeR(data2);
data2(end,:) = 0;
s3 = W1*data2;

%//////////////////////////////////////////////////////////////////////////
L = 1;
t = (0:N-1)/fs;
close all;
for i = 1:9,
    figure(floor((i-1)/L)+1);
    subplot(L,1,mod(i-1,L)+1);
    hold on;
    for j = 1:size(s0,1),
        %mtx = corrcoef(s3(i,:),s0(j,:));
        %c(j) = mtx(1,2);
        c(j) = max(xcorr(s3(i,:)/std(s3(i,:)),s0(j,:)/std(s0(j,:))));
    end
    [Y,J] = max(c);
    %J = find(max(c));
    plot(t, s0(J,:)/std(s0(J,:)),'r');
    plot(t, s3(i,:)/std(s3(i,:)));
    xlabel('time(s)');
    grid;
end
