% test VCG loop alignment
% UNDER TEST
% Reza Sameni, 2005

clear;
close all;
load patient165_s0322lre;data = val(13:15, 1:15500)';
%load FOETAL_ECG.dat; data = FOETAL_ECG(:,2:9);clear FOETAL_ECG;

N = 70;
%N = 200;
delta = 5;
%//////////////////////////////////////////////////////////////////////////
% filter data
fs = 1000;
L = length(data);
S = fft(data,L,1);
S(1,:) = 0;
k1 = 2:ceil(0.5*L/fs);
S(k1,:) = 0; S(L-k1+2,:) = 0;
for i = 1:9,
    k = floor((i*50-.3)*L/fs):ceil((i*50+.3)*L/fs);
    S(k,:) = 0; S(L-k+2,:) = 0;
end
data = real(ifft(S,L,1));
clear S*; clear k*
%//////////////////////////////////////////////////////////////////////////
% peak detection
%plot(data(:,1));grid
phase = zeros(1,length(data));
cntr = 0;
for i = 2:length(data)-1
    if(data(i,1)>1 & max(data(i-1:i+1,1))==data(i,1))
    %if(data(i,2)>400 & max(data(i-1:i+1,2))==data(i,2))
    %if(data(i,7)>400 & max(data(i-1:i+1,7))==data(i,7))
        phase(i) = 1;
        cntr = cntr + 1;
    end
end
%cntr
%//////////////////////////////////////////////////////////////////////////
% extracting the ECGs
I = find(phase==1);
AllECG = zeros(N+2*delta,size(data,2),cntr);
for i = 2:cntr,
    AllECG(:,:,i-1) = data(I(i)-N/2-delta:I(i)+N/2+delta-1,:);
end
% figure
% plot3(squeeze(AllECG(:,1,3:end-3)),squeeze(AllECG(:,2,3:end-3)),squeeze(AllECG(:,3,3:end-3)));grid
figure
plot3(squeeze(AllECG(:,1,:)),squeeze(AllECG(:,2,:)),squeeze(AllECG(:,3,:)));grid;title('before alignment');
cntr = cntr-1;
%//////////////////////////////////////////////////////////////////////////
%
%ZR = squeeze(AllECG(:,:,30))';
%ZR = squeeze(AllECG(:,:,2))';
ZR = squeeze(mean(AllECG(:,:,3:end-3),3))';
er = zeros(2*delta,1);
for i = 1:cntr,
    Z = AllECG(delta:end-delta-1,:,i)';
    for tau = -delta+1:delta-1,
        J = [zeros(delta+tau-1,N) ; eye(N) ; zeros(delta-tau+1,N)];
        C = Z*J'*ZR';
        [U,S,V] = svd(C);
        Q = U*V';
        alpha = trace(Z'*Q*ZR*J)/trace(J'*ZR'*ZR*J);
        X = Z - alpha*Q*ZR*J;
        er(tau+delta) = sqrt(sum(diag(X'*X)));
    end
    [mn,I] = min(er);
    tauhat(i) = I-delta;
    J = [zeros(delta+tauhat(i)-1,N) ; eye(N) ; zeros(delta-tauhat(i)+1,N)];
    C = Z*J'*ZR';
    [U,S,V] = svd(C);
    Qhat(:,:,i) = U*V';
    Alphahat(i) = trace(Z'*Qhat(:,:,i)*ZR*J)/trace(J'*ZR'*ZR*J);
end
%//////////////////////////////////////////////////////////////////////////
for i = 1:cntr,
    Z(:,:,i) = Alphahat(i)*Qhat(:,:,i)*ZR*[zeros(delta+tauhat(i)-1,N) ; eye(N) ; zeros(delta-tauhat(i)+1,N)];
end

figure
plot3(squeeze(Z(1,:,:)),squeeze(Z(2,:,:)),squeeze(Z(3,:,:)));grid;title('after alignment');
% figure
% plot3(squeeze(Z(1,:,:)),squeeze(Z(2,:,:)),squeeze(Z(3,:,:)));grid

%plot(squeeze(Z(3,:,3:end-3)));grid

figure
for i = 1:3,
    subplot(1,3,i)
    plot(std(Z(i,:,3:end-3),[],3))
    hold on;
    plot(std(AllECG(:,i,3:end-3),[],3),'r')
    grid
end

