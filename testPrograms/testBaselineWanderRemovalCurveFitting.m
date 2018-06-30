% test baseline wander removal by curve fitting
% Reza Sameni, 2005

clear;
close all;
% load LOUETTE_07-Oct-2005_40s.txt -mat;data = LOUETTE_40s';clear LOUETTE_40s
% load fortict_20-Jan-2006.datfortict_20-Jan-2006.dat_acq.mat;
% load fortict5_20-Jan-2006.datfortict5_20-Jan-2006.dat_acq.mat;data = data';
%load Bonhomme3_20-Jan-2006.datBonhomme3_20-Jan-2006.dat_acq.mat;
load patient165_s0323lre; data = data(1:40000,2:6)';
% load FOETAL_ECG.dat; data = FOETAL_ECG(:,2:end)';clear FOETAL_ECG;

%//////////////////////////////////////////////////////////////////////////
% ICA 2
% W0 =  jadeR(data)
% s0 = W0*data;
% A0 = pinv(W0);
%//////////////////////////////////////////////////////////////////////////
% First stage of baseline wander removal
%th = [4.8,3.8,5,3.2,0,0,3.1,2.683]; % from s0
%th = [-30,80,-60,-24,-80,-600,600,600]; % from data of FOETAL_ECG
th = [.465 1.5 1 -.8 -.3]; % from patient165_s0323lre
datalen = length(data);

mn = mean(data(:,400:600),2);

AA = ones(size(data));
BB = zeros(size(data));
data2 = zeros(size(data));
for ch = 1:5;
    Channel = data(ch,:);

    % baseline removal 2
    peaks = zeros(datalen,1);
    cntr = 0;
    for i = 1:datalen
        index = max(i-10,1):min(i+10,datalen);
        if( (th(ch)>0 && Channel(i)>th(ch) && max(Channel(index))==Channel(i)) || (th(ch)<0 && Channel(i)<th(ch) && min(Channel(index))==Channel(i)) )
            peaks(i) = 1;
            cntr = cntr + 1;
        end
    end
    I = find(peaks==1);
    % remove possible double close peaks
    df = diff(I);
    df_true = find(df>20);
    I = I(df_true);

    % QRS
    shift=-3:3;
    Alpha = ones(size(I));
    Beta = zeros(size(I));
    winlen = 40/2;

    index0 = I(1) - winlen : I(1) + winlen;
    xref = Channel(index0);
    for i = 1:length(I),
        for j = 1:length(shift);
            index = I(i) - winlen + shift(j) : I(i) + winlen + shift(j);
            alpha(j) = (mean(xref.*Channel(index)) - mean(xref)*mean(Channel(index)))/(mean(Channel(index).^2)-mean(Channel(index))^2);
            beta(j) = mean(xref - alpha(j)*Channel(index));
            er(j) = sum((xref-alpha(j)*Channel(index)-beta(j)*ones(1,length(index))).^2);
        end
        [y,k] = min(er);
        Alpha(i) = alpha(k);
        Beta(i) = beta(k);
        %xref = (xref*(i-1) + Alpha(i)*Channel(index) + Beta(i))/i;

        I(i) = I(i) + shift(k);
        peaks(I) = 1;
    end

    %if(ch==8)
    %   PeakIndex = I;
    %end
    AA(ch,:) = interp1([1 ; I ; datalen],[1 ; Alpha ; 1],1:datalen,'cubic');
    BB(ch,:) = interp1([1 ; I ; datalen],[0 ; Beta ; 0],1:datalen,'cubic');
    data2(ch,:) = AA(ch,:).*Channel + BB(ch,:) - mn(ch);
    figure
    plot(Channel);
    grid;
    hold on;
    plot(data2(ch,:),'r','LineWidth',1);

    %plot(-BB,'m','LineWidth',1);
    %plot(1./AA,'g','LineWidth',1);
end

r0 = sqrt(sum(data.^2,1));
r = sqrt(sum(data2.^2,1));
figure
plot(r0);
grid;
hold on;
plot(r,'r');
