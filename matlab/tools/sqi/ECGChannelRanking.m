function [ind, rank] = ECGChannelRanking(x, fs)
% Ranking channels based on multiple indexes
% Reza Sameni, Copyright
% December 2008
%
L1 = size(x,1);
L2 = size(x,2);
t = (0:L2-1)/fs;

%//////////////////////////////////////////////////////////////////////////
% preprocessing

% mean removal (not really necessary!)
% % % x = x - mean(x,2)*ones(1,L2);
% % % 
% % % % baseline wander removal
% % % for i = 1:L1,
% % %     x(i,:) = x(i,:) - nonlinear_filter(x(i,:),'trmean',round(.300*fs),round(.10*fs),round(.10*fs));
% % % end
% % % 
% % % % out-of-band noise removal
% % % x = x - LPFilter(x,.7/fs);
% % % x = LPFilter(x,200/fs);
% % % 
% % % % notch filtering
% % % Wo = 60/(fs/2);  BW = Wo/35;
% % % [b,a] = iirnotch(Wo,BW);
% % % x = filter(b,a,x')';
% % % 
% % % % channel centralization & normalization
% % % for i = 1:L1,
% % %     x(i,:) = (x(i,:)-mean(x(i,:)))/std(x(i,:));
% % % end

%//////////////////////////////////////////////////////////////////////////
% channel evaluation routines

% reject channels that are close to saturation (the evaluation is done on the raw data before preprocessing)
ind00 = ChannelIndex0(x,.95*max(x,[],2),max(x,[],2)); % test for upwards saturation 
ind01 = ChannelIndex0(x,min(x,[],2),.95*min(x,[],2)); % test for downwards saturation
ind0 = ind00 + ind01;

% finding the channels that are within the frequency band of the fetal ECG
% (the internal filter of this function has been designed based on wavelet decomposition using Coiflet mother wavelets)
ind1  = ChannelIndex1(x);

% finding the channels that have energy deviations within short windows (test of nonstationarity over short windows)
% (since a good ECG channel should have high energy deviations in windows with a fractional length of the average heart beat)
ind2 = ChannelIndex2(x,.2*fs);

% finding the channels that have energy deviations within long windows (test of stationarity over long windows)
% (since a good ECG channel should have low energy deviations in long windows)
ind3 = 1./ChannelIndex2(x,3*fs);

% finding nongaussian channels using negentropy estimation of gaussianity
ind4 = ChannelIndex4(x);

% channel evaluation based on zero-crossings; a clean signal should not have too many zero-crossings
ind5 = 1./ChannelIndex5(x);

% channel evaluation based on number of peaks within windows with a length of the average fetal heart rate
ind6 = ChannelIndex6(x,2.0/fs,.1);

% channel evaluation based on number of level crossings of signal energy
% % % [ind7 yy] = ChannelIndex7(x,.2*fs,.2);

%//////////////////////////////////////////////////////////////////////////
% channel selection based on 'weighted voting' among the different extracted measures

rank0 = ones(L1,1);
rank0(ind0>30) = L1;   % hard thresholding

[temp, I] = sort(ind1,'descend');
[temp rank1] = sort(I,'ascend');

[temp, I] = sort(ind2,'descend');
[temp rank2] = sort(I,'ascend');

[temp, I] = sort(ind3,'descend');
[temp rank3] = sort(I,'ascend');

[temp, I] = sort(ind4,'descend');
[temp rank4] = sort(I,'ascend');

[temp, I] = sort(ind5,'descend');
[temp rank5] = sort(I,'ascend');

[temp, I] = sort(ind6,'descend');
[temp rank6] = sort(I,'ascend');

% [temp, I] = sort(ind7,'descend');
% [temp rank7] = sort(I,'ascend');

% the overall index (without voting)
% ind = (ind0.*ind1.*ind2.*ind3.*ind4.*ind5.*ind6).^(1/7);
% ind = ind7; % single measure analysis (for testing)

% weighted voting between the ranks of each measure (the lower the rank, the better!)
% we can play with the weights for fine-tuning the performance

% rank = ((rank0-1) + 2*(rank1-1) + (rank2-1) + 0.5*(rank3-1) + 0.5*(rank4-1) + 0.5*(rank5-1) + 2*(rank6-1))/7.5;
rank = rank6-1; % single measure analysis (for testing)
ind = ind6;
[temp,I] = sort(rank,'ascend');

% % % % % %//////////////////////////////////////////////////////////////////////////
% % % % % % plotting the ranked channels
% % % % % L = 4;
% % % % % t = (0:L2-1)/fs;
% % % % % for i = 1:L1,
% % % % %     if(mod(i,L)==1 || L==1)
% % % % %         figure;
% % % % %     end
% % % % %     subplot(L,1,mod(i-1,L)+1);
% % % % %     plot(t,x(I(i),:));
% % % % %     hold on;
% % % % % %     plot(t,yy(I(i),:),'r');
% % % % %     ylabel([num2str(rank(I(i)),'%2.2f') ' (' num2str(ind(I(i)),'%2.2f') ')' ' [' num2str(I(i),'%d') ']']);
% % % % %     grid;
% % % % %     if(mod(i,L)==0 || L==1)
% % % % %         xlabel('time(s)');
% % % % %     end
% % % % % end
% % % % % 
% % % % % toc