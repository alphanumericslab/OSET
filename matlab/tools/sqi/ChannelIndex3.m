function ind = ChannelIndex3(x,w,bin,N)
% Channel selection based on energy deviations followed by a bandpass
% filter using a DFT filter bank
% 
% Reza Sameni
% December 2008

% L1 = size(x,1);
L2 = size(x,2);

% y = zeros(size(x));
% % % ind = zeros(L1,1);

x = (x - mean(x,2)*ones(1,L2)) ./ (std(x,[],2)*ones(1,L2));
u = filter(ones(1,round(w)),round(w),(x.^2)')';
u = sqrt(u);
% % % y = BPFilterComplex(u,bin/N,2*N);

ind = std(u,[],2);

% ind = 100*sum(abs(y).^2,2)./sum(abs(u).^2,2);
% ind = 100*sum(abs(y(:,L2/2:L2)).^2,2)./sum(abs(u(:,L2/2:L2)).^2,2);

% % % % % for i = 1:L1,
% % % % % % % %     x(i,:) = x(i,:)/std(x(i,:));
% % % % % % % %     y(i,:) = sqrt(filter(ones(1,round(w)),round(w),x(i,:).^2));
% % % % % % % %     Y = abs(fft(y(i,:),2*N));
% % % % %    
% % % % %     figure;
% % % % % % % %     subplot(211);
% % % % %     plot(x(i,:));
% % % % %     hold on;
% % % % %     plot(real(u(i,:)),'r');
% % % % % % % %     plot(real(y(i,:)),'m');
% % % % %     grid;
% % % % % % % %     subplot(212);
% % % % % % % %     plot(20*log10(abs(Y)));
% % % % % % % %     ind(i) = 100*Y(bin)^2/sum(Y(1:N).^2);
% % % % %     grid;
% % % % %     ylabel(num2str(ind(i)));
% % % % % end