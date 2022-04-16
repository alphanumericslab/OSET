function ind = ChannelIndex2(x,w)
% Channel selection based on energy deviations
% 
% Reza Sameni
% December 2008

L1 = size(x,1);
% L2 = size(x,2);

y = zeros(size(x));
ind = zeros(L1,1);
for i = 1:L1
    y(i,:) = sqrt(filter(ones(1,round(w)),round(w),x(i,:).^2));
    ind(i) = 100*std(y(i,:));
    
% % %     figure;
% % %     plot(x(i,:));
% % %     hold on;
% % %     grid;
% % %     plot(y(i,:),'r');
end