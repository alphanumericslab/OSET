function ind = ChannelIndex4(x)
% Channel selection based on negentropy estimation
%
% Reza Sameni
% December 2008

L1 = size(x,1);
% L2 = size(x,2);

% y = zeros(size(x));
ind = zeros(L1,1);

for i = 1:L1
    % channel centralization & normalization
    %     x(i,:) = (x(i,:)-mean(x(i,:)))/std(x(i,:)); % already done before this module....

    k1 = 36/(8*sqrt(3)-9);
    k2 = 24/(16*sqrt(3)-27);
    ind(i) = k1*mean(x(i,:).*exp(-x(i,:).^2/2))^2 + k2*(mean(exp(-x(i,:).^2/2))-sqrt(.5))^2;

    % alternative method:
    % % %     k1 = 36/(8*sqrt(3)-9);
    % % %     k2 = 1/(2-6/pi);
    % % %     ind(i) = k1*mean(x(i,:).*exp(-x(i,:).^2/2))^2 + k2*(mean(abs(x(i,:)))-sqrt(2/pi))^2;
end


