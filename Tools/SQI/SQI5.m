function ind = SQI5(x, winLen)
% Channel selection based on channel power in sliding window
%
% Fahimeh Jamshidian Tehrani
% February 2020

% Comment out for individual using
% % %     Remove the mean
% % mn = mean(x,2)*ones(1,size(x,2));
% % x = x - mn;
% % 
% % %     Normalize the variance of x
% % x = x./var(x,0,2);

[CH, ~] = size(x);

for ch = 1: CH
    yk(ch,:) = movmean(x(ch, :).^2, winLen);
    ind(ch) = var(yk(ch,:));
end