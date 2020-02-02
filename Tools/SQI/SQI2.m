function ind = SQI2(x, peaks)
% Channel selection based on R-peak amplitude constancy
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
    Peakloc{ch} = find(peaks{ch});
    ind(ch) = var(x(ch,Peakloc{ch}));
end