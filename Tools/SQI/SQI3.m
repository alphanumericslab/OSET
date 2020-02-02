function ind = SQI3(x)
% Channel selection based on measure of non-gaussianity
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
k1 = 36/(8*sqrt(3)-9);
k2 = 24/(16*sqrt(3)-27);

for ch = 1: CH
   term1 = k1 * mean(x(ch,:).*exp(-(x(ch,:).^2)/2))^2; 
   term2 = k2 * (mean(exp(-(x(ch,:).^2)/2))-sqrt(2)/2)^2;
   ind(ch) = abs(term1 + term2);
end