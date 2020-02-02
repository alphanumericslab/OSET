function ind = SQI6(x, peaks)
% Channel selection based on measure of PICA periodicity
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
    Peakloc = find(peaks{ch});
    T = mean(diff(Peakloc));
    [T0,T1] = SynchPhaseTimes2(peaks{ch});
    be4peak = 1: T0(1)-1;
    T0 = [be4peak T0];
    T1 = [be4peak+round(T) T1];
    
    term1 = mean(x(ch, T0).*x(ch, T1));
    term2 = mean(x(ch,:).^2);
    ind(ch) = term1 / term2;
end