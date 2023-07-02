function ind = ChannelIndex1(x)
% Channel selection based on wavelet decomposition onto ECG bands, or
% equivalently using a filter that is customized to the fetal ECG frequency
% band
% 
% Reza Sameni
% December 2008

L1 = size(x,1);
% L2 = size(x,2);

load FiltCoefs Num

y = zeros(size(x));
ind = zeros(L1,1); 
for i = 1:L1
    y(i,:) = filter(Num,1,x(i,:));
    ind(i) = 100*var(y(i,:))/var(x(i,:));
end