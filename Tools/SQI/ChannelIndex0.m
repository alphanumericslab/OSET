function ind = ChannelIndex0(x,th1,th2)
% Channel selection based on signal range (thresholding)
% 
% Reza Sameni
% December 2008

L1 = size(x,1);
L2 = size(x,2);

ind = zeros(L1,1); 
for i = 1:L1
    I = find(x(i,:) >= th1(i) & x(i,:) <= th2(i));
    ind(i) = 100*length(I)/L2;
end
