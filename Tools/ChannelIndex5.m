function ind = ChannelIndex5(x)
% Channel selection based on zero-crossings
% 
% Reza Sameni
% December 2008

L1 = size(x,1);
L2 = size(x,2);

ind = zeros(L1,1); 
for i = 1:L1
    sgn = x(i,1:end-1) .* x(i,2:end);
    I = find(sgn<=0);
    ind(i) = 100*length(I)/(L2-1);
end