function [ind, yy] = ChannelIndex7(x,w,th)
% Channel selection based on number of level crossings of signal energy
% 
% Reza Sameni
% December 2008

L1 = size(x,1);
% L2 = size(x,2);

yy = zeros(size(x));
ind = zeros(L1,1);
for i = 1:L1
    s = sqrt(filter(ones(1,round(w)),round(w),x(i,:).^2));
    y = s - max(s)*th;
    yy(i,:) = y;
    u = y(1:end-1).*y(2:end);
    I = find(u<=0);
%     ind(i) = 100*length(I)/L2;
    ind(i) = length(I);
end