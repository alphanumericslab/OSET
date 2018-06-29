function [ind, xcore] = ChannelIndex11(x,fs,f0)
% [ind xcore] = ChannelIndex11(x,fs,f0)
% Channel saturation and sine wave modulation detector
% Inputs:
%   x : input data matrix
%   fs : input sampling frequency
%   f0 : desired sine (or square) wave frequency (approximate value is ok, e.g., f0 = 4)
% Outputs:
%   ind : the evaluation index in seconds. The smaller the value, the
%   higher the prossibility of having a sine wave or saturation. By setting
%   a threshold on this value, one can detect the saturated channels.
%   According to my simulations, 0.05s works for most simulated cases.
%
%   xcore: the correlation matrix calculated in this function (for displaying and research purposes only)
%
% Algorithm Logic: The logic behind the technique is that channels
% saturated by periodic noise have regular peaks in the auto-correlation
% function. The more regular the peaks are (smaller values of ind), the
% more probable it is that the channel is corrupted with periodic noise.
%
% Copyright (C) Reza Sameni
% January 2013
% Shiraz University, Shiraz, Iran
% reza.sameni@gmail.com

L1 = size(x,1);
L2 = size(x,2);
ind = zeros(L1,1);
xcore = zeros(L1,2*L2-1);

for i = 1:L1
    xcore(i,:) = xcorr(x(i,:))./[1:L2 L2-1:-1:1]; % calculate the auto-correlation of each channel normalized by the number of overlaping samples
    px = PeakDetection(xcore(i,:),f0/fs,1); % find positive peaks
    ix = find(px(round(L2/2):round(3*L2/2)));
    ind(i) = std(diff(ix))/fs;
    %     figure
    %     hold on
    %     plot(xcore(i,:));
    %     ii = find(px==1);
    %     plot(ii,xcore(i,ii),'ro');
    %     grid
end