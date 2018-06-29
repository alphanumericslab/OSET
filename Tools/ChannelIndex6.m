function ind = ChannelIndex6(x,ff,th)
% Channel selection based on peak-detection
%
% Reza Sameni
% December 2008

L1 = size(x,1);
L2 = size(x,2);


rng = floor(1/(2*ff));

% prefiltering
% load FiltCoefs Num
% x = filter(Num,1,x')';

% x = x.^2; % + and - peaks are both acceptable
x = abs(x); % + and - peaks are both acceptable

ind = zeros(L1,1);
for i = 1:L1
%     peaks = zeros(1,L2);
    MAX = max(x(i,:));
% % %     for j = 1:L2,
% % %         index = max(j-rng,1):min(j+rng,L2);
% % %         if(max(x(i,index))==x(i,j))
% % %             peaks(j) = 1;
% % %         end
% % %     end
% % % 
    peaks = PeakDetection(x(i,:),ff);
    I = find(peaks);
%     ind(i) = length(find(x(i,I)>=th*MAX));
    HR = diff(I);
%     ff;
    
% % %     % remove fake peaks
% % %     d = diff(I);
% % %     peaks(I(d<rng)) = 0;
% % % 
% % %     ind(i) = exp(-(ind(i) - (ff*L2))^2); % ff*L2 is the expected number of peaks. the closer to this value, the better!
% % %     ind(i) = length(find(peaks));
end



