function [score1, score2, peaks1, peaks2] = ChannelIndex12(x,f,fs)
% Channel selection based on fixed normalized template matched filter
% scatter and average waveform standard deviation
%
% Reza Sameni
% December 2015

% % % th = 0.5;
% % % h0 = load('BPFilter01.txt');

t1 = -.02:1/fs:.02; %ms
template1 = exp(-t1.^2/(0.007)^2);
% % % figure;plot(t1,template1);grid;

t2 = -.02:1/fs:.02+1/fs; %ms
template2 = diff(exp(-t2.^2/(0.007)^2));
% % % figure;plot(t2(1:end-1),template2);grid;

t3 = -.02:1/fs:.02; %ms
template3 = -.2*exp(-(t3+.007).^2/(0.005)^2) + exp(-t3.^2/(0.005)^2) + -.2*exp(-(t3-.007).^2/(0.005)^2);
% % % figure;plot(t3,template3);grid;

L = size(x,1);
N = size(x,2);

score1 = zeros(1,L);
score2 = zeros(1,L);
ppeaks = zeros(L, N);
for i = 1:L
    xx = x(i,:);
    
    L1 = length(template1);
    h1 = template1(end:-1:1);
    w1 = floor(L1/2);
    %     r1 = filter(h1,sqrt(sum(h1.^2)),[xx zeros(1,w1-1)])./sqrt(filter(ones(1, L1),1,[xx.^2 zeros(1,w1-1)]));
    r1 = filter(h1,sqrt(sum(h1.^2)*sum(xx.^2)),[xx zeros(1,w1-1)]);
    r1 = r1(w1:N+w1-1);
    
    L2 = length(template2);
    h2 = template2(end:-1:1);
    w2 = floor(L2/2);
    %     r2 = filter(h2,sqrt(sum(h2.^2)),[xx zeros(1,w2-1)])./sqrt(filter(ones(1, L2),1,[xx.^2 zeros(1,w2-1)]));
    r2 = filter(h2,sqrt(sum(h2.^2)*sum(xx.^2)),[xx zeros(1,w2-1)]);
    r2 = r2(w2:N+w2-1);
    
    L3 = length(template3);
    h3 = template3(end:-1:1);
    w3 = floor(L3/2);
    %     r3 = filter(h3,sqrt(sum(h3.^2)),[xx zeros(1,w3-1)])./sqrt(filter(ones(1, L3),1,[xx.^2 zeros(1,w3-1)]));
    r3 = filter(h3,sqrt(sum(h3.^2)*sum(xx.^2)),[xx zeros(1,w3-1)]);
    r3 = r3(w3:N+w3-1);
    
    r = sqrt(r1.^2 + r2.^2 + r3.^2);
    ppeaks(i, :) = PeakDetection(r, f/fs, 1);
    LL = max([L1 L2 L3]);
    peakpulses = conv(ppeaks(i, :), ones(1, LL), 'same');
    
    I_spikes = peakpulses == 1;
    I_nonspikes = peakpulses == 0;
    
    % first score based on the assumption of local peaks around the R-peak
    score1(i) = median(r(I_spikes))/median(r(I_nonspikes));
    
    % second score based on small beat variances
    I = find(ppeaks(i, :) == 1);
    I(1) = []; I(end) = []; % exclude the first and last beats from this score calculation to avoid boundary effects
    blocks = zeros(length(I), 2*LL);
    for k = 1 : length(I)
        blocks(k, :) = r(I(k) - LL : I(k) + LL - 1)/std(r); 
    end
    sample_std = std(blocks, [], 1);
    score2(i) = 1/mean(sample_std);
end

[~, ch] = max(score1);
peaks1 = ppeaks(ch, :);

[~, ch] = max(score2);
peaks2 = ppeaks(ch, :);

