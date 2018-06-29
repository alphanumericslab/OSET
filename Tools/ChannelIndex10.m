function [ind, peaks] = ChannelIndex10(x,f,fs)
% Channel selection based on fixed template matched filter
%
% Reza Sameni
% September 2010

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

score = zeros(6,L);
ppeaks = zeros(6,L,N);
for i = 1:L
    xx = x(i,:);

    L1 = length(template1);
    h1 = template1(end:-1:1);
    w1 = floor(L1/2);
    r1 = filter(h1,sqrt(sum(h1.^2)),[xx zeros(1,w1-1)]);
    r1 = r1(w1:N+w1-1);

    L2 = length(template2);
    h2 = template2(end:-1:1);
    w2 = floor(L2/2);
    r2 = filter(h2,sqrt(sum(h2.^2)),[xx zeros(1,w2-1)]);
    r2 = r2(w2:N+w2-1);

    L3 = length(template3);
    h3 = template3(end:-1:1);
    w3 = floor(L3/2);
    r3 = filter(h3,sqrt(sum(h3.^2)),[xx zeros(1,w3-1)]);
    r3 = r3(w3:N+w3-1);

    r = [r1 ; -r1 ; r2 ; -r2 ; r3 ; -r3];
    for kk = 1:size(r,1)
        % % %         ppeaks(kk,i,:) = PeakDetection6(r(kk,:),f/fs,.01,1); % search for positive peaks
        ppeaks(kk,i,:) = PeakDetectionFromC3(r(kk,:), f, fs, 1);
        I = find(ppeaks(kk,i,:));
        d = diff(I);
        score(kk,i) = std(d(3:end-2));
    end
    % % % % % % r(r<.1*std(r)) = 0;
    % % % % % % r(r>10*std(r)) = 0;

    % % %     r = [r1 ; r2 ; r3];
    % % %     r(r<.2*max(r)) = 0;
    % % %

    % % %     [tm k] = min(score(:,i));

    % % %     figure;
    % % %     hold on;
    % % %     plot(n,r);
    % % %     plot(n(I),r(I),'ro');
    % % %     grid
end

[ind, I] = min(score,[],1);
ind = ind';
[val, ch] = min(ind);
peaks = ppeaks(I(ch),ch,:);

% % % peaks = squeeze(peaks); %%% ADDED TO WORK FOR SINGLE CHANNEL DATA
% % % peaks = peaks(:)';
