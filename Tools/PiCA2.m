function [s,peaks] = PiCA2(dat,ref,refframe,th,wlen1,wlen2,Itr)


% % % t0 = 90;    % reference frame center
% % % w = 5;      % reference frame half width
% % % wlen1 = 10; % exclusion averaging frame half width (should be more than the fetal expected effective width)
% % % wlen2 = 50; % peak exclusion window half length (should be less than the minimum expected fetal RR-interval)
% % %
% % % ref = dat(2,:);
% % % refframe = dat(2,t0-w:t0+w);

N = length(ref);
w = floor(length(refframe)/2);
for k = 1:Itr,
    % % %     if k == 1,
    % % %         peaks = zeros(1,N);
    % % %         peaks(1:250/2:end) = 1;
    % % %     else
    
    r = xcorr(refframe,ref);
    r = r(:)';
    r = r(w+1:N+w);
    r = r(end:-1:1);
    r = r - Median(r,N,wlen1)';
    r(r < th*max(r)) = 0;
    peaks = zeros(1,N);
    I = find(r>0);
    for i = I,
        ind = max(1,i-wlen2):min(N,i+wlen2);
        if (max(r(ind))==r(i)),
            peaks(i) = 1;
        end
    end
    
    % % %     end

    [s,W,A] = PiCA(dat,peaks);
    ref = s(1,:);
end