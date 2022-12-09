function [s,peaksm,peaksf] = PiCA3(dat,refm,refframem,thm,wlen1m,wlen2m,reff,refframef,thf,wlen1f,wlen2f,Itr)


% % % t0 = 90;    % reference frame center
% % % w = 5;      % reference frame half width
% % % wlen1 = 10; % exclusion averaging frame half width (should be more than the fetal expected effective width)
% % % wlen2 = 50; % peak exclusion window half length (should be less than the minimum expected fetal RR-interval)
% % %
% % % ref = dat(2,:);
% % % refframe = dat(2,t0-w:t0+w);

N = size(dat,2);
for k = 1:Itr
    % first subspace
    w = floor(length(refframem)/2);
    r = xcorr(refframem,refm);
    r = r(:)';
    r = r(w+1:N+w);
    r = r(end:-1:1);
    r = r - Median(r,N,wlen1m)';
    r(r < thm*max(r)) = 0;
    peaksm = zeros(1,N);
    I = find(r>0);
    for i = I
        ind = max(1,i-wlen2m):min(N,i+wlen2m);
        if (max(r(ind))==r(i))
            peaksm(i) = 1;
        end
    end

    % second subspace
    w = floor(length(refframef)/2);
    r = xcorr(refframef,reff);
    r = r(:)';
    r = r(w+1:N+w);
    r = r(end:-1:1);
    r = r - Median(r,N,wlen1f)';
    r(r < thf*max(r)) = 0;
    peaksf = zeros(1,N);
    I = find(r>0);
    for i = I
        ind = max(1,i-wlen2f):min(N,i+wlen2f);
        if (max(r(ind))==r(i))
            peaksf(i) = 1;
        end
    end
    
    [s, ~, ~] = PiCA(dat,peaksm,peaksf);
    refm = s(1,:);
    reff = s(end,:);
end