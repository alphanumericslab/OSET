function [s,peaks] = PiCA5(dat,refslice,th,wlen,Itr)

N = size(dat,2);
w = floor(size(refslice,2)/2);
ref = dat(1,:);

% % % peaks = zeros(1,N);
% % % peaks(160:round(1025/2.4):end) = 1;
% % % 
% % % [s,W,A] = PiCA(dat,peaks);
% % % ref = s(1,:);
% % % refslice = W*refslice;

for k = 1:Itr
    refframe = refslice(1,:);
    r = zeros(1,N);
    for i = 1:N
        ind = max(1,i-w):min(N,i+w);
        cr = corrcoef(refframe(1:length(ind)),ref(ind));
        r(i) = cr(1,2);
    end

    r(r<th) = min(r);
    I = find(r>=th);

    peaks = zeros(1,N);
    for i = I
        ind = max(1,i-wlen):min(N,i+wlen);
        if (max(r(ind))==r(i))
            peaks(i) = 1;
        else
            peaks(i) = 0;
        end
    end

    [s, W, ~] = PiCA(dat,peaks);
    ref = s(1,:);
    refslice = W*refslice;
end