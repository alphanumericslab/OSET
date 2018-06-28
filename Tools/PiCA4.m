function [s,peaks] = PiCA4(dat,ref,t0,w,th,wlen,Itr)


N = length(ref);
for k = 1:Itr,
    refframe = ref(max(1,t0-w):min(length(ref),t0+w));
    error = zeros(size(ref));
    for i = 1:N,
        ind = max(1,i-w):min(length(ref),i+w);

        %error(i) = std(refframe(1:length(ind)) - ref(ind),1);

        cr = corrcoef(refframe(1:length(ind)),ref(ind));
        error(i) = cr(1,2);
    end

    r = error;

    %     th2 = th*std(refframe,1);
    %     th2 = th;

    %     r(r>=th2) = max(r);
    r(r<th) = min(r);

    %     I = find(r<th2);
    I = find(r>=th);

    peaks = zeros(size(ref));
    for i = I,
        ind = max(1,i-wlen):min(length(ref),i+wlen);
        %         if (min(r(ind))==r(i)),
        if (max(r(ind))==r(i)),
            peaks(i) = 1;
        else
            peaks(i) = 0;
        end
    end

    [s,W,A] = PiCA(dat,peaks);
    ref = s(1,:);
end