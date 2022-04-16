function peaks1 = LocalPeakDetection(x,w,number)

N = length(x);
peaks = zeros(N,1);
peaks1 = zeros(N,1);
win = round(w);

x = abs(x);
for j = 1:N,
    index = max(j-win,1):min(j+win,N);
    if((x(j)>=0 & max(x(index))==x(j)) | (x(j)<0 & min(x(index))==x(j)))
        peaks(j) = 1;
    end
end

J = find(peaks);
[Y,I] = sort(x(J),'descend');
for i = 1:number,
    ind = find(x==Y(i));
    for j = 1:length(ind),
        if(peaks(ind(j))==1)
            peaks1(ind(j)) = 1;
        end
    end
end