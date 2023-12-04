function peaks1 = LocalPeakDetection(x, half_win_len, number)

N = length(x);
peaks = zeros(N,1);
peaks1 = zeros(N,1);
half_win_len = round(half_win_len);

x = abs(x);
for j = 1:N
    index = max(j-half_win_len,1) : min(j+half_win_len,N);
    % if (x(j)>=0 && max(x(index))==x(j)) || (x(j)<0 && min(x(index))==x(j)) % find all local max or local min
    if max(x(index))==x(j)% find all local maxima
        peaks(j) = 1;
    end
end

% sort the signal local peaks in order or amplitude
[Y, ~] = sort(x(peaks == 1),'descend');
for i = 1 : number % search for up to 'number' of peaks
    ind = find(x==Y(i)); % find all the local peaks that match that ampitude
    for j = 1:length(ind)
        if peaks(ind(j)) == 1 
            peaks1(ind(j)) = 1;
        end
    end
end