function [T0,T1] = CalculateTimeLags(peaks,phase)


% PM time calculation
J = find(peaks);
n1 = diff(J);
prd = round(mean(n1));
wlen = max(n1)-min(n1);

T1 = zeros(length(peaks)-prd-wlen,1);
NN = length(T1);
for t = 1:NN
    df = abs(phase(t) - phase(t+prd-wlen:t+prd+wlen));
    [~, I] = min(df);
    T1(t) = t + prd + I - wlen -1;
end

T0 = 1:NN;