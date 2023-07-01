function [T0,T1] = SynchPhaseTimes(peaks,phase)

L3 = length(peaks);
J = find(peaks);
n1 = diff(J);
prd = round(mean(n1));
wlen = max(n1)-min(n1);

T1 = zeros(1,L3-prd-wlen);
NN = length(T1);
for t = 1:NN,
    df = abs(phase(t) - phase(max(t+prd-wlen,1):min(t+prd+wlen,L3)));
    [Y,I] = min(df);
    T1(t) = t + prd + I - wlen -1;
end
T1 = max(T1,1);
T1 = min(T1,NN);
T0 = 1:NN;
