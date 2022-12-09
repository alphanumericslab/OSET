function [A,B] = PiCAMatrixes1(x,peaks)

[phase, ~] = PhaseCalculation(peaks);


% PM time calculation
J = find(peaks);
n1 = diff(J);
prd = round(mean(n1));
wlen = max(n1)-min(n1);

N0 = length(x);
T1 = zeros(length(peaks)-prd-wlen,1);
NN = length(T1);
for t = 1:NN
    df = abs(phase(t) - phase(max(t+prd-wlen,1):min(t+prd+wlen,N0)));
    [~, I] = min(df);
    T1(t) = t + prd + I - wlen -1;
end
T1 = max(T1,1);
T1 = min(T1,N0);

T0 = 1:NN;

A = x(:,T0)*x(:,T1)';
B = x(:,T0)*x(:,T0)';
A = (A+A')/2;
B = (B+B')/2;
