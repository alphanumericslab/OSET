function r = InstCorr (x,tau)
% Instantaneous signal correlation calculator
% Reza Sameni, 2010

x = x(:);
N = length(x);
M = ceil(max(abs(tau))/2);
xx = [zeros(1,M) x' zeros(1,M)];

r = zeros(N,length(tau));
t = 1:N;
for i = 1:length(tau)
    tt1 = t - round(tau(i)/2);
    tt2 = t + round(tau(i)/2);
    r(:,i) = xx(tt1+M).*xx(tt2+M);
end