function [X teta]= ECGGenerator(N,fs,f,alphai,bi,tetai,teta0);

w = 2*pi*f;
dt = 1/fs;

X = zeros(1,N);
teta = zeros(1,N);
teta(1:end) = teta0:w*dt:(w*dt)*(N-1)+teta0;

for j = 1:length(alphai),
    dtetai = rem(teta - tetai(j) + pi,2*pi)-pi;
    X = X + alphai(j) .* exp(-dtetai .^2 ./ (2*bi(j) .^ 2));   % x state variable
end

teta = rem(teta+ pi,2*pi)-pi;
