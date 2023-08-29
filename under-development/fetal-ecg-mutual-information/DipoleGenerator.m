function [DIP teta]= DipoleGenerator(N,fs,f,alphai,bi,tetai,teta0);

w = 2*pi*f;
dt = 1/fs;

X = zeros(1,N);
Y = zeros(1,N);
Z = zeros(1,N);
teta = zeros(1,N);
teta(1:end) = teta0:w*dt:(w*dt)*(N-1)+teta0;

for j = 1:length(alphai.x),
    dtetaix = rem(teta - tetai.x(j) + pi,2*pi)-pi;
    X = X + alphai.x(j) .* exp(-dtetaix .^2 ./ (2*bi.x(j) .^ 2));   % x state variable
end

for j = 1:length(alphai.y),
    dtetaiy = rem(teta - tetai.y(j) + pi,2*pi)-pi;
    Y = Y + alphai.y(j) .* exp(-dtetaiy .^2 ./ (2*bi.y(j) .^ 2));   % x state variable
end

for j = 1:length(alphai.z),
    dtetaiz = rem(teta - tetai.z(j) + pi,2*pi)-pi;
    Z = Z + alphai.z(j) .* exp(-dtetaiz .^2 ./ (2*bi.z(j) .^ 2));   % x state variable
end

DIP.x = X;
DIP.y = Y;
DIP.z = Z;
teta = rem(teta+ pi,2*pi)-pi;
